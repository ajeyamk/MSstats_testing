library(data.table)
library(microbenchmark)
raw_unc = readRDS("./processed_data/DDA-iPRG-OpenMS/ABRF2015_OpenMS_raw.RDS")
raw = MSstats::OpenMStoMSstatsFormat(raw_unc)
logTrans = 2
normalization = "equalizeMedians"
nameStandards = NULL
address = ""
fillIncompleteRows = TRUE
featureSubset = "all"
remove_uninformative_feature_outlier = FALSE
n_top_feature = 3
summaryMethod = "TMP"
equalFeatureVar = TRUE
censoredInt = "NA"
cutoffCensored = "minFeature"
MBimpute = TRUE
remove50missing = FALSE
maxQuantileforCensored = 0.999
clusters = NULL
# Old code: setup logging ----
## save process output in each step
allfiles <- list.files()
num <- 0
filenaming <- "msstats"
finalfile <- "msstats.log"
while(is.element(finalfile, allfiles)) {
    num <- num + 1
    finalfile <- paste0(paste(filenaming, num, sep="-"), ".log")
}
session <- sessionInfo()
sink("sessionInfo.txt")
print(session)
sink()

processout <- as.matrix(read.table("sessionInfo.txt", header=TRUE, sep="\t"))
write.table(processout, file=finalfile, row.names=FALSE)
processout <- rbind(processout, as.matrix(c(" "," ","MSstats - dataProcess function"," "), ncol=1))

# New code: setup logging ----
MSstatsdev:::.saveSessionInfo()
getOption("MSstatsLog")("INFO", "MSstats - dataProcess function")
# Old code: check parameters ----
processout = .checkParamsOld(raw, processout, normalization, fillIncompleteRows,
                             logTrans, summaryMethod, cutoffCensored, MBimpute,
                             censoredInt, nameStandards, finalfile)
# New code: check parameters ----

# 1. benchmark: Check parameters ----
microbenchmark::microbenchmark(
    new = MSstatsdev:::.checkDataProcessParams(logTrans, normalization, nameStandards,
                                         address, fillIncompleteRows,
                                         list(remove_uninformative = remove_uninformative_feature_outlier),
                                         list(method = summaryMethod),
                                         list(cutoff = cutoffCensored, symbol = censoredInt,
                                              MB = MBimpute),
                                         n_clusters = clusters),
    old = .checkParamsOld(raw, processout, normalization, fillIncompleteRows,
                    logTrans, summaryMethod, cutoffCensored, MBimpute,
                    censoredInt, nameStandards, finalfile),
    times = 100
)
# 291/564 = 0.51
# Old code: peptides for standards normalization ----
tempPeptide <- unique(raw[, c("PeptideSequence", "PrecursorCharge")])
tempPeptide$PEPTIDE <- paste(tempPeptide$PeptideSequence, tempPeptide$PrecursorCharge, sep="_")
# New code:  peptides for standards normalization ----
input = data.table::as.data.table(as(raw, "data.frame"))
peptides_dict = unique(input[, list(PeptideSequence, PrecursorCharge)])
peptides_dict[, PEPTIDE := paste(PeptideSequence, PrecursorCharge, sep = "_")]
# Very small benchmark: ----
setDTthreads(1)
microbenchmark::microbenchmark(
    new = {
        peptides_dict = unique(input[, list(PeptideSequence, PrecursorCharge)])
        peptides_dict[, PEPTIDE := paste(PeptideSequence, PrecursorCharge, sep = "_")]
    },
    old = {
        tempPeptide <- unique(raw[, c("PeptideSequence", "PrecursorCharge")])
        tempPeptide$PEPTIDE <- paste(tempPeptide$PeptideSequence, tempPeptide$PrecursorCharge, sep="_")
    },
    times = 10
)
# 909 / 1603 = 0.56
# 2. benchmark: Check data validity ----
microbenchmark::microbenchmark(
    new = MSstatsdev:::.checkDataValidity(as.data.table(raw)),
    old = .checkDataOld(raw, processout)
)
# Old code: check data validity ----
further = .checkDataOld(raw, processout, finalfile)
raw.temp = further[[1]]
processout = further[[2]]
# New code: check data validity ----
input = MSstatsdev:::.checkDataValidity(raw)
# 1.80 / 1.82 = 0.99
# 3. benchmark: reformating
microbenchmark(
    new = MSstatsdev:::.updateColumnsForProcessing(input),
    old = .reformatOld(raw.temp, processout, finalfile),
    times = 10
)
# 96 / 487 = 0.19
# Old code: reformat data ----
further = .reformatOld(raw.temp, processout, finalfile)
work = further[[1]]
processout = further[[2]]
# New code: reformat data ----
input = MSstatsdev:::.updateColumnsForProcessing(input)
# 2.88 / 6.57 = 0.43
# 4. microbenchmark: intensities ----
microbenchmark(
    new = MSstatsdev:::.preProcessIntensities(input, logTrans),
    old = .intensitiesOld(work, processout, finalfile),
    times = 100
)
# 3.9 / 4.6 = 0.84
# Old code: modify intensities -----
further = .intensitiesOld(work, processout, finalfile)
work = further[[1]]
processout = further[[2]]
# New code: modify intensities  ----
MSstatsdev:::.preProcessIntensities(input, logTrans)    # rm(raw) # here?
input
# 5. benchmark: check multi run
microbenchmark(
    new = MSstatsdev:::.checkMultiRun(input),
    old = MSstats:::.countMultiRun(work),
    times = 100
)
# 38 / 367 = 0.1
# Old code: check multi run
## Check multi-method or not : multiple run for a replicate
work$RUN <- factor(work$RUN)
checkMultirun <- MSstats:::.countMultiRun(work)
# New code:
check_multi_run = MSstatsdev:::.checkMultiRun(input)
# 6. benchmark: fractions
microbenchmark(
    new = MSstatsdev:::.handleFractions(input, check_multi_run),
    old = .fractionsOld(work, checkMultirun, processout, finalfile),
    times = 10
)
# WTF????
# Old code: fractions:
futher = .fractionsOld(work, checkMultirun, processout, finalfile)
work = further[[1]]
processout = further[[2]]
# New code: fractions
input = MSstatsdev:::.handleFractions(input, check_multi_run)
# 7. microbenchmark: balanced design ----
microbenchmark(
    new = {
        input = MSstatsdev:::.makeBalancedDesign(input, fillIncompleteRows)
        MSstatsdev:::.checkDuplicatedMeasurements(input)
    },
    old = .balancedOld(work, checkMultirun, fillIncompleteRows, processout, finalfile),
    times = 10
)
# 46/198 = 0.23
# Old code: balanced design ----
further = .balancedOld(work, checkMultirun, fillIncompleteRows, processout, finalfile)
work = further[[1]]
processout = further[[2]]
# New code: balanced design ----
input = MSstatsdev:::.makeBalancedDesign(input, fillIncompleteRows)
MSstatsdev:::.checkDuplicatedMeasurements(input)
# Small benchmark - pre-normalization
microbenchmark(
    new = MSstatsdev:::.makeFactorColumns(input),
    old = .prenormOld(work, processout, finalfile),
    times = 10
)
# 168 / 315 = 0.53
# Old code - pre-normalization
further = .prenormOld(work, processout, finalfile)
work = further[[1]]
processout = further[[2]]
# New code - pre-normalization
input = MSstatsdev:::.makeFactorColumns(input)
# # 8. benchmark: normalization
microbenchmark(
    new = MSstatsdev:::.normalize(input, normalization, nameStandards),
    old = .normalizeOld(work, normalization, nameStandards, processout, finalfile),
    times = 10
)
# # WTF???
# profvis::profvis(MSstatsdev:::.normalize(input, normalization, nameStandards))
