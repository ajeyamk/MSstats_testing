# 1. Datasets
library(MSstats)
library(MSstatsConvert)
library(data.table)
source("./08_data_process_modules.R")
# DDA-Controlledmix-MaxQuant
maxquant_input_path = "./processed_data/DDA-Controlledmix-MaxQuant/ControlMixture_DDA_MaxQuant_evidence.RDS"
maxquant_pg_path = "./processed_data/DDA-Controlledmix-MaxQuant/ControlMixture_DDA_MaxQuant_proteinGroups.RDS"
maxquant_annot = "./processed_data/DDA-Controlledmix-MaxQuant/ControlMixture_DDA_MaxQuant_annotation.RDS"
mq_ev = readRDS(maxquant_input_path)
mq_pg = readRDS(maxquant_pg_path)
annot = readRDS(maxquant_annot)
maxquant_v4 = MSstatsConvert::MaxQtoMSstatsFormat(mq_ev, annot, mq_pg)
# DDA-Controlledmix-PD
pd_input_path = "./processed_data/DDA-Controlledmix-PD/ControlMixture_DDA_ProteomeDiscoverer_input.RDS"
pd_annot_path = "./processed_data/DDA-Controlledmix-PD/ControlMixture_DDA_ProteomeDiscoverer_annotation.RDS"
pd_input = readRDS(pd_input_path)
pd_annot = readRDS(pd_annot_path)
pd_v4 = MSstatsConvert::PDtoMSstatsFormat(pd_input, pd_annot)
# DDA-Controlledmix-Progenesis
pg_input_path = "./processed_data/DDA-Controlledmix-Progenesis/ControlMixture_DDA_Progenesis_input.RDS"
pg_annot_path = "./processed_data/DDA-Controlledmix-Progenesis/ControlMixture_DDA_Progenesis_annotation.RDS"
pg_input = readRDS(pg_input_path)
pg_annot = readRDS(pg_annot_path)
pg_v4 = MSstatsConvert::ProgenesistoMSstatsFormat(pg_input, pg_annot)
# DDA-Controlledmix-Skyline
sl_input_path = "./processed_data/DDA-Controlledmix-Skyline/ControlMixture_DDA_Skyline_input.RDS"
sl_input = readRDS(sl_input_path)
sl_v4 = MSstatsConvert::SkylinetoMSstatsFormat(sl_input)
# DDA-iPRG-OpenMS
om_input_path = "./processed_data/DDA-iPRG-OpenMS/ABRF2015_OpenMS_raw.RDS"
om_input = readRDS(om_input_path)
om_v4 = MSstatsConvert::OpenMStoMSstatsFormat(om_input)
# DIA-Navarro2016-DIAUmpire
dia_frag_path = "./processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_FragSummary.RDS"
dia_pept_path = "./processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_PeptideSummary.RDS"
dia_prot_path = "./processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_ProtSummary.RDS"
dia_annot_path = "./processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_annotation.RDS"
dia_frag = readRDS(dia_frag_path)
dia_pept = readRDS(dia_pept_path)
dia_prot = readRDS(dia_prot_path)
dia_annot = readRDS(dia_annot_path)
dia_v4 = MSstatsConvert::DIAUmpiretoMSstatsFormat(dia_frag, dia_pept, dia_prot, dia_annot)
# DIA-Navarro2016-OpenSWATH
os_input_path = "./processed_data/DIA-Navarro2016-OpenSWATH/Navarro2016_DIA_OpenSWATH_input.RDS"
os_annot_path = "./processed_data/DIA-Navarro2016-OpenSWATH/Navarro2016_DIA_OpenSWATH_annotation.RDS"
os_input = readRDS(os_input_path)
os_annot = readRDS(os_annot_path)
os_v4 = MSstatsConvert::OpenSWATHtoMSstatsFormat(os_input, os_annot)
# DIA-Navarro2016-Spectronaut
sn_input_path = "./processed_data/DIA-Navarro2016-Spectronaut/Navarro2016_DIA_Spectronaut_input.RDS"
sn_annot_path = "./processed_data/DIA-Navarro2016-Spectronaut/Navarro2016_DIA_Spectronaut_annotation.RDS"
sn_input = readRDS(sn_input_path)
sn_annot = readRDS(sn_annot_path)
sn_v4 = MSstatsConvert::SpectronauttoMSstatsFormat(sn_input, sn_annot)

# 2. Preparation + fractions + balanced design
## v3
dp_mq = .checkDataOld(maxquant_v4, matrix(""), "log.log")
dp_mq = .reformatOld(dp_mq, matrix(c("")), "log.log")
logTrans = 2
dp_mq = .intensitiesOld(dp_mq, matrix(c("")), "log.log")
multirun = MSstats:::.countMultiRun(dp_mq)
dp_mq = .fractionsOld(dp_mq, multirun, matrix(""), "log.log")
dp_mq = .balancedOld(dp_mq, multirun, TRUE, matrix(""), "log.log")
## v4
dp_mq_v4 = MSstatsdev:::.checkDataValidity(maxquant_v4)
dp_mq_v4 = MSstatsdev:::.updateColumnsForProcessing(dp_mq_v4)
MSstatsdev:::.preProcessIntensities(dp_mq_v4, 2)
library(data.table)
colnames(dp_mq)
colnames(dp_mq_v4)

dp_mq_v4[, GROUP_ORIGINAL := as.factor(as.character(GROUP_ORIGINAL))]
dp_mq_v4[, SUBJECT_ORIGINAL := as.factor(as.character(SUBJECT_ORIGINAL))]

comp = merge(dp_mq_v4, as.data.table(dp_mq),
             by = setdiff(colnames(dp_mq), c("INTENSITY", "ABUNDANCE")),
             all.x = TRUE, all.y = TRUE)
head(comp)
comp[ABUNDANCE.x != ABUNDANCE.y]
# 3. Pre-norm + Normalization
dp_mq = .prenormOld(dp_mq, matrix(""), "log.log")
dp_mq = .normalizeOld(dp_mq, "EQUALIZEMEDIANS", NULL, matrix(""), "log.log")
dp_mq_v4 = MSstatsdev:::.makeFactorColumns(dp_mq_v4)
dp_mq_v4 = MSstatsdev:::.normalize(dp_mq_v4, "EQUALIZEMEDIANS")
comp2 = merge(dp_mq_v4, as.data.table(dp_mq),
              by = setdiff(colnames(dp_mq), c("INTENSITY", "ABUNDANCE")),
              all.x = TRUE, all.y = TRUE)
head(comp2)
# 4. Feature selection (top N)
library(dplyr)
.selectFeaturesOld(dp_mq, matrix(""), "topN", 3, "log.log")

# 5. Summarization


dp_output_full = MSstats::dataProcess(maxquant_v4, MBimpute = F)
dp_output_full_v4 = MSstatsdev::dataProcess(maxquant_v4, MBimpute = F)


# dp_output_full = MSstats::dataProcess(maxquant_v4, summaryMethod = "linear")
# dp_output_full_v4 = MSstatsdev::dataProcess(maxquant_v4, summaryMethod = "linear")

dp_output_full = MSstats::dataProcess(maxquant_v4)
dp_output_full_v4 = MSstatsdev::dataProcess(maxquant_v4)


# MSstatsdev::dataProcess(maxquant_v4)

summarized_v3 = dp_output_full$RunlevelData
summarized_v4 = dp_output_full_v4$RunlevelData

colnames(summarized_v3)
colnames(summarized_v4)

summarized_v4 = as.data.table(as.data.table(summarized_v4))
summarized_v4[, GROUP_ORIGINAL := as.factor(as.character(GROUP_ORIGINAL))]
summarized_v4[, SUBJECT_ORIGINAL := as.factor(as.character(SUBJECT_ORIGINAL))]
compare_summarized = merge(
    as.data.table(summarized_v3),
    as.data.table(summarized_v4),
    # by = c("RUN", "Protein"),
    by = c("RUN", "Protein", "originalRUN", "GROUP", "GROUP_ORIGINAL",
           "SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT"),
    # by = intersect(colnames(summarized_v3), colnames(summarized_v4)),
    all.x = T, all.y = T
)
# setdiff(colnames(summarized_v3), colnames(summarized_v4))
# setdiff(colnames(summarized_v4), colnames(summarized_v3))
# colnames(summarized_v4)
# colnames(compare_summarized)
head(compare_summarized)

compare_summarized[LogIntensities.x != LogIntensities.y]
compare_summarized[NumImputedFeature.x != NumImputedFeature.y]

compare_summarized$diff = compare_summarized$NumImputedFeature.x - compare_summarized$NumImputedFeature.y

compare_summarized[MissingPercentage.x != MissingPercentage.y]
compare_summarized[NumMeasuredFeature.x != NumMeasuredFeature.y]


compare_summarized[NumImputedFeature.x == 0 & (LogIntensities.x != LogIntensities.y)]

head(compare_summarized)

table(compare_summarized$diff)

processed_v3 = as.data.table(dp_output_full$ProcessedData)
processed_v4 = as.data.table(dp_output_full_v4$ProcessedData)

colnames(processed_v3)
processed_v4[, GROUP_ORIGINAL := as.factor(as.character(GROUP_ORIGINAL))]
processed_v4[, SUBJECT_ORIGINAL := as.factor(as.character(SUBJECT_ORIGINAL))]

compare_processed = merge(processed_v3, processed_v4,
                          by = setdiff(colnames(processed_v3),
                                       c("INTENSITY", "ABUNDANCE", "censored")),
                          all.x = T, all.y = T)
compare_processed
compare_processed[ABUNDANCE.x != ABUNDANCE.y]

problematic = processed_v4[PROTEIN == "P02929"]

MSstatsdev:::.runTukey(problematic, FALSE, "NA", FALSE, FALSE, 2)

processed_v3[PROTEIN == "P02929"]
as.data.table(summarized_v3)[Protein == "P02929"]



