library(data.table)
source("./02_utility_functions.R")
getStatsSingleVersion = function(df, version) {
    data.table(
        version = version,
        n_features = count_features(df),
        n_proteins = count_proteins(df),
        n_peptides = count_peptides(df),
        n_infinite = count_infinite(df),
        n_missing = count_missing_values(df),
        n_zero = count_zero_values(df),
        n_exactly_zero = count_exactly_zero(df),
        n_rows = count_rows(df),
        n_cols = count_cols(df),
        n_conditions = count_design(df, "Condition"),
        n_bioreps = count_design(df, "BioReplicate"),
        n_runs = count_design(df, "Run"),
        n_tech_replicates = count_design(df, "TechRep"),
        n_tech_rep_mixture = count_design(df, "TechRepMixture"),
        n_fractions = count_fractions(df),
        n_channels = count_design(df, "Channel"))
}
getStats = function(v3, v4) {
    rbindlist(
        list(getStatsSingleVersion(as.data.table(v3), "v3"),
             getStatsSingleVersion(as.data.table(as(v4, "data.frame")), "v4"))
    )
    # list(class(v3), class(v4))
}
# TMT-Plubell2016-OpenMS
om_input_path = "./processed_data/TMT-Plubell2016-OpenMS/20200225_MSstatsTMT_OpenMS_Export.RDS"
om_input = readRDS(om_input_path)
om_v3 = readRDS("./single_outputs/TMT-Plubell2016-OpenMS_v3.RDS")
head(om_v3)
om_v4 = MSstatsConvert::OpenMStoMSstatsTMTFormat(om_input)
head(om_v4)
# TMT-exampledata-SpectroMine
sm_input_path = "./processed_data/TMT-exampledata-SpectroMine/exampledata_input.RDS"
sm_annot_path = "./processed_data/TMT-exampledata-SpectroMine/SpectroMine_annotation.RDS"
sm_input = readRDS(sm_input_path)
sm_annot = readRDS(sm_annot_path)
sm_v3 = readRDS("./single_outputs/TMT-exampledata-SpectroMine_v3.RDS")
sm_v4 = MSstatsConvert::SpectroMinetoMSstatsTMTFormat(sm_input, sm_annot)
# TMT-Controlledmix-MS3-MaxQuant
mq_ev_path = "./processed_data/TMT-Controlledmix-MS3-MaxQuant/evidence.RDS"
mq_pg_path = "./processed_data/TMT-Controlledmix-MS3-MaxQuant/proteinGroups.RDS"
mq_annot_path = "./processed_data/TMT-Controlledmix-MS3-MaxQuant/MaxQuant_annotation.RDS"
mq_ev = readRDS(mq_ev_path)
mq_pg = readRDS(mq_pg_path)
mq_annot = readRDS(mq_annot_path)
mq_v3 = readRDS("./single_outputs/TMT-Controlledmix-MS3-MaxQuant_v3.RDS")
mq_v4 = MSstatsConvert::MaxQtoMSstatsTMTFormat(mq_ev, mq_pg, mq_annot)
# TMT-Controlledmix-SP3-PD
pd_input_path = "./processed_data/TMT-Controlledmix-SP3-PD/spike_input.RDS"
pd_annot_path = "./processed_data/TMT-Controlledmix-SP3-PD/SpikeIn5mix_PD_annotation.RDS"
pd_input = readRDS(pd_input_path)
pd_annot = readRDS(pd_annot_path)
pd_v3 = readRDS("./single_outputs/TMT-Controlledmix-SP3-PD_v3.RDS")
pd_v4 = MSstatsConvert::PDtoMSstatsTMTFormat(pd_input, pd_annot)

ipsc_v3 = readRDS("./single_outputs/TMT_iPSC-PD_v3.RDS")
ipsc_input_path = "./processed_data/TMT-iPSC-PD/ipsc_input.RDS"
ipsc_annot_path = "./processed_data/TMT-iPSC-PD/ipsc_annotation.RDS"
ipsc_input = readRDS(ipsc_input_path)
ipsc_annot = readRDS(ipsc_annot_path)
ipsc_v4 = MSstatsConvert::PDtoMSstatsTMTFormat(ipsc_input, ipsc_annot)

all_comps_tmt = rbindlist(list(
    cbind(tool = "MaxQuant", getStats(mq_v3, mq_v4)),
    cbind(tool = "PD", getStats(pd_v3, pd_v4)),
    cbind(tool = "SpectroMine", getStats(sm_v3, sm_v4)),
    cbind(tool = "OpenMS", getStats(om_v3, om_v4)),
    cbind(tool = "PD_ipsc", getStats(ipsc_v3, ipsc_v4))
), fill = TRUE)
all_comps_tmt

ipsc_v4_check = as.data.table(ipsc_v4)
ipsc_v4_check[, Mixture := as.factor(as.character(Mixture))]
ipsc_v4_check[, TechRepMixture := as.factor(as.character(TechRepMixture))]
ipsc_v4_check[, BioReplicate := as.factor(as.character(BioReplicate))]
ipsc_compare = merge(ipsc_v4_check, as.data.table(ipsc_v3),
                     by = setdiff(colnames(ipsc_v4), "Intensity"),
                     all.x = TRUE)
unique(ipsc_compare[is.na(Intensity.y) & !is.na(Intensity.x), .(PeptideSequence, Charge)])

