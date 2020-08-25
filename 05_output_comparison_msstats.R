library(MSstats)
library(MSstatsConvert)
library(data.table)
source("./functions.R")
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
        n_fractions = count_fractions(df),
        n_conditions = count_design(df, "Condition"),
        n_bioreps = count_design(df, "BioReplicate"),
        n_runs = count_design(df, "Run"),
        n_tech_replicates = count_design(df, "TechRep"),
        n_tech_rep_mixture = count_design(df, "TechRepMixture"),
        n_channels = count_design(df, "Channel"))
}
getStats = function(v3, v4) {
    rbindlist(
        list(getStatsSingleVersion(as.data.table(v3), "v3"),
             getStatsSingleVersion(as.data.table(as(v4, "data.frame")), "v4"))
    )
    # list(class(v3), class(v4))
}
# DDA-Controlledmix-MaxQuant
maxquant_input_path = "./processed_data/DDA-Controlledmix-MaxQuant/ControlMixture_DDA_MaxQuant_evidence.RDS"
maxquant_pg_path = "./processed_data/DDA-Controlledmix-MaxQuant/ControlMixture_DDA_MaxQuant_proteinGroups.RDS"
maxquant_annot = "./processed_data/DDA-Controlledmix-MaxQuant/ControlMixture_DDA_MaxQuant_annotation.RDS"
mq_ev = readRDS(maxquant_input_path)
mq_pg = readRDS(maxquant_pg_path)
annot = readRDS(maxquant_annot)

maxquant_v3 = MSstats::MaxQtoMSstatsFormat(mq_ev, annot, mq_pg)
maxquant_v4 = MSstatsConvert::MaxQtoMSstatsFormat(mq_ev, annot, mq_pg)

mq_compare = getStats(maxquant_v3, maxquant_v4)
mq_compare

# DDA-Controlledmix-PD
pd_input_path = "./processed_data/DDA-Controlledmix-PD/ControlMixture_DDA_ProteomeDiscoverer_input.RDS"
pd_annot_path = "./processed_data/DDA-Controlledmix-PD/ControlMixture_DDA_ProteomeDiscoverer_annotation.RDS"
pd_input = readRDS(pd_input_path)
pd_annot = readRDS(pd_annot_path)
pd_v3 = MSstats::PDtoMSstatsFormat(pd_input, pd_annot)
pd_v4 = MSstatsConvert::PDtoMSstatsFormat(pd_input, pd_annot)
pd_compare = getStats(pd_v3, as(pd_v4,  "data.frame"))
pd_compare

# DDA-Controlledmix-Progenesis
pg_input_path = "./processed_data/DDA-Controlledmix-Progenesis/ControlMixture_DDA_Progenesis_input.RDS"
pg_annot_path = "./processed_data/DDA-Controlledmix-Progenesis/ControlMixture_DDA_Progenesis_annotation.RDS"
pg_input = readRDS(pg_input_path)
pg_annot = readRDS(pg_annot_path)
pg_v3 = MSstats::ProgenesistoMSstatsFormat(pg_input, pg_annot)
pg_v4 = MSstatsConvert::ProgenesistoMSstatsFormat(pg_input, pg_annot)
warnings()
pg_compare = getStats(pg_v3, pg_v4)
pg_compare
# DDA-Controlledmix-Skyline
sl_input_path = "./processed_data/DDA-Controlledmix-Skyline/ControlMixture_DDA_Skyline_input.RDS"
sl_input = readRDS(sl_input_path)
sl_v3 = MSstats::SkylinetoMSstatsFormat(sl_input)
sl_v4 = MSstatsConvert::SkylinetoMSstatsFormat(sl_input)
sl_compare = getStats(sl_v3, sl_v4)
sl_compare
# DDA-iPRG-OpenMS
om_input_path = "./processed_data/DDA-iPRG-OpenMS/ABRF2015_OpenMS_raw.RDS"
om_input = readRDS(om_input_path)
om_v3 = MSstats::OpenMStoMSstatsFormat(om_input)
om_v4 = MSstatsConvert::OpenMStoMSstatsFormat(om_input)
om_compare = getStats(om_v3, om_v4)
om_compare
# DIA-Navarro2016-DIAUmpire
dia_frag_path = "./processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_FragSummary.RDS"
dia_pept_path = "./processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_PeptideSummary.RDS"
dia_prot_path = "./processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_ProtSummary.RDS"
dia_annot_path = "./processed_data/DIA-Navarro2016-DIAUmpire/Navarro2016_DIA_DIAumpire_input_annotation.RDS"
dia_frag = readRDS(dia_frag_path)
dia_pept = readRDS(dia_pept_path)
dia_prot = readRDS(dia_prot_path)
dia_annot = readRDS(dia_annot_path)
dia_v3 = MSstats::DIAUmpiretoMSstatsFormat(dia_frag, dia_pept, dia_prot, dia_annot)
dia_v4 = MSstatsConvert::DIAUmpiretoMSstatsFormat(dia_frag, dia_pept, dia_prot, dia_annot)
dia_compare = getStats(dia_v3, dia_v4)
dia_compare
# DIA-Navarro2016-OpenSWATH
os_input_path = "./processed_data/DIA-Navarro2016-OpenSWATH/Navarro2016_DIA_OpenSWATH_input.RDS"
os_annot_path = "./processed_data/DIA-Navarro2016-OpenSWATH/Navarro2016_DIA_OpenSWATH_annotation.RDS"
os_input = readRDS(os_input_path)
# os_input$aggr_Peak_Area = as.character(os_input$aggr_Peak_Area)
os_annot = readRDS(os_annot_path)
os_v3 = MSstats::OpenSWATHtoMSstatsFormat(os_input, os_annot)
os_v4 = MSstatsConvert::OpenSWATHtoMSstatsFormat(os_input, os_annot)
os_compare = getStats(os_v3, os_v4)
os_compare
# DIA-Navarro2016-Spectronaut
sn_input_path = "./processed_data/DIA-Navarro2016-Spectronaut/Navarro2016_DIA_Spectronaut_input.RDS"
sn_annot_path = "./processed_data/DIA-Navarro2016-Spectronaut/Navarro2016_DIA_Spectronaut_annotation.RDS"
sn_input = readRDS(sn_input_path)
sn_annot = readRDS(sn_annot_path)
sn_v3 = MSstats::SpectronauttoMSstatsFormat(sn_input, sn_annot)
sn_v4 = MSstatsConvert::SpectronauttoMSstatsFormat(sn_input, sn_annot)
sn_compare = getStats(sn_v3, sn_v4)
sn_compare
# TMT-Plubell2016-OpenMS
# TMT-exampledata-SpectroMine
# TMT-Controlledmix-MS3-MaxQuant
# TMT-Controlledmix-SP3-PD

all_comps = rbindlist(list(
    cbind(tool = "MaxQuant", mq_compare),
    cbind(tool = "PD", pd_compare),
    cbind(tool = "Progenesis", pg_compare),
    cbind(tool = "Skyline", sl_compare),
    cbind(tool = "OpenMS", om_compare),
    cbind(tool = "DIAUmpire", dia_compare),
    cbind(tool = "OpenSWATH", os_compare),
    cbind(tool = "Spectronaut", sn_compare)
), fill = TRUE)
all_comps
