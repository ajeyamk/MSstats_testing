library(data.table)
source("./02_utility_functions.R")
# TMT-Plubell2016-OpenMS
om_input_path = "./processed_data/TMT-Plubell2016-OpenMS/20200225_MSstatsTMT_OpenMS_Export.RDS"
om_input = readRDS(om_input_path)
om_v3 = MSstatsTMT::OpenMStoMSstatsTMTFormat(om_input)
om_v4 = MSstatsConvert::OpenMStoMSstatsTMTFormat(om_input)
# TMT-exampledata-SpectroMine
sm_input_path = "./processed_data/TMT-exampledata-SpectroMine/exampledata_input.RDS"
sm_annot_path = "./processed_data/TMT-exampledata-SpectroMine/SpectroMine_annotation.RDS"
sm_input = readRDS(sm_input_path)
sm_annot = readRDS(sm_annot_path)
sm_v3 = MSstatsTMT::SpectroMinetoMSstatsTMTFormat(sm_input, sm_annot)
sm_v4 = MSstatsConvert::SpectroMinetoMSstatsTMTFormat(sm_input, sm_annot)
# TMT-Controlledmix-MS3-MaxQuant
mq_ev_path = "./processed_data/TMT-Controlledmix_MS3-MaxQuant/evidence.RDS"
mq_pg_path = "./processed_data/TMT-Controlledmix_MS3-MaxQuant/proteinGroups.RDS"
mq_annot_path = "./processed_data/TMT-Controlledmix_MS3-MaxQuant/MaxQuant_annotation.RDS"
mq_ev = readRDS(mq_ev_path)
mq_pg = readRDS(mq_pg_path)
mq_annot = readRDS(mq_annot_path)
mq_v3 = MSstatsTMT::MaxQtoMSstatsTMTFormat(mq_ev, mq_pg, mq_annot)
mq_v4 = MSstatsConvert::MaxQtoMSstatsTMTFormat(mq_ev, mq_pg, mq_annot)
# TMT-Controlledmix-SP3-PD
pd_input_path = "./processed_data/TMT-Controlledmix_SP3-PD/spike_input.RDS"
pd_annot_path = "./processed_data/TMT-Controlledmix_SP3-PD/SpikeIn5mix_PD_annotation.RDS"
pd_input = readRDS(pd_input_path)
pd_annot = readRDS(pd_annot_path)
pd_v3 = MSstatsTMT::PDtoMSstatsTMTFormat(pd_input, pd_annot)
pd_v4 = MSstatsConvert::PDtoMSstatsTMTFormat(pd_input, pd_annot)

ipsc_input_path = "./processed_data/TMT-iPSC-PD/ipsc_input.RDS"
ipsc_annot_path = "./processed_data/TMT-iPSC-PD/ipsc_annotation.RDS"
ipsc_input = readRDS(ipsc_input_path)
ipsc_annot = readRDS(ipsc_annot_path)
ipsc_v3 = MSstatsTMT::PDtoMSstatsTMTFormat(ipsc_input, ipsc_annot)
ipsc_v4 = MSstatsConvert::PDtoMSstatsTMTFormat(ipsc_input, ipsc_annot)

all_comps_tmt = rbindlist(list(
    cbind(tool = "MaxQuant", getStats(mq_v3, mq_v4)),
    cbind(tool = "PD", getStats(pd_v3, pd_v4)),
    cbind(tool = "SpectroMine", getStats(sm_v3, sm_v4)),
    cbind(tool = "OpenMS", getStats(om_v3, om_v4)),
    cbind(tool = "PD_ipsc", getStats(ipsc_v3, ipsc_v4))
), fill = TRUE)
all_comps_tmt
