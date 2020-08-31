library(microbenchmark)
library(parallel)
library(data.table)
source("./02_utility_functions.R")
metadata = readRDS("metadata.RDS")
metadata = metadata[sapply(metadata, function(x) x$type) != "TMT"]

problems = list()
save_single_output = function(metadata) {
    for (dataset in metadata) {
        tryCatch({
            dataset = lapply(dataset, function(x) gsub("./datasets/", "", x))
            if (dataset$tool == "MaxQuant") {
                evidence = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$evidence_path)))
                protein_groups = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$protein_groups)))
                annotation = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation)))
                if (dataset$type == "TMT") {
                    v3 = MSstatsTMT::MaxQtoMSstatsTMTFormat(evidence, protein_groups, annotation)
                    v4 = MSstatsConvert::MaxQtoMSstatsTMTFormat(evidence, protein_groups, annotation)
                } else {
                    v3 = MSstats::MaxQtoMSstatsFormat(evidence, annotation, protein_groups)
                    v4 = MSstatsConvert::MaxQtoMSstatsFormat(evidence, annotation, protein_groups)
                }
            } else if (dataset$tool == "DIAUmpire") {
                fragments = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$fragment_path)))
                peptides = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$peptide_path)))
                proteins = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$protein_path)))
                annotation = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation)))
                v3 = MSstats::DIAUmpiretoMSstatsFormat(fragments, peptides, proteins, annotation)
                v4 = MSstatsConvert::DIAUmpiretoMSstatsFormat(fragments, peptides, proteins, annotation)
            } else if (dataset$tool == "OpenMS") {
                input = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path)))
                if (dataset$type == "TMT") {
                    v3 = MSstatsTMT::OpenMStoMSstatsTMTFormat(input)
                    v4 = MSstatsConvert::OpenMStoMSstatsTMTFormat(input)
                } else {
                    v3 = MSstats::OpenMStoMSstatsFormat(input)
                    v4 = MSstatsConvert::OpenMStoMSstatsFormat(input)
                }
            } else {
                input = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path)))
                annotation = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation_path)))
                input = as.data.frame(input)
                annotation = as.data.frame(annotation)

                if (dataset$tool == "PD") {
                    if (dataset$type == "TMT") {
                        v3 = MSstatsTMT::PDtoMSstatsTMTFormat(input, annotation)
                        v4 = MSstatsConvert::PDtoMSstatsTMTFormat(input, annotation)
                    } else {
                        v3 = MSstats::PDtoMSstatsFormat(input, annotation)
                        v4 = MSstatsConvert::PDtoMSstatsFormat(input, annotation)
                    }
                }

                if (dataset$tool == "OpenSWATH") {
                    v3 = MSstats::OpenSWATHtoMSstatsFormat(input, annotation)
                    v4 = MSstatsConvert::OpenSWATHtoMSstatsFormat(input, annotation)
                }

                if (dataset$tool == "Progenesis") {
                    v3 = MSstats::ProgenesistoMSstatsFormat(input, annotation)
                    v4 = MSstatsConvert::ProgenesistoMSstatsFormat(input, annotation)
                }

                if (dataset$tool == "Skyline") {
                    v3 = MSstats::SkylinetoMSstatsFormat(input, annotation)
                    v4 = MSstatsConvert::SkylinetoMSstatsFormat(input, annotation)
                }

                if (dataset$tool == "SpectroMine") {
                    v3 = MSstatsTMT::SpectroMinetoMSstatsTMTFormat(input, annotation)
                    v4 = MSstatsConvert::SpectroMinetoMSstatsTMTFormat(input, annotation)
                }

                if (dataset$tool == "Spectronaut") {
                    v3 = MSstats::SpectronauttoMSstatsFormat(input, annotation)
                    v4 = MSstatsConvert::SpectronauttoMSstatsFormat(input, annotation)
                }
            }
            saveRDS(v3, file = paste0("./single_outputs/", dataset$folder_path, "_", "v3", ".RDS"))
            saveRDS(v4, file = paste0("./single_outputs/", dataset$folder_path, "_", "v4", ".RDS"))
        }, error = function(e) {
            problems <<- c(problems, dataset$folder_path)
            print(e)
        })
    }
}

save_single_output(metadata)



all_files = list.files("./single_outputs/", full.names = TRUE)
files = data.table(file = all_files)
files[, dataset := gsub("v[0-9].RDS", "", file)]
files[, v := stringr::str_extract(file, "v[0-9]")]
splitted = split(files, files$dataset)
stats_lf = rbindlist(lapply(splitted, function(test) {
    cbind(dataset = unique(gsub("./single_outputs//", "", test$dataset)),
          getStats(readRDS(test[v == "v3", file]),
             readRDS(test[v == "v4", file])))
}))
stats_lf

