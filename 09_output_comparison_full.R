library(microbenchmark)
library(parallel)
library(data.table)
source("./02_utility_functions.R")
metadata = readRDS("metadata.RDS")
# metadata = metadata[sapply(metadata, function(x) x$type) != "TMT"]

# Additional options to check:
# Remove 1 feature proteins
# Do not remove few measurements
get_comparison_df = function(path) {
    all_files = list.files(path, full.names = TRUE)
    files = data.table(file = all_files)
    files[, dataset := gsub("v[0-9].RDS", "", file)]
    files[, v := stringr::str_extract(file, "v[0-9]")]
    splitted = split(files, files$dataset)
    stats_lf = rbindlist(lapply(splitted, function(test) {
        cbind(dataset = unique(gsub(paste0(path, "/"), "", test$dataset)),
              getStats(readRDS(test[v == "v3", file]),
                       readRDS(test[v == "v4", file])))
    }))
    stats_lf
}


problems = list()
save_single_output = function(metadata, output_path = "./single_outputs/", remove_single_feature = FALSE, remove_few = TRUE) {
    remove_few_lf = ifelse(remove_few, "remove", "keep")
    for (dataset in metadata) {
        tryCatch({
            dataset = lapply(dataset, function(x) gsub("./datasets/", "", x))
            if (dataset$tool == "MaxQuant") {
                evidence = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$evidence_path)))
                protein_groups = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$protein_groups)))
                annotation = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation)))
                if (dataset$type == "TMT") {
                    v3 = MSstatsTMT::MaxQtoMSstatsTMTFormat(evidence, protein_groups, annotation,
                                                            rmPSM_withfewMea_withinRun = remove_few,
                                                            rmProtein_with1Feature = remove_single_feature)
                    v4 = MSstatsConvert::MaxQtoMSstatsTMTFormat(evidence, protein_groups, annotation,
                                                                rmPSM_withfewMea_withinRun = remove_few,
                                                                rmProtein_with1Feature = remove_single_feature)
                } else {
                    v3 = MSstats::MaxQtoMSstatsFormat(evidence, annotation, protein_groups,
                                                      removeProtein_with1Peptide = remove_single_feature,
                                                      fewMeasurements = remove_few_lf)
                    v4 = MSstatsConvert::MaxQtoMSstatsFormat(evidence, annotation, protein_groups,
                                                             removeProtein_with1Peptide = remove_single_feature,
                                                             fewMeasurements = remove_few_lf)
                }
            } else if (dataset$tool == "DIAUmpire") {
                fragments = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$fragment_path)))
                peptides = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$peptide_path)))
                proteins = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$protein_path)))
                annotation = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation)))
                v3 = MSstats::DIAUmpiretoMSstatsFormat(fragments, peptides, proteins, annotation,
                                                       removeProtein_with1Feature = remove_single_feature,
                                                       fewMeasurements = remove_few_lf)
                v4 = MSstatsConvert::DIAUmpiretoMSstatsFormat(fragments, peptides, proteins, annotation,
                                                              removeProtein_with1Feature = remove_single_feature,
                                                              fewMeasurements = remove_few_lf)
            } else if (dataset$tool == "OpenMS") {
                input = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path)))
                if (dataset$type == "TMT") {
                    v3 = MSstatsTMT::OpenMStoMSstatsTMTFormat(input,
                                                              rmPSM_withfewMea_withinRun = remove_few,
                                                              rmProtein_with1Feature = remove_single_feature)
                    v4 = MSstatsConvert::OpenMStoMSstatsTMTFormat(input,
                                                                  rmPSM_withfewMea_withinRun = remove_few,
                                                                  rmProtein_with1Feature = remove_single_feature)
                } else {
                    v3 = MSstats::OpenMStoMSstatsFormat(input,
                                                        removeProtein_with1Peptide = remove_single_feature,
                                                        fewMeasurements = remove_few_lf)
                    v4 = MSstatsConvert::OpenMStoMSstatsFormat(input,
                                                               removeProtein_with1Peptide = remove_single_feature,
                                                               fewMeasurements = remove_few_lf)
                }
            } else {
                input = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path)))
                annotation = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation_path)))
                input = as.data.frame(input)
                annotation = as.data.frame(annotation)

                if (dataset$tool == "PD") {
                    if (dataset$type == "TMT") {
                        v3 = MSstatsTMT::PDtoMSstatsTMTFormat(input, annotation,
                                                              rmPSM_withfewMea_withinRun = remove_few,
                                                              rmProtein_with1Feature = remove_single_feature)
                        v4 = MSstatsConvert::PDtoMSstatsTMTFormat(input, annotation,
                                                                  rmPSM_withfewMea_withinRun = remove_few,
                                                                  rmProtein_with1Feature = remove_single_feature)
                    } else {
                        v3 = MSstats::PDtoMSstatsFormat(input, annotation,
                                                        removeProtein_with1Peptide = remove_single_feature,
                                                        fewMeasurements = remove_few_lf)
                        v4 = MSstatsConvert::PDtoMSstatsFormat(input, annotation,
                                                               removeProtein_with1Peptide = remove_single_feature,
                                                               fewMeasurements = remove_few_lf)
                    }
                }

                if (dataset$tool == "OpenSWATH") {
                    v3 = MSstats::OpenSWATHtoMSstatsFormat(input, annotation,
                                                           removeProtein_with1Feature = remove_single_feature,
                                                           fewMeasurements = remove_few_lf)
                    v4 = MSstatsConvert::OpenSWATHtoMSstatsFormat(input, annotation,
                                                                  removeProtein_with1Feature = remove_single_feature,
                                                                  fewMeasurements = remove_few_lf)
                }

                if (dataset$tool == "Progenesis") {
                    v3 = MSstats::ProgenesistoMSstatsFormat(input, annotation,
                                                            removeProtein_with1Peptide = remove_single_feature,
                                                            fewMeasurements = remove_few_lf)
                    v4 = MSstatsConvert::ProgenesistoMSstatsFormat(input, annotation,
                                                                   removeProtein_with1Peptide = remove_single_feature,
                                                                   fewMeasurements = remove_few_lf)
                }

                if (dataset$tool == "Skyline") {
                    v3 = MSstats::SkylinetoMSstatsFormat(input, annotation,
                                                         removeProtein_with1Feature = remove_single_feature,
                                                         fewMeasurements = remove_few_lf)
                    v4 = MSstatsConvert::SkylinetoMSstatsFormat(input, annotation,
                                                                removeProtein_with1Feature = remove_single_feature,
                                                                fewMeasurements = remove_few_lf)
                }

                if (dataset$tool == "SpectroMine") {
                    v3 = MSstatsTMT::SpectroMinetoMSstatsTMTFormat(input, annotation,
                                                                   rmPSM_withfewMea_withinRun = remove_few,
                                                                   rmProtein_with1Feature = remove_single_feature)
                    v4 = MSstatsConvert::SpectroMinetoMSstatsTMTFormat(input, annotation,
                                                                       rmPSM_withfewMea_withinRun = remove_few,
                                                                       rmProtein_with1Feature = remove_single_feature)
                }

                if (dataset$tool == "Spectronaut") {
                    v3 = MSstats::SpectronauttoMSstatsFormat(input, annotation,
                                                             removeProtein_with1Feature = remove_single_feature,
                                                             fewMeasurements = remove_few_lf)
                    v4 = MSstatsConvert::SpectronauttoMSstatsFormat(input, annotation,
                                                                    removeProtein_with1Feature = remove_single_feature,
                                                                    fewMeasurements = remove_few_lf)
                }
            }
            saveRDS(v3, file = paste0(output_path, dataset$folder_path, "_", "v3", ".RDS"))
            saveRDS(v4, file = paste0(output_path, dataset$folder_path, "_", "v4", ".RDS"))
        }, error = function(e) {
            problems <<- c(problems, dataset$folder_path)
            print(e)
        })
    }
}

# metadata = metadata[sapply(metadata, function(x) x$tool) == "Skyline"]
save_single_output(metadata, "./single_outputs/")
save_single_output(metadata, output_path = "./single_outputs_no_few/", remove_single_feature = FALSE, remove_few = FALSE)
save_single_output(metadata, output_path = "./single_outputs_no_singles/", remove_single_feature = TRUE, remove_few = TRUE)

compare_default = get_comparison_df("./single_outputs/")
compare_no_few = get_comparison_df("./single_outputs_no_few/")
compare_no_singles = get_comparison_df("./single_outputs_no_singles/")
