library(microbenchmark)
library(parallel)
library(data.table)
metadata = readRDS("metadata.RDS")
measure_time = function(metadata, output_path, n_reps = 10, n_threads = 1,
                        path_to_datasets = "./datasets") {
    data.table::setDTthreads(n_threads)
    results = list()
    for (dataset in metadata) {
        dataset = lapply(dataset, function(x) gsub(path_to_datasets, "", x))
        if (dataset$tool == "MaxQuant") {
            evidence = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$evidence_path)))
            protein_groups = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$protein_groups)))
            annotation = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation)))
            if (dataset$type == "TMT") {
                result = tryCatch(as.data.frame(microbenchmark(
                    v3 = MSstatsTMT::MaxQtoMSstatsTMTFormat(evidence, protein_groups, annotation),
                    v4 = MSstatsConvert::MaxQtoMSstatsTMTFormat(evidence, protein_groups, annotation),
                    times = n_reps)), error = function(e) data.frame())
            } else {
                result = tryCatch(as.data.frame(microbenchmark(
                    v3 = MSstats::MaxQtoMSstatsFormat(evidence, annotation, protein_groups),
                    v4 = MSstatsConvert::MaxQtoMSstatsFormat(evidence, annotation, protein_groups),
                    times = n_reps)), error = function(e) data.frame())
            }
        } else if (dataset$tool == "DIAUmpire") {
            fragments = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$fragment_path)))
            peptides = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$peptide_path)))
            proteins = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$protein_path)))
            annotation = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation)))
            result = tryCatch(as.data.frame(microbenchmark(
                v3 = MSstats::DIAUmpiretoMSstatsFormat(fragments, peptides, proteins, annotation),
                v4 = MSstatsConvert::DIAUmpiretoMSstatsFormat(fragments, peptides, proteins, annotation),
                times = n_reps)), error = function(e) data.frame())
        } else if (dataset$tool == "OpenMS") {
            input = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path)))
            if (dataset$type == "TMT") {
                result = tryCatch(as.data.frame(microbenchmark(
                    v3 = MSstatsTMT::OpenMStoMSstatsTMTFormat(input),
                    v4 = MSstatsConvert::OpenMStoMSstatsTMTFormat(input),
                    times = n_reps
                )), error = function(e) data.frame())
            } else {
                result = tryCatch(as.data.frame(microbenchmark(
                    v3 = MSstats::OpenMStoMSstatsFormat(input),
                    v4 = MSstatsConvert::OpenMStoMSstatsFormat(input),
                    times = n_reps
                )), error = function(e) data.frame())
            }
        } else {
            input = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path)))
            annotation = readRDS(paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation_path)))

            if (dataset$tool == "PD") {
                if (dataset$type == "TMT") {
                    result = tryCatch(as.data.frame(microbenchmark(
                        v3 = MSstatsTMT::PDtoMSstatsTMTFormat(input, annotation),
                        v4 = MSstatsConvert::PDtoMSstatsTMTFormat(input, annotation),
                        times = n_reps
                    )), error = function(e) data.frame())
                } else {
                    result = tryCatch(as.data.frame(microbenchmark(
                        v3 = MSstats::PDtoMSstatsFormat(input, annotation),
                        v4 = MSstatsConvert::PDtoMSstatsFormat(input, annotation),
                        times = n_reps
                    )), error = function(e) data.frame())
                }
            }

            if (dataset$tool == "OpenSWATH") {
                result = tryCatch(as.data.frame(microbenchmark(
                    v3 = MSstats::OpenSWATHtoMSstatsFormat(input, annotation),
                    v4 = MSstatsConvert::OpenSWATHtoMSstatsFormat(input, annotation),
                    times = n_reps
                )), error = function(e) data.frame())
            }

            if (dataset$tool == "Progenesis") {
                result = tryCatch(as.data.frame(microbenchmark(
                    v3 = MSstats::ProgenesistoMSstatsFormat(input, annotation),
                    v4 = MSstatsConvert::ProgenesistoMSstatsFormat(input, annotation),
                    times = n_reps
                )), error = function(e) data.frame())
            }

            if (dataset$tool == "Skyline") {
                result = tryCatch(as.data.frame(microbenchmark(
                    v3 = MSstats::SkylinetoMSstatsFormat(input, annotation),
                    v4 = MSstatsConvert::SkylinetoMSstatsFormat(input, annotation),
                    times = n_reps
                )), error = function(e) data.frame())
            }

            if (dataset$tool == "SpectroMine") {
                result = tryCatch(as.data.frame(microbenchmark(
                    v3 = MSstatsTMT::SpectroMinetoMSstatsTMTFormat(input, annotation),
                    v4 = MSstatsConvert::SpectroMinetoMSstatsTMTFormat(input, annotation),
                    times = n_reps
                )), error = function(e) data.frame())
            }

            if (dataset$tool == "Spectronaut") {
                result = tryCatch(as.data.frame(microbenchmark(
                    v3 = MSstats::SpectronauttoMSstatsFormat(input, annotation),
                    v4 = MSstatsConvert::SpectronauttoMSstatsFormat(input, annotation),
                    times = n_reps
                )), error = function(e) data.frame())
            }
        }
    saveRDS(result, file = paste0(output_path, gsub(path_to_datasets, "", dataset$folder_path), ".RDS"))
    results[[dataset$folder_path]] = result
    }
    results
}

# Test
# measure_time(metadata[1], ".", 1)
