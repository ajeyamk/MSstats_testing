library(data.table)
library(microbenchmark)
# 1. Find paths to all datasets
path_to_datasets = "./datasets"
all_paths = list.dirs(path_to_datasets, full.names = FALSE)[-1]
datasets = lapply(all_paths, function(path) {
    x = strsplit(path, "-")[[1]]
    list(folder_path = path,
         type = x[1],
         name = paste0(x[2:(length(x) - 1)], collapse = ""),
         tool = x[length(x)])
})
# 2. Generate metadata
metadata = lapply(datasets, function(dataset) {
    dataset$folder_path = paste(path_to_datasets, dataset$folder_path, sep = "/")
    if (dataset$tool == "MaxQuant") {
        result = list(
            evidence_path = list.files(dataset$folder_path, pattern = "evidence",
                                       full.names = TRUE, include.dirs = TRUE),
            protein_groups = list.files(dataset$folder_path, pattern = "proteinGroups",
                                        full.names = TRUE, include.dirs = TRUE),
            annotation_path = list.files(dataset$folder_path, pattern = "annotation",
                                    full.names = TRUE, include.dirs = TRUE)
        )
    } else if (dataset$tool == "DIAUmpire") {
        result = list(
            fragment_path = list.files(dataset$folder_path, pattern = "FragSummary",
                                       full.names = TRUE, include.dirs = TRUE),
            peptide_path = list.files(dataset$folder_path, pattern = "PeptideSummary",
                                      full.names = TRUE, include.dirs = TRUE),
            protein_path = list.files(dataset$folder_path, pattern = "ProtSummary",
                                      full.names = TRUE, include.dirs = TRUE),
            annotation_path = list.files(dataset$folder_path, pattern = "annotation",
                                    full.names = TRUE, include.dirs = TRUE)
        )
    } else if (dataset$tool == "OpenMS") {
        result = list(
            input_path = list.files(dataset$folder_path, pattern = ".csv",
                                    full.names = TRUE, include.dirs = TRUE)
        )
    } else {
        result = list(
            input_path = list.files(dataset$folder_path, pattern = "input",
                                    full.names = TRUE, include.dirs = TRUE),
            annotation_path = list.files(dataset$folder_path, pattern = "annotation",
                                         full.names = TRUE, include.dirs = TRUE)
        )
    }
    c(dataset, result)
})
saveRDS(metadata, "metadata.RDS")
# 3. Save files to a new folder
lapply(metadata, function(x) {
    system(paste("mkdir", paste0("./processed_data/", gsub(path_to_datasets, "", x$folder_path), collapse = "")))
    for (file_path in intersect(names(x), c("annotation_path", "input_path",
                                            "evidence_path", "protein_groups",
                                            "fragment_path", "peptide_path",
                                            "protein_path", "annotation"))) {
        new_path = gsub("tsv|csv|txt|xls", "RDS", paste0("./processed_data/", gsub(path_to_datasets, "", x[[file_path]])))
        if (!file.exists(new_path)) {
            print(paste0("writing ", new_path, "\n"))
            if (grepl("csv", x[[file_path]])) {
                file_read = read.csv(x[[file_path]])
            } else {
                file_read = read.delim(x[[file_path]])
            }
            saveRDS(file_read, new_path)
        }
    }
})
# 3. Add metadata to datasets
datasets = rbindlist(lapply(metadata, function(x) {
    data.table(
        data = x$name,
        type = x$type,
        tool = x$tool,
        path = x$folder_path,
        file = gsub(".+/", "", unlist(x[!names(x) %in% c("name", "type", "tool", "folder_path")],
                                      use.names = FALSE)),
        file_size_on_disk = sapply(x[!names(x) %in% c("name", "type", "tool", "folder_path")],
                                   function(x) file.info(x)$size/1e6)
    )
}))
datasets$size_in_ram = 0
datasets$number_columns = 0
datasets$number_rows = 0
for (i in 1:nrow(datasets)) {
    new_path = paste0("./processed_data/", gsub(path_to_datasets, "", datasets$path[i]), "/", datasets$file[i], collapse = "")
    new_path = gsub("tsv|csv|xls|txt", "RDS", new_path)
    print(paste("processing file", new_path))
    this_file = readRDS(new_path)
    datasets$size_in_ram[i] = pryr::object_size(this_file)
    datasets$number_columns[i] = ncol(this_file)
    datasets$number_rows[i] = nrow(this_file)
}
saveRDS(datasets, file = "datasets.RDS")
# 4. Tests runs
## 4.1. Files exist
files_exist = lapply(metadata, function(dataset) {
    if (dataset$tool == "MaxQuant") {
        file.exists(paste0("./processed_data", gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS", dataset$evidence_path))))
        file.exists(paste0("./processed_data/", gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS", dataset$protein_groups))))
        file.exists(paste0("./processed_data/", gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation))))
    } else if (dataset$tool == "DIAUmpire") {
        file.exists(paste0("./processed_data/", gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS", dataset$fragment_path))))
        file.exists(paste0("./processed_data/", gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS", dataset$peptide_path))))
        file.exists(paste0("./processed_data/",  gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS",dataset$protein_path))))
        file.exists(paste0("./processed_data/", gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation))))
    } else if (dataset$tool == "OpenMS") {
        file.exists(paste0("./processed_data/", gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path))))
    } else {
        file.exists(paste0("./processed_data/", gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path))))
        file.exists(paste0("./processed_data/", gsub(path_to_datasets, "", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation_path))))
    }
})
all(unlist(files_exist))
## 4.2 converters run
test_run = function(metadata, n_reps = 1) {
    lapply(metadata, function(dataset) {
        dataset = lapply(dataset, function(x) gsub(path_to_datasets, "", x))
        if (dataset$tool == "MaxQuant") {
            evidence = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$evidence_path))))
            protein_groups = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$protein_groups))))
            annotation = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation))))
            if (dataset$type == "TMT") {
                try({microbenchmark(
                    v3 = MSstatsTMT::MaxQtoMSstatsTMTFormat(evidence, protein_groups, annotation),
                    v4 = MSstatsConvert::MaxQtoMSstatsTMTFormat(evidence, protein_groups, annotation),
                    times = n_reps)
                print("MaxQuantTMT - OK")})
            } else {
                try({microbenchmark(
                    v3 = MSstats::MaxQtoMSstatsFormat(evidence, annotation, protein_groups),
                    v4 = MSstatsConvert::MaxQtoMSstatsFormat(evidence, annotation, protein_groups),
                    times = n_reps)
                print("MaxQuant - OK")})
            }
        } else if (dataset$tool == "DIAUmpire") {
            fragments = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$fragment_path))))
            peptides = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$peptide_path))))
            proteins = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$protein_path))))
            annotation = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation))))
            try({microbenchmark(
                v3 = MSstats::DIAUmpiretoMSstatsFormat(fragments, peptides, proteins, annotation),
                v4 = MSstatsConvert::DIAUmpiretoMSstatsFormat(fragments, peptides, proteins, annotation),
                times = n_reps)
            print("DIAUmpire - OK")})
        } else if (dataset$tool == "OpenMS") {
            input = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path))))
            if (dataset$type == "TMT") {
                try({microbenchmark(
                    v3 = MSstatsTMT::OpenMStoMSstatsTMTFormat(input),
                    v4 = MSstatsConvert::OpenMStoMSstatsTMTFormat(input),
                    times = n_reps
                )
                print("OpenMSTMT - OK")})
            } else {
                try({microbenchmark(
                    v3 = MSstats::OpenMStoMSstatsFormat(input),
                    v4 = MSstatsConvert::OpenMStoMSstatsFormat(input),
                    times = n_reps
                )
                print("OpenMS - OK")})
            }
        } else {
            input = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$input_path))))
            annotation = readRDS(gsub(path_to_datasets, "", paste0("./processed_data/", gsub("tsv|csv|xls|txt", "RDS", dataset$annotation_path))))

            if (dataset$tool == "PD") {
                if (dataset$type == "TMT") {
                    try({microbenchmark(
                        v3 = MSstatsTMT::PDtoMSstatsTMTFormat(input, annotation),
                        v4 = MSstatsConvert::PDtoMSstatsTMTFormat(input, annotation),
                        times = n_reps
                    )
                    print("PDTMT - OK")})
                } else {
                    try({microbenchmark(
                        v3 = MSstats::PDtoMSstatsFormat(input, annotation),
                        v4 = MSstatsConvert::PDtoMSstatsFormat(input, annotation),
                        times = n_reps
                    )
                    print("PD - OK")})
                }
            }

            if (dataset$tool == "OpenSWATH") {
                try({microbenchmark(
                    v3 = MSstats::OpenSWATHtoMSstatsFormat(input, annotation),
                    v4 = MSstatsConvert::OpenSWATHtoMSstatsFormat(input, annotation),
                    times = n_reps
                )
                print("OpenSWATH - OK")})
            }

            if (dataset$tool == "Progenesis") {
                try({microbenchmark(
                    v3 = MSstats::ProgenesistoMSstatsFormat(input, annotation),
                    v4 = MSstatsConvert::ProgenesistoMSstatsFormat(input, annotation),
                    times = n_reps
                )
                print("Progenesis - OK")})
            }

            if (dataset$tool == "Skyline") {
                try({microbenchmark(
                    v3 = MSstats::SkylinetoMSstatsFormat(input, annotation),
                    v4 = MSstatsConvert::SkylinetoMSstatsFormat(input, annotation),
                    times = n_reps
                )
                print("Skyline - OK")})
            }

            if (dataset$tool == "SpectroMine") {
                try({microbenchmark(
                    v3 = MSstatsTMT::SpectroMinetoMSstatsTMTFormat(input, annotation),
                    v4 = MSstatsConvert::SpectroMinetoMSstatsTMTFormat(input, annotation),
                    times = n_reps
                )
                print("SpectroMine - OK")})
            }

            if (dataset$tool == "Spectronaut") {
                try({microbenchmark(
                    v3 = MSstats::SpectronauttoMSstatsFormat(input, annotation),
                    v4 = MSstatsConvert::SpectronauttoMSstatsFormat(input, annotation),
                    times = n_reps
                )
                print("Spectronaut - OK")})
            }
        }
    })
}
test_run(metadata[1])
