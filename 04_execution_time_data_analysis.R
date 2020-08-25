library(data.table)
library(ggplot2)

path_to_results = "./time_results/"
path_to_datasets = "./datasets"

datasets = as.data.table(readRDS("datasets.RDS"))
datasets[, type := ifelse(grepl("TMT", path), "TMT", "LF")]
datasets_summary = datasets[, .(total_size_disk = sum(file_size_on_disk),
                                total_size_ram = sum(size_in_ram),
                                total_number_rows = sum(number_rows),
                                total_number_cols = sum(number_columns)),
                            by = c("path", "type", "tool")]
datasets_summary[, path := gsub(path_to_datasets, "", path)]
datasets_summary[, path := gsub("MaxQuant|PD|Progenesis|Skyline|OpenMS|DIAUmpire|Spectronaut|OpenSWATH|SpectroMine",
                                "", path)]
datasets_summary[, mode := ifelse(grepl("DIA", path), "DIA", "DDA")]
datasets_summary[, path := gsub("DDA|DIA|TMT", "", path)]
datasets_summary[, path := gsub("-", "", path)]

time_files = list.files(path_to_results, full.names = TRUE)
time_files_processed = lapply(time_files, function(file) {
    tryCatch({
        time_results = readRDS(file)
        colnames(time_results) = c("version", "time")
        mode = ifelse(grepl("DIA", file), "DIA", "DDA")
        tmt = ifelse(grepl("TMT", file), "TMT", "LF")
        tool = stringr::str_extract(file, "MaxQuant|PD|Progenesis|Skyline|OpenMS|DIAUmpire|Spectronaut|OpenSWATH|SpectroMine")
        time_results$mode = mode
        time_results$type = tmt
        time_results$tool = tool
        time_results$path = file
        time_results
    }, error = function(e) data.table())
})
times = rbindlist(time_files_processed)
times
times[, path := gsub(path_to_results, "", path)]
times[, path := gsub("MaxQuant|PD|Progenesis|Skyline|OpenMS|DIAUmpire|Spectronaut|OpenSWATH|SpectroMine",
                     "", path)]
times[, mode := ifelse(grepl("DIA", path), "DIA", "DDA")]
times[, path := gsub("DDA|DIA|TMT", "", path)]
times[, path := gsub("-", "", path)]
times[, path := gsub(".RDS", "", path)]

times = merge(datasets_summary, times)
times[, time := time / 1e9]
times[, total_size_disk := scale(total_size_disk)]
times[, total_size_ram := scale(total_size_ram)]
times[, total_number_rows := scale(total_number_rows)]
times[, total_number_cols := scale(total_number_cols)]
times[, path := gsub("_", "", path)]
times[, name_short := substr(path, 1, 5)]
ggplot(times, aes(x = name_short, y = time, fill = version)) +
    geom_boxplot() +
    facet_wrap(type~tool, scales = "free") +
    theme_bw()

times_id = copy(times)
times_id = times_id[, id := 1:.N, by = c("path", "type", "tool", "mode", "version")]
times_id = dcast(times_id, path + type + tool + mode + id ~ version, value.var = "time")
times_id[, fold := (v3 - v4) / v3 ]
times_id[, ratio := v4 / v3]
times_id

times_id[tool == "OpenSWATH", ]

ggplot(times_id, aes(x = path, y = ratio)) +
    geom_boxplot() +
    facet_wrap(tool~type, scales = "free") +
    theme_bw()

# TODO: make a plot by dataset

