# 1. Functions to calculate statistics
library(data.table)

possible_feature_columns = c("PeptideSequence", "PeptideModifiedSequence",
                             "PrecursorCharge", "Charge", "FragmentIon", "ProductCharge")

count_features = function(df) {
    cols = intersect(colnames(df), possible_feature_columns)
    uniqueN(df[, cols, with = FALSE])
}

count_proteins = function(df) {
    uniqueN(df[, "ProteinName", with = FALSE])
}

count_peptides = function(df) {
    pep_col = intersect(colnames(df), c("PeptideSequence", "PeptideModifiedSequence"))
    uniqueN(df[, pep_col, with = FALSE])
}

get_features_per_protein = function(df) {
    paste_ = function(...) paste(..., sep = "_")
    # df = as.data.table(as(as(df, "MSstatsValidated"), "data.frame"))
    cols = intersect(colnames(df), possible_feature_columns)
    df$feature = do.call("paste_", as.list(df[, cols, with = FALSE]))
    df[, .(n_features = uniqueN(feature)), by = "ProteinName"]
}

get_mean_features_per_protein = function(features_df) {
    mean(features_df$n_features, na.rm = TRUE)
}

count_missing_values = function(df) {
    sum(is.na(df$Intensity), na.rm = TRUE)
}

count_infinite = function(df) {
    sum(!is.finite(df$Intensity))
}

count_zero_values = function(df) {
    sum(abs(df$Intensity) < 1e-6, na.rm = TRUE)
}

count_exactly_zero = function(df) {
    sum(df$Intensity == 0, na.rm = TRUE)
}

get_ram_size = function(df) {
    pryr::object_size(df)
}

get_disk_size = function(file_path) {
    file.info(file_path)$size / 1e6
}

count_rows = function(df) {
    nrow(df)
}

count_cols = function(df) {
    ncol(df)
}

count_measurements = function(df) {
    paste_ = function(...) paste(..., sep = "_")
    cols = intersect(colnames(df), possible_feature_columns)
    df$feature = do.call("paste_", as.list(df[, cols, with = FALSE]))
    count_by = intersect(colnames(df), c("ProteinName", cols, "Run", "Channel"))
    df[, .(n_measurement = .N), by = count_by]
}

count_fractions = function(df) {
    if (is.element("Fraction", colnames(df))) {
        uniqueN(df$Fraction)
    } else {
        NA
    }
}

count_design = function(df, col) {
    if (is.element(col, colnames(df))) {
        uniqueN(df[, col, with = FALSE])
    } else {
        NA
    }
}

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

# files = list.files("./single_outputs", full.names = TRUE)
# statistics = lapply(files, function(path) {
#     result = tryCatch({
#         df = as.data.table(readRDS(path))
#         list(stats = data.table(
#             n_features = count_features(df),
#             n_proteins = count_proteins(df),
#             n_peptides = count_peptides(df),
#             n_missing = count_missing_values(df),
#             n_rows = count_rows(df),
#             n_cols = count_cols(df),
#             n_conditions = count_design(df, "Condition"),
#             n_conditions = count_design(df, "BioReplicate"),
#             n_runs = count_design(df, "Run"),
#             n_tech_replicates = count_design(df, "TechRep"),
#             n_tech_rep_mixture = count_design(df, "TechRepMixture"),
#             n_channels = count_design(df, "Channel")),
#         n_features_per_protein = get_features_per_protein(df),
#         n_measurements = count_measurements(df)
#         )
#
#     }, error = function(e) data.table())
#     return(result)
# })
# saveRDS(statistics, "statistics.RDS")
# statistics
#
# simple_stats = lapply(statistics, function(x) {
#     if(length(x) == 3) {
#         x[[1]]
#     } else {
#         data.table()
#     }
# })
#
# simple_stats_full = lapply(1:length(simple_stats), function(x) {
#     df = simple_stats[[x]]
#     df$path = files[x]
#     df
# })
#
# rbindlist(simple_stats_full, fill = TRUE)
#
# path = files[1]
# df = as.data.table(readRDS(path))
# data.table(
#     n_features = count_features(df),
#     n_proteins = count_proteins(df),
#     n_peptides = count_peptides(df),
#     n_features_per_protein = get_features_per_protein(df),
#     n_missing = count_missing_values(df),
#     n_rows = count_rows(df),
#     n_cols = count_cols(df),
#     n_measurements = count_measurements(df),
#     n_conditions = count_design(df, "Condition"),
#     n_conditions = count_design(df, "BioReplicate"),
#     n_runs = count_design(df, "Run"),
#     n_tech_replicates = count_design(df, "TechRep")
# )
