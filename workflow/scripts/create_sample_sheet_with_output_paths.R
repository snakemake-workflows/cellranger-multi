log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(rlang)
rlang::global_entrace()

library(tidyverse)

samples <- enframe(
    unlist(snakemake@input),
    name = NULL,
    value = "per_sample_counts_path"
  ) |>
  filter(
    str_detect(
      per_sample_counts_path,
      "sample_filtered_feature_bc_matrix/matrix.mtx.gz"
    )
  ) |>
  separate_wider_regex(
    per_sample_counts_path,
    c(
      ".+/cellranger/",
      pool_id = "[^/]+",
      "/outs/per_sample_outs/",
      sample_id = "[^/]+",
      "/.*sample_filtered_feature_bc_matrix/matrix\\.mtx\\.gz"
    ),
    cols_remove = FALSE
  ) |>
  mutate(
    per_sample_counts_path = str_remove(per_sample_counts_path, fixed("sample_filtered_feature_bc_matrix/matrix.mtx.gz")),
    per_sample_counts_path = str_c("../", per_sample_counts_path)
  ) |>
  select(
    sample_id,
    pool_id,
    per_sample_counts_path
  )

write_csv(
  samples,
  file = snakemake@output[["sample_sheet"]],
)
