log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

rlang::global_entrace()

library(tidyverse)

library_table <- read_tsv(snakemake@input[["sample_sheet"]]) |>
  filter(
    sample == snakemake@wildcards[["sample"]]
  ) |>
  # add NA column if the sample sheet does not list library_type
  full_join(
    tibble(
      library_type = NA_character_
    )
  ) |>
  drop_na(sample) |>
  # provide useful default, if library_type is not specified
  mutate(
    library_type = replace_na(library_type, "Gene Expression")
  ) |>
  add_column(
    fastqs = snakemake@params[["fastqs_dir"]]
  ) |>
  select(
    fastqs,
    sample,
    library_type
  ) |>
  # we might have multiple lanes per sample in the main sample
  # sheet, but only need one entry per sample here
  distinct()

library_table

write_csv(
  library_table,
  file = snakemake@output[["library_csv"]]
)