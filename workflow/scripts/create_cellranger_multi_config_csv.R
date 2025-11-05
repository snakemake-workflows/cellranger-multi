log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(rlang)
rlang::global_entrace()

library(tidyverse)
library(cli)

pool_id <- snakemake@wildcards[["pool_id"]]
pool_sheet <- snakemake@input[["pool_sheet"]]

cellranger_fastq_dirs <- enframe(
  snakemake@input[["fq1"]],
  name = NULL,
  value = "filename"
) |>
  separate_wider_regex(
    filename,
    c(
      "results/input/",
      pool_id,
      "_",
      feature_types = "[^/]+",
      "/",
      pool_id,
      "_.+_R1_001.fastq.gz"
    ),
    cols_remove = FALSE
  ) |>
  add_column(
    id = pool_id
  ) |>
  mutate(
    fastqs = normalizePath(dirname(filename)),
    feature_types = str_replace(feature_types, "_", " ")
  ) |>
  select(
    id,
    feature_types,
    fastqs
  ) |>
  distinct()

libraries_table <- read_tsv(
  pool_sheet,
  col_types = cols(.default = col_character())
) |>
  filter(
    id == pool_id
  ) |>
  left_join(
    cellranger_fastq_dirs,
    by = c("id", "feature_types")
  ) |>
  rename(
    fastq_id = id
  ) |>
  bind_rows(
    # ensure the lane_number column exists
    tibble(lane_number = character())
  ) |>
  replace_na(
    list(lane_number = "1")
  ) |>
  # we might have multiple lanes per sample in the main sample
  # sheet, but only need one entry per sample here
  summarize(
    lanes = str_flatten(
      lane_number,
      collapse = "|"
    ),
    .by = any_of(
      c(
        "sample",
        "fastq_id",
        "fastqs",
        "feature_types",
        "physical_library_id",
        "subsample_rate",
        "chemistry"
      )
    )
  ) |>
  select(
    any_of(
      c(
        "sample",
        "fastq_id",
        "fastqs",
        "feature_types",
        "lanes",
        "physical_library_id",
        "subsample_rate",
        "chemistry"
      )
    )
  ) |>
  mutate(
    across(
      everything(),
      ~ replace_na(.x, "")
    )
  )


specified_feature_types <- libraries_table |>
  pull(feature_types)

# Only start writing anything after we have done the libraries parsing.
# Otherwise, we need to debug the sample sheet, anyways.

write_lines(
  "",
  file = snakemake@output[["multi_config_csv"]],
  sep = "", # ensure that the file gets overwritten with every script execution
  append = FALSE
)

parse_and_write_section_if_required <- function(
  feature_types,
  section_heading
) {
  if (any(feature_types %in% specified_feature_types)) {
    section_table <- enframe(
      snakemake@params[["multi_config_csv_sections"]][[section_heading]]
    ) |>
      mutate(
        value = unlist(value)
      ) |>
      filter(value != "") # remove any empty entries, to keep csv succinct
    write_lines(
      str_c("[", section_heading, "]"),
      file = snakemake@output[["multi_config_csv"]],
      append = TRUE
    )
    write_csv(
      section_table,
      file = snakemake@output[["multi_config_csv"]],
      col_names = FALSE,
      append = TRUE
    )
    write_lines(
      "",
      file = snakemake@output[["multi_config_csv"]],
      append = TRUE
    )
  }
}

parse_and_write_section_if_required(
  c(
    "Gene Expression"
  ),
  "gene-expression"
)

parse_and_write_section_if_required(
  c(
    "VDJ",
    "VDJ-T",
    "VDJ-T-GD",
    "VDJ-B"
  ),
  "vdj"
)

parse_and_write_section_if_required(
  c(
    "Antibody Capture",
    "Antigen Capture",
    "CRISPR Guide Capture"
  ),
  "feature"
)

# parsing for antigen-specificity section is different, so we do it without the helper

if ("Antigen Capture" %in% specified_feature_types) {
  # needed for matching up negative control ids with MHC alleles
  feature_reference <- read_csv(
    snakemake@input[["feature_reference"]]
  ) |>
    select(
      any_of(
        c(
          "id",
          "mhc_allele"
        )
      )
    )

  antigen_specificity_table <- enframe(
    snakemake@params[["multi_config_csv_sections"]][["antigen-specificity"]][[
      "control_ids"
    ]],
    name = NULL,
    value = "control_id"
  ) |>
    mutate(
      control_id = as.character(unlist(control_id))
    ) |>
    left_join(
      feature_reference,
      by = join_by(control_id == id)
    ) |>
    bind_rows(
      # ensure the mhc_allele column exists
      tibble(mhc_allele = character())
    ) |>
    mutate(
      # make sure any unavailable value is ""
      mhc_allele = replace_na(mhc_allele, "")
    )

  end_of_line = "\n"
  if (all(antigen_specificity_table |> pull(mhc_allele) == "")) {
    # remove the mhc_allele column, if no sample has an entry in it
    antigen_specificity_table <- antigen_specificity_table |>
      select(-mhc_allele)
    # make sure we get trailing commas after the header and sample lines, for
    # an example multi config csv, see:
    # https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-antigen#bcr
    end_of_line = ",\n"
  }

  write_lines(
    str_c("[antigen-specificity]"),
    file = snakemake@output[["multi_config_csv"]],
    append = TRUE
  )

  write_csv(
    antigen_specificity_table,
    file = snakemake@output[["multi_config_csv"]],
    eol = end_of_line,
    append = TRUE,
    col_names = TRUE
  )

  write_lines(
    "",
    file = snakemake@output[["multi_config_csv"]],
    append = TRUE
  )
}

# parsing for the libraries section is different, so we write without the helper

write_lines(
  "[libraries]",
  file = snakemake@output[["multi_config_csv"]],
  append = TRUE
)

write_csv(
  libraries_table,
  file = snakemake@output[["multi_config_csv"]],
  append = TRUE,
  col_names = TRUE
)

# parsing for the samples section is different, so we write without the helper
if (
  snakemake@params[["multi_config_csv_sections"]][["multiplexing"]][[
    "activate"
  ]]
) {
  multiplexing_sheet <- snakemake@input[["multiplexing"]]

  multiplexing_barcodes <- read_tsv(
    multiplexing_sheet,
    col_types = cols(.default = col_character())
  ) |>
    filter(
      id == pool_id
    )

  if (n_distinct(multiplexing_barcodes) >= 1) {
    write_lines(
      c("", "[samples]"),
      file = snakemake@output[["multi_config_csv"]],
      append = TRUE
    )

    write_csv(
      multiplexing_barcodes |> select(-id),
      file = snakemake@output[["multi_config_csv"]],
      append = TRUE,
      col_names = TRUE
    )
  }
}
