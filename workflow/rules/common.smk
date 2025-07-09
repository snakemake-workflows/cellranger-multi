# import basic packages
import pandas as pd
from os import path
from snakemake.utils import validate


# read sample sheet
sample_sheet = (
    pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)
    .set_index("sample", drop=False)
    .sort_index()
)


# validate sample sheet and config file
validate(sample_sheet, schema="../../config/schemas/sample_sheet.schema.yaml")
validate(config, schema="../../config/schemas/config.schema.yaml")


def get_input_file(wildcards, read_number):
    if "lane_number" in sample_sheet.columns:
        return lookup(
            within=sample_sheet,
            query="sample == '{wildcards.sample}' & lane_number == '{wildcards.lane_number}'",
            cols=read_number,
        )
    else:
        return lookup(
            within=sample_sheet,
            query="sample == '{wildcards.sample}'",
            cols=read_number,
        )


def get_sample_fastqs(wildcards, read_number):
    lane_numbers = lookup(
        within=sample_sheet,
        query="sample == '{wildcards.sample}'",
        cols="lane_number",
    )
    lane_number = (
        1 if not lane_numbers else lane_numbers
    )
    return expand(
            "results/input/{sample}_S1_L00{lane_number}_{read_number}_001.fastq.gz",
            sample=wildcards.sample,
            lane_number=lane_numbers,
            read_number=read_number,
    )