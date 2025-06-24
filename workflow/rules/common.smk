# import basic packages
import pandas as pd
from os import path
from snakemake.utils import validate


# read sample sheet
sample_sheet = (
    pd.read_csv(config["sample_sheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)


# validate sample sheet and config file
validate(sample_sheet, schema="../../config/schemas/sample_sheet.schema.yaml")
validate(config, schema="../../config/schemas/config.schema.yaml")


def get_input_files(wildcards, read_number):
    if "lane_number" in samples.columns:
        return lookup(
            within=sample_sheet,
            query="sample == '{sample_name}' & lane_number = '{lane_number}'",
            cols=read_number,
        )
    else:
        return lookup(
            within=sample_sheet,
            query="sample == '{sample_name}'",
            cols=read_number,
        )

def get_sample_fastqs(wildcards, read_number):
    sample_row = lookup(
        within=sample_sheet,
        query="sample == '{wildcards.sample}'",
        cols="sample",
    )
    lane_number = 1 if not hasattr(sample_row, "lane_number") else sample_row.lane_number
    return f"results/input/{wildcards.sample}_S1_L00{lane_number}_{read_number}_001.fastq.gz"


def get_all_sample_fastqs(wildcards):
    all_fastqs = []
    for row in sample_sheet.itertuples():
        # provide default lane_number of 1, if no lane_number is specified in sample sheet
        lane_number = 1 if not hasattr(row, "lane_number") else row.lane_number
        all_fastqs.extend(
            [
                f"results/input/{row.sample}_S1_L00{lane_number}_R1_001.fastq.gz",
                f"results/input/{row.sample}_S1_L00{lane_number}_R2_001.fastq.gz",
            ],
        )
    return all_fastqs
