# import basic packages
import pandas as pd
from os import path
import math
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


wildcard_constraints:
    sample="|".join(sample_sheet["sample"]),
    feature_types="|".join(sample_sheet["feature_types"].replace(" ", "_")),


def get_input_file(wildcards, read_number):
    ft = wildcards.feature_type.replace("_", " ")
    if "lane_number" in sample_sheet.columns:
        return lookup(
            within=sample_sheet,
            query=f"sample == '{wildcards.sample}' & feature_types == '{ft}' & lane_number == '{wildcards.lane_number}'",
            cols=read_number,
        )
    else:
        return lookup(
            within=sample_sheet,
            query=f"sample == '{wildcards.sample}' & feature_types == '{ft}'",
            cols=read_number,
        )


def get_sample_fastqs(wildcards, read_number):
    feature_types = sample_sheet.loc[
        sample_sheet["sample"] == wildcards.sample, "feature_types"
    ].str.replace(" ", "_")
    files = []
    for ft in feature_types:
        # default value to use, if no lane number specified
        lane_numbers = [
            "1",
        ]
        if "lane_number" in sample_sheet.columns:
            lane_numbers = sample_sheet.loc[
                (sample_sheet["sample"] == wildcards.sample)
                & (sample_sheet["feature_types"] == ft.replace("_", " ")),
                "lane_number",
            ]
        files.extend(
            expand(
                "results/input/{sample}_{feature_type}/{sample}_{feature_type}_S1_L00{lane_number}_{read_number}_001.fastq.gz",
                sample=wildcards.sample,
                feature_type=ft,
                lane_number=lane_numbers,
                read_number=read_number,
            )
        )
    return files
