# import basic packages
import pandas as pd
from os import path
import math
from snakemake.utils import validate


# read pool sheet
pool_sheet = (
    pd.read_csv(config["pool_sheet"], sep="\t", dtype=str)
    .set_index("id", drop=False)
    .sort_index()
)

# validate pool sheet and config file
validate(pool_sheet, schema="../schemas/pool_sheet.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")

# load and validate multiplexing sheet if activated
if config["multi_config_csv_sections"]["multiplexing"]["activate"]:
    multiplexing_sheet = pd.read_csv(
        config["multi_config_csv_sections"]["multiplexing"]["tsv"],
        sep="\t",
        dtype=str,
    )
    validate(multiplexing_sheet, schema="../schemas/multiplexing_sheet.schema.yaml")


# set global variables

ALL_IDS = pool_sheet["id"].unique()


wildcard_constraints:
    pool_id="|".join(pool_sheet["id"]),
    feature_types="|".join(pool_sheet["feature_types"].replace(" ", "_")),
    vdj_type="|".join(["vdj_b", "vdj_t", "vdj_t_gd"]),


def determine_final_output(wildcards):
    # automatically generated copy of the multi config CSV, which is the only
    # consistent global output of cellranger  (independent of assay types)
    final_output = expand(
        "results/cellranger/{pool_id}/outs/config.csv",
        pool_id=ALL_IDS,
    )

    for pool in ALL_IDS:

        samples = [pool]
        if config["multi_config_csv_sections"]["multiplexing"]["activate"]:
            samples = multiplexing_sheet.loc[
                multiplexing_sheet["id"] == pool, "sample_id"
            ].unique()
            # for non-multiplexed pools in a mixed setup, use pool_id as the
            # sample_id (these pools won't have entries in the multiplexing
            # sheet, while other pools will)
            if len(samples) == 0:
                samples = [pool]

        # request per-sample summaries, which are the only consistent per-sample
        # output of cellranger  (independent of assay types)
        final_output.extend(
            expand(
                [
                    "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/metrics_summary.csv",
                    "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/web_summary.html",
                ],
                pool_id=pool,
                sample_id=samples,
            )
        )

        # handle feature_types defined for this pool of samples
        feature_types = pool_sheet.loc[
            pool_sheet["id"] == pool, "feature_types"
        ].unique()

        # request multiplexing output, if any has been done
        if config["multi_config_csv_sections"]["multiplexing"]["activate"]:
            final_output.extend(
                expand(
                    [
                        "results/cellranger/{pool_id}/outs/multi/count/feature_reference.csv",
                        "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
                        "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/cells_per_tag.json",
                        "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv",
                        "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/tag_calls_summary.csv",
                        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/feature_reference.csv",
                    ],
                    pool_id=pool,
                    sample_id=samples,
                )
            )
            if (
                "Antibody Capture" in feature_types
                or "Antigen Capture" in feature_types
            ):
                final_output.extend(
                    expand(
                        [
                            "results/cellranger/{pool_id}/outs/multi/count/antibody_analysis/aggregate_barcodes.csv",
                        ],
                        pool_id=pool,
                    )
                )
            if "CRISPR Guide Capture" in feature_types:
                final_output.extend(
                    expand(
                        [
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/cells_per_protospacer.json",
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/feature_reference.csv",
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/perturbation_effects_by_feature",
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/perturbation_effects_by_target",
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/perturbation_efficiencies_by_feature.csv",
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/perturbation_efficiencies_by_target.csv",
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/protospacer_calls_per_cell.csv",
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/protospacer_calls_summary.csv",
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/protospacer_umi_thresholds.csv",
                            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/protospacer_umi_thresholds.json",
                        ],
                        pool_id=pool,
                    )
                )

        for ft in feature_types:

            if ft == "Gene Expression":
                final_output.extend(
                    expand(
                        [
                            "results/cellranger/{pool_id}/outs/multi/count/raw_molecule_info.h5",
                            "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz",
                            "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz",
                            "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix/matrix.mtx.gz",
                            "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix.h5",
                            "results/cellranger/{pool_id}/outs/multi/count/unassigned_alignments.bam",
                            "results/cellranger/{pool_id}/outs/multi/count/unassigned_alignments.bam.bai",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_barcodes.csv",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_alignments.bam",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_alignments.bam.bai",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/features.tsv.gz",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix.h5",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_molecule_info.h5",
                        ],
                        pool_id=pool,
                        sample_id=samples,
                    )
                )

            if ft in ["VDJ-B", "VDJ-T", "VDJ-T-GD"]:
                final_output.extend(
                    expand(
                        [
                            "results/cellranger/{pool_id}/outs/vdj_reference/reference.json",
                            "results/cellranger/{pool_id}/outs/vdj_reference/fasta/regions.fa",
                            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig_annotations.bed",
                            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig_annotations.csv",
                            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig_annotations.json",
                            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.bam",
                            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.bam.bai",
                            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.fasta",
                            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.fasta.fai",
                            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.fastq",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/airr_rearrangement.tsv",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/cell_barcodes.json",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/clonotypes.csv",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/concat_ref.bam",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/concat_ref.bam.bai",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/concat_ref.fasta",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/concat_ref.fasta.fai",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus.bam",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus.bam.bai",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus.fasta",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus.fasta.fai",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus_annotations.csv",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/donor_regions.fa",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/filtered_contig_annotations.csv",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/filtered_contig.fasta",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/filtered_contig.fastq",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/vdj_contig_info.pb",
                            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/vloupe.vloupe",
                        ],
                        pool_id=pool,
                        vdj_type=ft.lower().replace("-", "_"),
                        sample_id=samples,
                    )
                )

    return final_output


def get_input_file(wildcards, read_number):
    ft = wildcards.feature_type.replace("_", " ")
    if "lane_number" in pool_sheet.columns:
        return pool_sheet.loc[
            (pool_sheet["id"] == wildcards.pool_id)
            & (pool_sheet["feature_types"] == ft)
            & (pool_sheet["lane_number"] == wildcards.lane_number),
            read_number,
        ].unique()
    else:
        return pool_sheet.loc[
            (pool_sheet["id"] == wildcards.pool_id)
            & (pool_sheet["feature_types"] == ft),
            read_number,
        ].unique()


def get_sample_fastqs(wildcards, read_number):
    feature_types = (
        pool_sheet.loc[pool_sheet["id"] == wildcards.pool_id, "feature_types"]
        .str.replace(" ", "_")
        .unique()
    )
    files = []
    for ft in feature_types:
        # default value to use, if no lane number specified
        lane_numbers = [
            "1",
        ]
        if "lane_number" in pool_sheet.columns:
            lane_numbers = pool_sheet.loc[
                (pool_sheet["id"] == wildcards.pool_id)
                & (pool_sheet["feature_types"] == ft.replace("_", " ")),
                "lane_number",
            ].unique()
        files.extend(
            expand(
                "results/input/{pool_id}_{feature_type}/{pool_id}_S1_L00{lane_number}_{read_number}_001.fastq.gz",
                pool_id=wildcards.pool_id,
                feature_type=ft,
                lane_number=lane_numbers,
                read_number=read_number,
            )
        )
    return files
