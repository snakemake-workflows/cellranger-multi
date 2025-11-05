# Link the raw fastq files into an input directory and
# make sure they are named EXACTLY how Cell Ranger needs
# them. Sigh.
# -----------------------------------------------------
rule follow_pedantic_cell_ranger_naming_scheme:
    input:
        fq1=lambda wc: get_input_file(wc, "read1"),
        fq2=lambda wc: get_input_file(wc, "read2"),
    output:
        fq1="results/input/{pool_id}_{feature_type}/{pool_id}_S1_L00{lane_number}_R1_001.fastq.gz",
        fq2="results/input/{pool_id}_{feature_type}/{pool_id}_S1_L00{lane_number}_R2_001.fastq.gz",
    log:
        "logs/input/{pool_id}_{feature_type}/{pool_id}_{feature_type}_S1_L00{lane_number}_001.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    params:
        fq1=lambda wc, input, output: path.relpath(
            str(input.fq1), start=path.dirname(output.fq1)
        ),
        fq2=lambda wc, input, output: path.relpath(
            str(input.fq2), start=path.dirname(output.fq2)
        ),
    shell:
        "( ln --symbolic {params.fq1} {output.fq1}; "
        "  ln --symbolic {params.fq2} {output.fq2}; "
        ") >{log} 2>&1 "


# Create a multi config CSV file for Cell Ranger.
# -----------------------------------------------------
rule create_cellranger_multi_config_csv:
    input:
        pool_sheet=lookup(within=config, dpath="pool_sheet"),
        fq1=lambda wc: get_sample_fastqs(wc, "R1"),
        fq2=lambda wc: get_sample_fastqs(wc, "R2"),
        feature_reference=lookup(
            within=config,
            dpath="multi_config_csv_sections/feature/reference",
            default=[],
        ),
        multiplexing=branch(
            lookup(
                within=config,
                dpath="multi_config_csv_sections/multiplexing/activate",
            ),
            then=lookup(
                within=config,
                dpath="multi_config_csv_sections/multiplexing/tsv",
                default=[],
            ),
            otherwise=[],
        ),
    output:
        multi_config_csv="results/input/{pool_id}.cell_ranger_multi_config.csv",
    log:
        "logs/input/{pool_id}.cell_ranger_multi_config.log",
    conda:
        "../envs/tidyverse.yaml"
    localrule: True
    params:
        multi_config_csv_sections=lookup(
            within=config, dpath="multi_config_csv_sections"
        ),
    script:
        "../scripts/create_cellranger_multi_config_csv.R"


# Run cellranger multi on one sample.
# -----------------------------------------------------
rule cellranger_multi_run:
    input:
        multi_config_csv="results/input/{pool_id}.cell_ranger_multi_config.csv",
        fq1=lambda wc: get_sample_fastqs(wc, "R1"),
        fq2=lambda wc: get_sample_fastqs(wc, "R2"),
        reference=lookup(
            within=config,
            dpath="multi_config_csv_sections/gene-expression/reference",
            default=[],
        ),
    output:
        "results/cellranger/{pool_id}/outs/config.csv",
        out_dir=directory("results/cellranger/{pool_id}/"),
    log:
        "logs/cellranger/multi/multi_run_{pool_id}.log",
    conda:
        "../envs/cellranger.yaml"
    threads: 16
    resources:
        mem_mb=lambda wc, threads: threads * 4000,
    params:
        mem_gb=lambda wc, resources: math.floor(resources.mem_mb / 1000),
    shell:
        "(rm -rf {output.out_dir}; "
        " cellranger multi "
        "  --id={wildcards.pool_id} "
        "  --output-dir={output.out_dir} "
        "  --csv={input.multi_config_csv} "
        "  --localcores={threads} "
        "  --localmem={params.mem_gb}; "
        ") >{log} 2>&1 "


# Check all of the cellranger output files required by the setup


rule cellranger_multi_files_summaries:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/metrics_summary.csv"
        ),
        report(
            update(
                "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/web_summary.html"
            ),
            caption="../report/cellranger_count.rst",
            category="cellranger",
            subcategory="count report",
            labels={"sample": "{sample_id}"},
        ),
    log:
        "logs/cellranger/multi/summary_files/summaries_{pool_id}_{sample_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"


# multiplexing outputs, according to:
# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-3p-outputs-cellplex


rule cellranger_multi_files_multiplexing_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update("results/cellranger/{pool_id}/outs/multi/count/feature_reference.csv"),
        update(
            "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/assignment_confidence_table.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/cells_per_tag.json"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/tag_calls_summary.csv"
        ),
    log:
        "logs/cellranger/multi/multiplexing_files/multiplexing_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"


rule cellranger_multi_files_multiplexing_per_sample:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/feature_reference.csv"
        ),
    log:
        "logs/cellranger/multi/multiplexing_files/multiplexing_per_sample_{pool_id}_{sample_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"


rule cellranger_multi_files_multiplexing_antibody_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update(
            "results/cellranger/{pool_id}/outs/multi/count/antibody_analysis/aggregate_barcodes.csv"
        ),
    log:
        "logs/cellranger/multi/multiplexing_files/multiplexing_antibody_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"


rule cellranger_multi_files_multiplexing_crispr_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/cells_per_protospacer.json"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/feature_reference.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/perturbation_effects_by_feature"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/perturbation_effects_by_target"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/perturbation_efficiencies_by_feature.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/perturbation_efficiencies_by_target.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/protospacer_calls_per_cell.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/protospacer_calls_summary.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/protospacer_umi_thresholds.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/crispr_analysis/protospacer_umi_thresholds.json"
        ),
    log:
        "logs/cellranger/multi/multiplexing_files/multiplexing_crispr_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"


# "Gene Expression" output files


rule cellranger_multi_files_gene_expression_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update("results/cellranger/{pool_id}/outs/multi/count/raw_molecule_info.h5"),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix/matrix.mtx.gz"
        ),
        update("results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix.h5"),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/unassigned_alignments.bam"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/count/unassigned_alignments.bam.bai"
        ),
    log:
        "logs/cellranger/multi/gene_expression_files/gex_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"


rule cellranger_multi_files_gene_expression_per_sample:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_barcodes.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_alignments.bam"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_alignments.bam.bai"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/features.tsv.gz"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix.h5"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_molecule_info.h5"
        ),
    log:
        "logs/cellranger/multi/gene_expression_files/gex_per_sample_{pool_id}_{sample_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"


# VDJ output files


rule cellranger_multi_files_vdj_reference:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update("results/cellranger/{pool_id}/outs/vdj_reference/reference.json"),
        update("results/cellranger/{pool_id}/outs/vdj_reference/fasta/regions.fa"),
    log:
        "logs/cellranger/multi/vdj_reference_files_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"


rule cellranger_multi_files_vdj_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update(
            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig_annotations.bed"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig_annotations.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig_annotations.json"
        ),
        update("results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.bam"),
        update("results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.bam.bai"),
        update("results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.fasta"),
        update(
            "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.fasta.fai"
        ),
        update("results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.fastq"),
    log:
        "logs/cellranger/multi/{vdj_type}_files/{vdj_type}_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"


rule cellranger_multi_files_vdj_per_sample:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/airr_rearrangement.tsv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/cell_barcodes.json"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/clonotypes.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/concat_ref.bam"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/concat_ref.bam.bai"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/concat_ref.fasta"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/concat_ref.fasta.fai"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus.bam"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus.bam.bai"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus.fasta"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus.fasta.fai"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/consensus_annotations.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/donor_regions.fa"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/filtered_contig_annotations.csv"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/filtered_contig.fasta"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/filtered_contig.fastq"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/vdj_contig_info.pb"
        ),
        update(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/{vdj_type}/vloupe.vloupe"
        ),
    log:
        "logs/cellranger/multi/{vdj_type}_files/{vdj_type}_per_sample_{pool_id}_{sample_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    localrule: True
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.sh"
