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
    localrule: True
    conda:
        "../envs/bash_coreutils.yaml"
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
            default=None,
        ),
        multiplexing=branch(
            lookup(
                within=config,
                dpath="multi_config_csv_sections/multiplexing/activate",
            ),
            lookup(
                within=config,
                dpath="multi_config_csv_sections/multiplexing/tsv",
                default=None,
            ),
            None
        ),
    output:
        multi_config_csv="results/input/{pool_id}.cell_ranger_multi_config.csv",
    log:
        "logs/input/{pool_id}.cell_ranger_multi_config.log",
    localrule: True
    conda:
        "../envs/tidyverse.yaml"
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
            default=None,
        ),
    output:
        "results/cellranger/{pool_id}/outs/config.csv",
        out_dir=directory("results/cellranger/{pool_id}/outs/"),
    log:
        "logs/cellranger/multi/multi_run_{pool_id}.log",
    conda:
        "../envs/cellranger.yaml"
    threads: 8
    resources:
        mem_mb=lambda wc, threads: threads * 4000,
    params:
        mem_gb=lambda wc, resources: math.floor(resources.mem_mb / 1000),
        out_dir=lambda wc, output: path.abspath(
            path.dirname(output["out_dir"]).removesuffix("outs")
        ),
    shell:
        "(rm -rf {params.out_dir}; "
        " cellranger multi "
        "  --id={wildcards.pool_id} "
        "  --output-dir={params.out_dir} "
        "  --csv={input.multi_config_csv} "
        "  --localcores={threads} "
        "  --localmem={params.mem_gb}; "
        ") >{log} 2>&1 "


rule cellranger_multi_files_summaries:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/metrics_summary.csv",
        report(
            "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/web_summary.html",
            caption="../report/cellranger_count.rst",
            category="cellranger",
            subcategory="count report",
            labels={"sample": "{sample_id}"},
        ),
    log:
        "logs/cellranger/multi/summary_files/summaries_{pool_id}_{sample_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"


# multiplexing outputs, according to:
# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-3p-outputs-cellplex


rule cellranger_multi_files_multiplexing_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        "results/cellranger/{pool_id}/outs/multi/count/feature_reference.csv",
        "results/cellranger/{pool_id}/outs/multi/count/raw_cloupe.cloupe",
        "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/assignment_confidence_table.csv",
        "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/cells_per_tag.json",
        "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv",
        "results/cellranger/{pool_id}/outs/multi/multiplexing_analysis/tag_calls_summary.csv",
    log:
        "logs/cellranger/multi/multiplexing_files/multiplexing_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"


rule cellranger_multi_files_multiplexing_per_sample:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/feature_reference.csv",
    log:
        "logs/cellranger/multi/multiplexing_files/multiplexing_per_sample_{pool_id}_{sample_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"


rule cellranger_multi_files_multiplexing_antibody_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        "results/cellranger/{pool_id}/outs/multi/count/antibody_analysis/aggregate_barcodes.csv",
    log:
        "logs/cellranger/multi/multiplexing_files/multiplexing_antibody_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"


rule cellranger_multi_files_multiplexing_crispr_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
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
    log:
        "logs/cellranger/multi/multiplexing_files/multiplexing_crispr_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"


# "Gene Expression" output files


rule cellranger_multi_files_gene_expression_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        "results/cellranger/{pool_id}/outs/multi/count/raw_molecule_info.h5",
        "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz",
        "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz",
        "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix/matrix.mtx.gz",
        "results/cellranger/{pool_id}/outs/multi/count/raw_feature_bc_matrix.h5",
        "results/cellranger/{pool_id}/outs/multi/count/unassigned_alignments.bam",
        "results/cellranger/{pool_id}/outs/multi/count/unassigned_alignments.bam.bai",
    log:
        "logs/cellranger/multi/gene_expression_files/gex_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"


rule cellranger_multi_files_gene_expression_per_sample:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_barcodes.csv",
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_alignments.bam",
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_alignments.bam.bai",
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz",
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/features.tsv.gz",
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz",
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_filtered_feature_bc_matrix.h5",
        "results/cellranger/{pool_id}/outs/per_sample_outs/{sample_id}/count/sample_molecule_info.h5",
    log:
        "logs/cellranger/multi/gene_expression_files/gex_per_sample_{pool_id}_{sample_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"


# VDJ output files


rule cellranger_multi_files_vdj_reference:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        "results/cellranger/{pool_id}/outs/vdj_reference/reference.json",
        "results/cellranger/{pool_id}/outs/vdj_reference/fasta/regions.fa",
    log:
        "logs/cellranger/multi/vdj_reference_files_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"


rule cellranger_multi_files_vdj_global:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
        "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig_annotations.bed",
        "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig_annotations.csv",
        "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig_annotations.json",
        "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.bam",
        "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.bam.bai",
        "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.fasta",
        "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.fasta.fai",
        "results/cellranger/{pool_id}/outs/multi/{vdj_type}/all_contig.fastq",
    log:
        "logs/cellranger/multi/{vdj_type}_files/{vdj_type}_global_{pool_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"


rule cellranger_multi_files_vdj_per_sample:
    input:
        csv_copy="results/cellranger/{pool_id}/outs/config.csv",
    output:
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
    log:
        "logs/cellranger/multi/{vdj_type}_files/{vdj_type}_per_sample_{pool_id}_{sample_id}.log",
    conda:
        "../envs/bash_coreutils.yaml"
    threads: 1
    script:
        "../scripts/check_cellranger_outputs.bash"
