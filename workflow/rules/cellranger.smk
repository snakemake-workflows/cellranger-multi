# Link the raw fastq files into an input directory and
# make sure they are named EXACTLY how Cell Ranger needs
# them. Sigh.
# -----------------------------------------------------
rule follow_pedantic_cell_ranger_naming_scheme:
    input:
        fq1=lambda wc: get_input_file(wc, "read1"),
        fq2=lambda wc: get_input_file(wc, "read2"),
    output:
        fq1="results/input/{sample}_{feature_type}/{sample}_S1_L00{lane_number}_R1_001.fastq.gz",
        fq2="results/input/{sample}_{feature_type}/{sample}_S1_L00{lane_number}_R2_001.fastq.gz",
    log:
        "logs/input/{sample}_{feature_type}/{sample}_{feature_type}_S1_L00{lane_number}_001.log",
    localrule: True
    conda:
        "../envs/coreutils.yaml"
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
        sample_sheet=lookup(within=config, dpath="sample_sheet"),
        fq1=lambda wc: get_sample_fastqs(wc, "R1"),
        fq2=lambda wc: get_sample_fastqs(wc, "R2"),
        # use the sample_sheet as an existing dummy placeholder, in case no
        # feature reference file is specified for this analysis
        feature_reference=lookup(
            within=config,
            dpath="multi_config_csv_sections/feature/reference",
            default=lookup(within=config, dpath="sample_sheet"),
        ),
    output:
        multi_config_csv="results/input/{sample}.cell_ranger_multi_config.csv",
    log:
        "logs/input/{sample}.cell_ranger_multi_config.log",
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
        multi_config_csv="results/input/{sample}.cell_ranger_multi_config.csv",
        fq1=lambda wc: get_sample_fastqs(wc, "R1"),
        fq2=lambda wc: get_sample_fastqs(wc, "R2"),
        # use the multi_config_csv as an existing dummy placeholder, in case no
        # reference is needed here (if no Gene Expression samples present)
        reference=lookup(
            within=config,
            dpath="multi_config_csv_sections/gene-expression/reference",
            default="results/input/{sample}.cell_ranger_multi_config.csv",
        ),
    output:
        "results/cellranger/{sample}/outs/config.csv",
        directory("results/cellranger/{sample}/outs/"),
    log:
        "logs/cellranger/multi/multi_run_{sample}.log",
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
        "(rm -r {params.out_dir}; "
        " cellranger multi "
        "  --id={wildcards.sample} "
        "  --output-dir={params.out_dir} "
        "  --csv={input.multi_config_csv} "
        "  --localcores={threads} "
        "  --localmem={params.mem_gb}; "
        ") >{log} 2>&1 "


rule cellranger_multi_files_summaries:
    input:
        csv_copy="results/cellranger/{sample}/outs/config.csv",
    output:
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/metrics_summary.csv",
        report(
            "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/web_summary.html",
            caption="../report/cellranger_count.rst",
            category="cellranger",
            subcategory="count report",
            labels={"sample": "{multiplex_sample}"},
        ),
    log:
        "logs/cellranger/multi/summary_files/summaries_{sample}_{multiplex_sample}.log",
    threads: 1
    shell:
        "( touch {output} "
        ") >{log} 2>&1 "


rule cellranger_multi_files_gene_expression_global:
    input:
        csv_copy="results/cellranger/{sample}/outs/config.csv",
    output:
        "results/cellranger/{sample}/outs/multi/count/raw_molecule_info.h5",
        "results/cellranger/{sample}/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz",
        "results/cellranger/{sample}/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz",
        "results/cellranger/{sample}/outs/multi/count/raw_feature_bc_matrix/matrix.mtx.gz",
        "results/cellranger/{sample}/outs/multi/count/raw_feature_bc_matrix.h5",
        "results/cellranger/{sample}/outs/multi/count/unassigned_alignments.bam",
        "results/cellranger/{sample}/outs/multi/count/unassigned_alignments.bam.bai",
    log:
        "logs/cellranger/multi/gene_expression_files/gex_global_{sample}.log",
    threads: 1
    shell:
        "( touch {output} "
        ") >{log} 2>&1 "


rule cellranger_multi_files_gene_expression_per_sample:
    input:
        csv_copy="results/cellranger/{sample}/outs/config.csv",
    output:
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/count/sample_filtered_barcodes.csv",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/count/sample_alignments.bam",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/count/sample_alignments.bam.bai",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/count/sample_filtered_feature_bc_matrix/features.tsv.gz",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/count/sample_filtered_feature_bc_matrix.h5",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/count/sample_molecule_info.h5",
    log:
        "logs/cellranger/multi/gene_expression_files/gex_per_sample_{sample}_{multiplex_sample}.log",
    threads: 1
    shell:
        "( touch {output} "
        ") >{log} 2>&1 "


rule cellranger_multi_files_vdj_reference:
    input:
        csv_copy="results/cellranger/{sample}/outs/config.csv",
    output:
       "results/cellranger/{sample}/outs/vdj_reference/reference.json",
       "results/cellranger/{sample}/outs/vdj_reference/fasta/regions.fa",
    log:
        "logs/cellranger/multi/vdj_reference_files_{sample}.log",
    threads: 1
    shell:
        "( touch {output} "
        ") >{log} 2>&1 "


rule cellranger_multi_files_vdj_global:
    input:
        csv_copy="results/cellranger/{sample}/outs/config.csv",
    output:
       "results/cellranger/{sample}/outs/multi/{vdj_type}/all_contig_annotations.bed",
       "results/cellranger/{sample}/outs/multi/{vdj_type}/all_contig_annotations.csv",
       "results/cellranger/{sample}/outs/multi/{vdj_type}/all_contig_annotations.json",
       "results/cellranger/{sample}/outs/multi/{vdj_type}/all_contig.bam",
       "results/cellranger/{sample}/outs/multi/{vdj_type}/all_contig.bam.bai",
       "results/cellranger/{sample}/outs/multi/{vdj_type}/all_contig.fasta",
       "results/cellranger/{sample}/outs/multi/{vdj_type}/all_contig.fasta.fai",
       "results/cellranger/{sample}/outs/multi/{vdj_type}/all_contig.fastq",
    log:
        "logs/cellranger/multi/{vdj_type}_files/{vdj_type}_global_{sample}.log",
    threads: 1
    shell:
        "( touch {output} "
        ") >{log} 2>&1 "


rule cellranger_multi_files_gene_expression_per_sample:
    input:
        csv_copy="results/cellranger/{sample}/outs/config.csv",
    output:
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/airr_rearrangement.tsv",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/cell_barcodes.json",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/clonotypes.csv",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/concat_ref.bam",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/concat_ref.bam.bai",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/concat_ref.fasta",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/concat_ref.fasta.fai",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/consensus.bam",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/consensus.bam.bai",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/consensus.fasta",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/consensus.fasta.fai",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/consensus_annotations.csv",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/donor_regions.fa",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/filtered_contig_annotations.csv",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/filtered_contig.fasta",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/filtered_contig.fastq",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/vdj_contig_info.pb",
        "results/cellranger/{sample}/outs/per_sample_outs/{multiplex_sample}/{vdj_type}/vloupe.vloupe",
    log:
        "logs/cellranger/multi/{vdj_type}_files/{vdj_type}_per_sample_{sample}_{multiplex_sample}.log",
    threads: 1
    shell:
        "( touch {output} "
        ") >{log} 2>&1 "
