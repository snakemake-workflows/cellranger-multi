# Link the raw fastq files into an input directory and
# make sure they are named EXACTLY how Cell Ranger needs
# them. Sigh.
# -----------------------------------------------------
rule follow_pedantic_cell_ranger_naming_scheme:
    input:
        fq1=lambda wc: get_input_file(wc, "read1"),
        fq2=lambda wc: get_input_file(wc, "read2"),
    output:
        fq1="results/input/{sample}_S1_L00{lane_number}_R1_001.fastq.gz",
        fq2="results/input/{sample}_S1_L00{lane_number}_R2_001.fastq.gz",
    log:
        "logs/input/{sample}_S1_L00{lane_number}_001.log",
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


# Create a libary CSV file for Cell Ranger.
# -----------------------------------------------------
rule create_cellranger_library_csv:
    input:
        sample_sheet=lookup(within=config, dpath="sample_sheet"),
        fq1=lambda wc: get_sample_fastqs(wc, "R1"),
        fq2=lambda wc: get_sample_fastqs(wc, "R2"),
    output:
        library_csv="results/input/{sample}.cell_ranger_library.csv",
    log:
        "logs/input/{sample}.cell_ranger_library.log",
    localrule: True
    conda:
        "../envs/tidyverse.yaml"
    params:
        fastqs_dir=lambda wc, input: path.abspath(path.dirname(input.fq1[0])),
    script:
        "../scripts/create_cellranger_library_csv.R"


# Run cellranger count on one sample.
# -----------------------------------------------------
rule cellranger_count:
    input:
        library_csv="results/input/{sample}.cell_ranger_library.csv",
        fq1=lambda wc: get_sample_fastqs(wc, "R1"),
        fq2=lambda wc: get_sample_fastqs(wc, "R2"),
        ref_data=lookup(within=config, dpath="ref_data"),
    output:
        "results/cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        "results/cellranger/{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        "results/cellranger/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        "results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
        "results/cellranger/{sample}/outs/metrics_summary.csv",
        "results/cellranger/{sample}/outs/molecule_info.h5",
        "results/cellranger/{sample}/outs/possorted_genome_bam.bam",
        "results/cellranger/{sample}/outs/possorted_genome_bam.bam.bai",
        "results/cellranger/{sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz",
        "results/cellranger/{sample}/outs/raw_feature_bc_matrix/features.tsv.gz",
        "results/cellranger/{sample}/outs/raw_feature_bc_matrix/matrix.mtx.gz",
        "results/cellranger/{sample}/outs/raw_feature_bc_matrix.h5",
        report(
            "results/cellranger/{sample}/outs/web_summary.html",
            caption="../report/cellranger_count.rst",
            category="cellranger",
            subcategory="count report",
            labels={"sample": "{sample}"},
        ),
        out_dir=directory("results/cellranger/{sample}/outs/"),
    log:
        "logs/cellranger/{sample}.log",
    conda:
        "../envs/cellranger.yaml"
    threads: 8
    resources:
        mem_mb=lambda wc, threads: threads * 4000,
    params:
        mem_gb=lambda wc, resources: int(resources.mem_mb / 1000),
        out_dir=lambda wc, output: path.abspath(
            path.dirname(output["out_dir"]).removesuffix("outs")
        ),
    shell:
        "(rm -r {params.out_dir}; "
        " cellranger count "
        "  --id={wildcards.sample} "
        "  --output-dir={params.out_dir} "
        "  --transcriptome={input.ref_data} "
        "  --libraries={input.library_csv} "
        "  --nosecondary "
        "  --create-bam=true "
        "  --localcores={threads} "
        "  --localmem={params.mem_gb}; "
        ") >{log} 2>&1 "
