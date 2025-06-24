# Link the raw fastq files into an input directory and
# make sure they are named EXACTLY how Cell Ranger needs
# them. Sigh.
# -----------------------------------------------------
rule follow_pedantic_cell_ranger_naming_scheme:
    input:
        fq1=lambda wc: get_input_files(wc, "read1"),
        fq2=lambda wc: get_input_files(wc, "read2"),
    output:
        fq1="results/input/{sample}_S1_L00{lane_number}_R1_001.fastq.gz",
        fq2="results/input/{sample}_S1_L00{lane_number}_R2_001.fastq.gz",
    log:
        "logs/input/{sample}_S1_L00{lane_number}_001.log",
    shell:
        "( ln -S {input.fq1} {output.fq1}; "
        "  ln -S {input.fq2} {output.fq2}; "
        ") 2>&1 >{log}"


# Create a libary CSV file for Cell Ranger.
# -----------------------------------------------------
rule create_cellranger_library_csv:
    input:
        sample_sheet=lookup(within=config, dpath="sample_sheet"),
        fastqs=get_all_sample_fastqs,
    output:
        library_csv="results/input/cell_ranger_library.csv",
    log:
        "logs/input/cell_ranger_library.log",
    conda:
        "../envs/tidyverse.yaml"
    params:
        fastqs_dir=lambda wc, input: path.abspath(path.dirname(input.fastqs[0])),
    script:
        "../scripts/create_cellranger_library_csv.R"


# Run cellranger count on one sample.
# -----------------------------------------------------
rule cellranger_count:
    input:
        library_csv="results/input/cell_ranger_library.csv",
        fq1="results/input/{sample}_S1_L00{lane_number}_R1_001.fastq.gz",
        fq2="results/input/{sample}_S1_L00{lane_number}_R2_001.fastq.gz",
        ref_data=lookup(within=config, dpath="ref_data"),
    output:
        "results/cellranger/{sample}/outs/loupe.cloupe",
        "results/cellranger/{sample}/outs/iltered_feature_bc_matrix/barcodes.tsv.gz",
        "results/cellranger/{sample}/outs/iltered_feature_bc_matrix/features.tsv.gz",
        "results/cellranger/{sample}/outs/iltered_feature_bc_matrix/matrix.mtx.gz",
        "results/cellranger/{sample}/outs/iltered_feature_bc_matrix.h5",
        "results/cellranger/{sample}/outs/etrics_summary.csv",
        "results/cellranger/{sample}/outs/olecule_info.h5",
        "results/cellranger/{sample}/outs/ossorted_genome_bam.bam",
        "results/cellranger/{sample}/outs/ossorted_genome_bam.bam.bai",
        "results/cellranger/{sample}/outs/aw_feature_bc_matrix/barcodes.tsv.gz",
        "results/cellranger/{sample}/outs/aw_feature_bc_matrix/features.tsv.gz",
        "results/cellranger/{sample}/outs/aw_feature_bc_matrix/matrix.mtx.gz",
        "results/cellranger/{sample}/outs/aw_feature_bc_matrix.h5",
        "results/cellranger/{sample}/outs/eb_summary.html",
    log:
        "results/simulate_reads/{sample}.log",
    conda:
        "../envs/cellranger.yaml"
    params:
        read_number=lookup(within=config, dpath="simulate_reads/read_number"),
    threads: 8
    resources:
        mem_mb=lambda wc, threads: threads * 4000,
    params:
        mem_gb=lambda wc, resources: resources.mem_mb / 1000,
        out_dir=lambda wc, output: path.dirname(output[0]).removesuffix("outs/"),
    shell:
        "(cellranger count "
        "  --id={wildcards.sample} "
        "  --output=dir={params.out_dir} "
        "  --transcriptome={input.ref_data} "
        "  --libraries={input.library_csv} "
        "  --sample={wildcards.sample} "
        "  --create-bam=true "
        "  --localcores={threads} "
        "  --localmem={params.mem_gb}; "
        ") 2>&1 >{log} "
