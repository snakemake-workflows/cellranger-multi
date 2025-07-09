## Workflow overview

This workflow is a best-practice workflow for systematically running `cellranger count` on one or more samples.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Link in files to a new file name that follows cellranger requirements.
2. Create a per-sample cellranger library CSV sheet.
3. Run cellranger count, parallelizing over samples.
4. Create a snakemake report with the Web Summaries.

## Running the workflow

### Input data

The sample sheet has the following layout:

| sample  | lane_number | library_type    | read1                           | read2                           |
| ------- | ----------- | --------------- | ------------------------------- | ------------------------------- |
| sample1 |           1 | Gene Expression | sample1.bwa.L001.read1.fastq.gz | sample1.bwa.L001.read2.fastq.gz |
| sample1 |           2 | Gene Expression | sample1.bwa.L002.read1.fastq.gz | sample1.bwa.L002.read2.fastq.gz |
| sample2 |           1 | Gene Expression | sample2.bwa.read1.fastq.gz      | sample2.bwa.read2.fastq.gz      |

The `lane` column is optional, and only necessary if a single sample is sequenced across multiple lanes.
All other columns are required.

### Parameters

This table lists all parameters that can be used to run the workflow.

| parameter          | type | details                                      | default                        |
| ------------------ | ---- | -------------------------------------------- | ------------------------------ |
| **sample_sheet**   |      |                                              |                                |
| path               | str  | path to sample sheet, mandatory              | "config/samples.tsv"           |
| **ref_data**       |      |                                              |                                |
| path               | str  | path to downloaded reference data, mandatory |                                |
