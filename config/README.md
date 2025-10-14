## Workflow overview

This workflow is a best-practice workflow for systematically running `cellranger count` on one or more samples.
See the [10X documentation choosing a pipeline](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-choosing-a-pipeline) to see whether this is the preprocessing you need.
If your assay setup suggests `cellranger multi`, have a look at the [standardised workflow for `cellranger multi` instead](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/cellranger-multi).

The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Link in files to a new file name that follows cellranger requirements.
2. Create a per-sample cellranger library CSV sheet.
3. Run `cellranger count`, parallelizing over samples.
4. Create a snakemake report with the Web Summaries.

## Running the workflow

### cellranger download

As a pre-requisite for running the workflow, you need to download the `*.tar.gz` file with the Cell Ranger executable from the Cell Ranger Download center:
https://www.10xgenomics.com/support/software/cell-ranger/downloads

Afterwards, set the environment variable `CELLRANGER_TARBALL` to the full path of this executable, for example:
```{bash}
export CELLRANGER_TARBALL="/absolute/path/to/cellranger-8.1.1.tar.gz"
```
To make this a permanently set environment variable for your user on the respective system, add the (adapted) line from above to your `~/.bashrc` file and make sure this file is always loaded.

With this environment variable set, the workflow will automatically install `cellranger` into a conda environment that is then used for all cellranger steps.

### Input data

The sample sheet has the following layout:

| sample  | lane_number | library_type    | read1                                   | read2                                   |
| ------- | ----------- | --------------- | --------------------------------------- | --------------------------------------- |
| sample1 |           1 | Gene Expression | ../data/sample1.bwa.L001.read1.fastq.gz | ../data/sample1.bwa.L001.read2.fastq.gz |
| sample1 |           2 | Gene Expression | ../data/sample1.bwa.L002.read1.fastq.gz | ../data/sample1.bwa.L002.read2.fastq.gz |
| sample2 |           1 | Gene Expression | ../data/sample2.bwa.read1.fastq.gz      | ../data/sample2.bwa.read2.fastq.gz      |

The `lane_number` column is optional, and only necessary if a single sample is sequenced across multiple lanes.
All other columns are required:

* `library_type` can be any of the [values listed in the `cellranger count` documentation on Library CSVs](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-libraries-csv).
* `read1` and `read2` require file names with paths relative to the main workflow directory (where you run the `snakemake` command).

### Parameters

This table lists the most important configuration parameters that can be set in the `config/config.yaml` file.

The `ref_data` needs to be downloaded manually from the Cell Ranger Download Center:
https://www.10xgenomics.com/support/software/cell-ranger/downloads
After download, extract the tar file into a directory and provide the directory's path under `ref_data:`.

You can also check with your local compute cluster, if they have the reference data available already.
In that case, you can just point the `ref_data:` configuration variable to the respective path.

| parameter          | type | details                                      | default                        |
| ------------------ | ---- | -------------------------------------------- | ------------------------------ |
| **sample_sheet**   |      |                                              |                                |
| path               | str  | path to sample sheet, mandatory              | "config/samples.tsv"           |
| **ref_data**       |      |                                              |                                |
| path               | str  | path to downloaded reference data, mandatory |                                |
