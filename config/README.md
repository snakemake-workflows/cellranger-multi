## Workflow overview

This workflow is a best-practice workflow for systematically running `cellranger multi` on one or more samples.
See the [10X documentation choosing a pipeline](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-choosing-a-pipeline) to see whether this is the preprocessing you need.
If your assay setup suggests `cellranger multi`, have a look at the [standardised workflow for `cellranger multi` instead](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/cellranger-multi).

The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Link in files to a new file name that follows cellranger requirements.
2. Create a per-sample cellranger library CSV sheet.
3. Run `cellranger multi`, parallelizing over samples.
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
So once your specific analysis has created this conda environment, the cellranger version will stay at the version specified at that time.
Should you ever want to update the cellranger version for an analysis, you will have to update the `CELLRANGER_TARBALL` environment variable and delete the conda environment, to ensure that it gets re-generated.

### Input data

The sample sheet configures all the possible [columns for the `[libraries]` section of the multi config CSV file](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts#libraries):

| sample  | feature_types   | read1                                                     | read2                                                     | lane_number |
| ------- | --------------- | --------------------------------------------------------- | --------------------------------------------------------- | ----------- |
| sample1 | Gene Expression | ../data/sample1_gex/sample1_gex.bwa.L001.read1.fastq.gz   | ../data/sample1_gex/sample1_gex.bwa.L001.read2.fastq.gz   |           1 |
| sample1 | VDJ-T           | ../data/sample1_vdjt/sample1_vdjt.bwa.L003.read1.fastq.gz | ../data/sample1_vdjt/sample1_vdjt.bwa.L003.read2.fastq.gz |           1 |
| sample2 | Gene Expression | ../data/sample2_gex/sample2_gex.bwa.L001.read1.fastq.gz   | ../data/sample2_gex/sample2_gex.bwa.L001.read2.fastq.gz   |           1 |
| sample2 | Gene Expression | ../data/sample2_gex/sample2_gex.bwa.L002.read1.fastq.gz   | ../data/sample2_gex/sample2_gex.bwa.L002.read2.fastq.gz   |           2 |


For more details on these columns, refer to the [10X documentation for the `[libraries]` section of the multi config CSV file](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts#libraries).
We also provide specific subsection links wherever available.

These are **required columns**:

* `sample` can be any sample name you want to assign, where one sample groups all of the lane fastq pairs per assay that was performed on a particular biological sample.
  Thus, the name should usually contain an identifier for the biological sample and for the assay type, for example `replicate_1_5gex`.
* `feature_types` can be any of the [values listed in the `cellranger multi` documentation on multi config CSVs](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts#feature-types).
* `read1` and `read2` require file names with paths relative to the main workflow directory (where you run the `snakemake` command).
  From these (and the `lane_number` column), the raw read data files are linked into the folder and file name structure that cellranger expects, and the `fastq_id` and `fastqs` columns of the multi config CSV file are set up accordingly.

These are **optional columns**:
* `lane_number` is only necessary if a single sample is sequenced across multiple lanes.
  Usually, you will number starting from 1 and only up to a single digit number of lanes.
  As we specify one pair of fastq files per row, the `lane_number` column also only contains a single lane number, as we have one pair of files per lane.
  For the `lanes` column in the  final multi config CSV file, multiple lane numbers get parsed into the format `1|2|3` etc.
* `physical_library_id` is usually auto-detected, so just omit it if in doubt.
* `subsample_rate` is not usually needed.
* `chemistry` is `auto` per default and only applicable for Flex assays.
  If you think this applies to your setup, see the [`chemistry` options in the 10X documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-multi-config-csv-opts#chem-opts).

### Global analysis-level configuration

All global configuration settings for the whole analysis are specified in the `config/config.yaml` file.
This file is extensively commented to explain how to set which options.