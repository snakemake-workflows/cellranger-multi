## Workflow overview

This workflow is a best-practice snakemake workflow for systematically running `cellranger multi` on one or more samples.
See the [10X documentation for choosing a pipeline](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-choosing-a-pipeline) to see whether this is the preprocessing you need.
If your assay setup suggests `cellranger count`, have a look at the [standardised workflow for `cellranger count` instead](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/cellranger-count).

The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Link in files to a new file name that follows cellranger requirements.
2. Create a per-sample cellranger multi config CSV sheet.
3. Run `cellranger multi`, parallelizing over biological samples.
4. Create a snakemake report with the a Web Summary per biological sample.

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
The conda environments are stored in the hidden `.snakemake/conda/` folder.
You can usually identify the exact conda environment used by a rule from the `.snakemake/logs/` files or the respective cluster system log files.
Search for the execution of the respective rule (`cellranger_multi`) and then look for `Activating conda environment:` right below.
You can then delete the respective file and directory under `.snakemake/conda/` and rerun the workflow.

### Input data

#### pool sheet TSV

The sample sheet configures all the possible [columns for the `[libraries]` section of the multi config CSV file](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts#libraries):

| id     | feature_types   | read1                                                     | read2                                                     | lane_number |
| ------ | --------------- | --------------------------------------------------------- | --------------------------------------------------------- | ----------- |
| pool_1 | Gene Expression | ../data/pool1_gex/sample1_gex.bwa.L001.read1.fastq.gz   | ../data/pool1_gex/sample1_gex.bwa.L001.read2.fastq.gz   |           1 |
| pool_1 | Gene Expression | ../data/pool1_gex/pool1_gex.bwa.L002.read1.fastq.gz   | ../data/pool1_gex/pool1_gex.bwa.L002.read2.fastq.gz   |           2 |
| pool_1 | VDJ-T           | ../data/pool1_vdjt/pool1_vdjt.bwa.L003.read1.fastq.gz | ../data/pool1_vdjt/pool1_vdjt.bwa.L003.read2.fastq.gz |           1 |

For more details on these columns, refer to the [10X documentation for the `[libraries]` section of the multi config CSV file](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts#libraries).
We also provide specific subsection links wherever available.

These are **required columns**:

* `id` is an arbitrary name assigned to represent one pool of samples.
  One multi config CSV file per `id` value will be created by the workflow, and each `id` will be processed separately, so that the workflow can parallelize as much as possible.
  If you used just one sample to create your assay libraries (without any multiplexing), then the `id` value used here can just be a sample name.
  If you used multiplexing barcodes to pool multiple samples, this id groups all of the samples that were pooled before preparing assay libraries.
  For multiplexed experiments, you have to provide a [`multiplexing TSV` file (see below)](#multiplexing-tsv) that provides a barcode or id for each sample in a multiplexed pool.
* `feature_types` can be any of the [values listed in the `cellranger multi` documentation on multi config CSVs](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts#feature-types).
  The only disallowed feature type is `VDJ`, as this results in auto-detection of the actual `VDJ-` feature type and thus to unpredictable folder names in the output.
* `read1` and `read2` require file names with paths relative to the main workflow directory (the directory, where you run the `snakemake` command).
  From these (and the optional `lane_number` column), the raw read data files are linked into the folder and file name structure that cellranger expects, and the `fastq_id` and `fastqs` columns of the multi config CSV file are set up accordingly.

These are **optional columns**:

* `lane_number` is only necessary if a single sample is sequenced across multiple lanes.
  Usually, you will number lanes starting from 1 and only up to a single digit number of lanes.
  As we specify one pair of fastq files per row, the `lane_number` column also only contains a single lane number, as we have one pair of files per lane.
  For the `lanes` column in the  final multi config CSV file, multiple lane numbers get parsed into the format `1|2|3` etc.
* `physical_library_id` is usually auto-detected, so just omit it if in doubt.
* `subsample_rate` is not usually needed.
* `chemistry` is `auto` per default and only applicable for Flex assays.
  If you think this applies to your setup, see the [`chemistry` options in the 10X documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-multi-config-csv-opts#chem-opts).

#### multiplexing TSV

This file is only necessary if you used [multiplexing to pool multiple samples for sequencing](https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-3p-what-is-cellplex#overview).
You provide its relative path in the global `config/config.yaml` file under `multi_config_csv_sections: multiplexing:`.
It specifies which sample used which `ocm_barcode_ids` (or `hashtag_ids` or `cmo_ids` or `probe_barcode_ids`).

| id     | sample_id | ocm_barcode_ids |
| ------ | --------- | --------------- |
| pool_1 | sample_1  | OB1\|OB2        |
| pool_1 | sample_2  | OB3\|OB4        |

These are **required columns**:
* `sample_id` is an arbitrary name assigned to represent one biological sample that was sequenced with a known set of barcodes.
* `id` is an arbitrary name assigned to represent one pool of samples.
  It needs to match an `id` from the [pool sheet TSV file](#pool-sheet-tsv).
  A pool can contain one or more samples.
* One of the following columns, depending on the barcode type used for the multiplexing setup:
  * `ocm_barcode_ids`
  * `hashtag_ids`
  * `cmo_ids`
  * `probe_barcode_ids`

You can read the [cellranger documentation on multiplexing setups](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi) to determine which column to use.
In addition, the [samples section of the multi config csv file documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts#samples) explains which values you can use there and which additional columns may be available.


### Global analysis-level configuration

All global configuration settings for the whole analysis are specified in the `config/config.yaml` file.
This file is extensively commented to explain how to set which options.
You can delete any options you don't need, or set them to an empty string (`""`).
The only required sections are those for the feature types present in the `feature_types` column of the sample sheet and the flag whether `multiplexing` should be `activated`.
And the only required entry for a required section is usually the `reference:` path or file specification.