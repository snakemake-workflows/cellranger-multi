# Snakemake workflow: `cellranger-count`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/cellranger-count/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/snakemake-workflows/cellranger-count/actions/workflows/main.yml)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/cellranger-count)

A Snakemake workflow for counting single cell RNAseq (scRNA-seq) data with Cell Ranger (Cell Ranger licensing requires a manual download of the software). 

- [Snakemake workflow: `cellranger-count`](#snakemake-workflow-name)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/cellranger-count).

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options

To run the workflow from command line, change the working directory.

```bash
cd path/to/snakemake-workflow-name
```

Adjust options in the default config file `config/config.yml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the workflow with test files using **conda** (this will download several GB of example data):

```bash
snakemake --cores 2 --sdm conda --directory .test/ --snakefile .test/Snakefile
```

## Authors

- David Lähnemann
  - German Cancer Consortium (DKTK), partner site Essen-Düsseldorf, A partnership between DKFZ and University Hospital Essen
  - https://orcid.org/0000-0002-9138-4112
- Smiths Sengkwawoh Lueong
  - German Cancer Consortium (DKTK), partner site Essen-Düsseldorf, A partnership between DKFZ and University Hospital Essen
  - https://orcid.org/0000-0002-2776-6706

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.