<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-rnaseqdge_logo_dark.png">
    <img alt="nf-core/rnaseqdge" src="docs/images/nf-core-rnaseqdge_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/rnaseqdge/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/rnaseqdge/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/rnaseqdge/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/rnaseqdge/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/rnaseqdge/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/rnaseqdge)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rnaseqdge-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/rnaseqdge)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/rnaseqdge** is a bioinformatics pipeline that analyses RNA sequencing data and performes differential gene expression analysis. It takes a samplesheet and FASTQ files as input, performs quality control, (pseudo-)alignment, produces a gene expression matrix, a QC report, and differentially expressed genes between two groups.

1. combine multiple RNAseq experiment per sample together
2. read quality control ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. summarize quality control of raw reads ([`MultiQC`](http://multiqc.info/))
4. if reference fasta is not given, download latest human reference genome and transcriptome from Ensembl
5. multiple gene quantification routes:
   - aligner (slow):
     - [`STAR`](https://github.com/alexdobin/STAR) -> [`Salmon`](https://combine-lab.github.io/salmon/)
     - [`STAR`](https://github.com/alexdobin/STAR) -> [`RSEM`](https://github.com/deweylab/RSEM)
   - pseudo-aligner (fast):
     - [`Salmon`](https://combine-lab.github.io/salmon/)
     - [`Kallisto`](https://pachterlab.github.io/kallisto/)
6. differential gene expression analysis using [`DESeq2`](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) and [`edgeR`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,group
c1,Control1_1.fq.gz,Control1_2.fq.gz,control
c2,Control2_1.fq.gz,Control2_2.fq.gz,control
t1,Tumor1_1.fq.gz,Tumor1_2.fq.gz,tumor
t2,Tumor2_1.fq.gz,Tumor2_2.fq.gz,tumor
```

Each row represents a pair of fastq files (currently only paired end supported).

Now, you can run the pipeline using:

```bash
nextflow run nf-core/rnaseqdge \
   -profile docker \
   --aligner <star_rsem|star_salmon|kallisto|salmon> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

(Note: `star_rsem` and `star_salmon` index file creation can take multiple hours and is RAM intensive.)

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/rnaseqdge/usage) and the [parameter documentation](https://nf-co.re/rnaseqdge/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/rnaseqdge/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/rnaseqdge/output).

## Possible improvements:

- single-end reads
- include pre-made indeces of aligners to avoid construction
- more preprocessing: adapter trimming, removal ribosomal RNA, ...
- make salmon aligner decoy-aware
- add test run option
- add ext.args to config
- add covariables to DGE analysis
- add 3 or more group comparisons for DGE
- add singularity, conda, etc ..
- proper software versions output

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/rnaseqdge for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
