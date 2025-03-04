[![version - AsaruSim](https://img.shields.io/badge/dynamic/yaml?url=https%3A%2F%2Fraw.githubusercontent.com%2FGenomiqueENS%2FAsaruSim%2Frefs%2Fheads%2Fmain%2Fversion.yml&query=%24.version&prefix=V&label=AsaruSim)](https://github.com/GenomiqueENS/AsaruSim/releases)
[![Made with Docker](https://img.shields.io/badge/Made_with-Docker-blue?logo=docker&logoColor=white)](https://www.docker.com/ "Go to Docker homepage")
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/6cb8154a35ae419689c35eeedcd71dac)](https://app.codacy.com/gh/alihamraoui/AsaruSim/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)

[![DOI](https://img.shields.io/badge/DOI-10.1093/bioinformatics/btaf087%20-green)](https://doi.org/10.1093/bioinformatics/btaf087)
[![Mastodon Follow](https://img.shields.io/mastodon/follow/111487733409234458?domain=genomic.social)](https://genomic.social/@Genomique_ENS)

# Documentation
<a href="images/asarusim_v2.svg"><img src="images/asarusim_v2.svg" align="middle" height="100" width="290" >

AsaruSim is an automated Nextflow workflow designed for simulating 10x single-cell long read data from the count matrix level to the sequence level. It aimed at creating a gold standard dataset for the assessment and optimization of single-cell long-read methods.
Full [documentation](https://GenomiqueENS.github.io/AsaruSim/) is avialable [here](https://GenomiqueENS.github.io/AsaruSim/).

<a href="images/schema.png"><img src="images/schema.png" align="middle" height="650" width="920" >

## $\textcolor{#FF7F00}{Requirements}$

This pipeline is powered by Nextflow workflow manager. All dependencies are automatically managed by Nextflow through a preconfigured Docker container, ensuring a seamless and reproducible installation process.

Before starting, ensure the following tools are installed and properly set up on your system:

- $\textcolor{#4DAF4A}{Nextflow}$ >= $\textcolor{#4DAF4A}{24.04.4}$ : A workflow engine for complex data pipelines. [Installation guide for Nextflow](https://www.nextflow.io/docs/latest/getstarted.html).
- $\textcolor{#4DAF4A}{Docker}$ or $\textcolor{#4DAF4A}{Singularity}$ : Containers for packaging necessary software, ensuring reproducibility. [Docker installation guide](https://docs.docker.com/get-docker/), [Singularity installation guide](https://sylabs.io/guides/3.0/user-guide/installation.html).

## $\textcolor{#FF7F00}{Installation}$

Clone the `AsaruSim` GitHub repository:

```bash
git clone https://github.com/GenomiqueENS/AsaruSim.git
cd AsaruSim
```

## $\textcolor{#FF7F00}{Test}$

To test your installation, we provide an automated script to download reference annotations and simulate a subset of human PBMC dataset `run_test.sh`.

```bash
bash run_test.sh
```

## $\textcolor{#FF7F00}{Configuration}$

Customize runs by editing the `nextflow.config` file and/or specifying parameters at the command line.

### Pipeline Input Parameters

Here are the primary input parameters for configuring the workflow:

### Main Parameters

| Parameter          | Description                                         | Format   | Default Value                                 |
|--------------------|-----------------------------------------------------|----------|-----------------------------------------------|
| `matrix`           | Path to the count matrix csv file (required)        |   .CSV       | `test_data/matrix.csv`                        |
| `transcriptome`    | Path to the reference transcriptome file (required) |    FASTA      | `test_data/transcriptome.fa`                  |
| `bc_counts`        | Path to the barcode count file (if no matrix provided).                   |    .CSV       | `test_data/test_bc.csv`                       |

### Optional Parameters

| Parameter          | Description                                         | Format   | Default Value                                 |
|--------------------|-----------------------------------------------------|----------|-----------------------------------------------|
| `features`         | Matrix feature counts                               |    STR      | `transcript_id`                               |
| `cell_types_annotation`    | Path to cell type annotation .csv file      |     CSV     | `null`                                        |
| `gtf`              | Path to transcriptom annotation .gtf file           |    GTF      | `null`                                        |
| `umi_duplication`    | UMI duplication     |     INT     | `0`                                        |
| `intron_retention`    | Simulate intron retention proces      |     BOOL     | `false`                                        |
| `ir_model`    | Intron retention MC model .CSV file     |     CSV     | `bin/models/SC3pv3_GEX_Human_IR_markov_model`                                        |
| `unspliced_ratio`    | percentage of transcrits to be unspliced     |     FLOAT     | `0.0`                                        |
| `ref_genome`       | reference genome .fasta file (if IR)                       | FASTA   | `null`                                       |
| `full_length`    | Indicates if transcripts are full length      |     BOOL     | `false`                                        |
| `truncation_model`    | Path to truncation probabilities .csv file      |     CSV     | `bin/models/truncation_default_model.csv`                                        |

### PCR Parameters

| Parameter          | Description                                            |   Format    | Default Value                                 |
|--------------------|--------------------------------------------------------|-------|-----------------------------------------------|
| `pcr_cycles`              | Number of PCR amplification cycles                                 |   INT    | `0`                                           |
| `pcr_error_rate`           | PCR error rate                           |   FLOAT    | `"0.0000001"`                                   |
| `pcr_dup_rate`      | PCR duplication rate                                   |    FLOAT   | `0.7`                              |
| `pcr_total_reads`      | Name of the project                                    |    INT   | `1000000`                              |

### Error/Qscore Parameters

Configuration for error model:

| Parameter          | Description                                                   | format   | Default Value                                 |
|--------------------|---------------------------------------------------------------|-------|--------------------------------------------|
| `trained_model`    | Badread pre-trained error/Qscore model name                   |  STR  | `nanopore2023`                                |
| `badread_identity` | Comma-separated values for Badread identity parameters        |  STR  | `"98,2,99"`                                   |
| `error_model`      | Custom error model file (optional)                            |  .TXT  | `null`                                        |
| `qscore_model`     | Custom Q-score model file (optional)                          |  .TXT  | `null`                                        |
| `build_model`      | to build your own error/Qscor model                           |  STR  | `false`                                       |
| `fastq_model`      | reference real read (.fastq) to train error model   (optional) |   FASTQ      | `false`                                       |


### Additional Parameters

| Parameter          | Description                                            |   Format    | Default Value                                 |
|--------------------|--------------------------------------------------------|-------|-----------------------------------------------|
| `amp`              | Amplification factor                                   |   INT    | `1`                                           |
| `outdir`           | Output directory for results                           |   PATH    | `"results"`                                   |
| `projectName`      | Name of the project                                    |    STR   | `"test_project"`                              |

### Run Parameters

Configuration for running the workflow:

| Parameter         | Description                        |   Format    | Default Value             |
|-------------------|------------------------------------|-------------|---------------------------|
| `threads`         | Number of threads to use           |      INT       | `4`                       |
| `container`       | Docker container for the workflow  |     STR        | `'hamraouii/asarusim:0.1'`    |
| `docker.runOptions` | Docker run options to use       |    STR         | `'-u $(id -u):$(id -g)'`  |

For more details about workflow options see the [Input parameters](https://genomiqueens.github.io/AsaruSim/parameters/) section in the documentation.

### File format discription
#### `--bc_counts`
To simulate specific UMI counts per cell barcode with random transcripts, set the --bc_counts parameter to the path of a UMI counts .CSV file. This parameter eliminates the need for an input matrix, enabling the simulation of UMI counts where transcripts are chosed randomly.

example of UMI counts per CB file:
|CB 	|counts|
|--------------|------|
|ACGGCGATCGCGAGCC 	|1260|
|ACGGCGATCGCGAGCC 	|1104|

#### `--cell_types_annotation`
AsaruSim allows user to estimate this characteristic from an existing count table. To do so, the user need to set --sim_celltypes parameter to true and to provide the list of cell barcodes of each group (.CSV file) using --cell_types_annotation parameter:
|CB |	cell_type|
|--------------|------|
|ACGGCGATCGCGAGCC| 	type 1|
|ACGGCGATCGCGAGCC|	type 2|

AsaruSim will then use the provided matrix to estimate characteristic of each cell groups and generate a synthetic count matrix.

## $\textcolor{#FF7F00}{Usage}$
User can choose among 4 ways to simulate template reads.
- use a real count matrix
- estimated the parameter from a real count matrix to simulate synthetic count matrix 
- specified by his/her own the input parameter
- a combination of the above options

We use SPARSIM tools to simulate count matrix. for more information a bout synthetic count matrix, please read [SPARSIM](https://gitlab.com/sysbiobig/sparsim/-/blob/master/vignettes/sparsim.Rmd?ref_type=heads#Sec_Input_parameter_estimated_from_data) documentaion.

### EXAMPLES 

<img width="602" alt="Capture d’écran 2025-03-03 à 23 44 57" src="https://github.com/user-attachments/assets/4a8b97c9-ead7-4216-9a84-b36369b55702" />

##### Sample data
A demonstration dataset to initiate this workflow is accessible on zenodo DOI : [10.5281/zenodo.12731408](https://zenodo.org/records/12731409). This dataset is a subsample from a Nanopore run of the [10X 5k human pbmcs](https://www.10xgenomics.com/datasets/5k-human-pbmcs-3-v3-1-chromium-controller-3-1-standard).

The human GRCh38 [reference transcriptome](https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/), [gtf annotation](https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/) and [fasta referance genome](https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/) can be downloaded from Ensembl.


You can use the `run_test.sh` script to automatically download all required datasets.

##### BASIC WORKFLOW

```bash
 nextflow run main.nf --matrix dataset/sub_pbmc_matrice.csv \
                      --transcriptome dataset/Homo_sapiens.GRCh38.cdna.all.fa \
                      --features gene_name \
                      --gtf dataset/GRCh38-2020-A-genes.gtf
```

##### WITH PCR AMPLIFICTION

```bash
 nextflow run main.nf --matrix dataset/sub_pbmc_matrice.csv \
                      --transcriptome dataset/Homo_sapiens.GRCh38.cdna.all.fa \
                      --features gene_name \
                      --gtf dataset/GRCh38-2020-A-genes.gtf \
                      --pcr_cycles 2 \
                      --pcr_dup_rate 0.7 \
                      --pcr_error_rate 0.00003
```

##### WITH SIMULATED CELL TYPE COUNTS

```bash
 nextflow run main.nf --matrix dataset/sub_pbmc_matrice.csv \
                      --transcriptome dataset/Homo_sapiens.GRCh38.cdna.all.fa \
                      --features gene_name \
                      --gtf dataset/GRCh38-2020-A-genes.gtf \
                      --sim_celltypes true \
                      --cell_types_annotation dataset/sub_pbmc_cell_type.csv
```

##### USING A SPARSIM PRESET MATRIX (e.g Chu et al. 10X Genomics datasets)

```bash
nextflow run main.nf --matrix Chu_param_preset \
                      --transcriptome datasets/Homo_sapiens.GRCh38.cdna.all.fa \
                      --features gene_name \
                      --gtf datasets/Homo_sapiens.GRCh38.112.gtf
```

##### WITH PERSONALIZED ERROR MODEL

```bash
nextflow run main.nf --matrix dataset/sub_pbmc_matrice.csv \
                     --transcriptome dataset/Homo_sapiens.GRCh38.cdna.all.fa \
                     --features gene_name \
                     --gtf dataset/GRCh38-2020-A-genes.gtf \
                     --build_model true \
                     --fastq_model dataset/sub_pbmc_reads.fq \
                     --ref_genome dataset/GRCh38-2020-A-genome.fa 
```

##### COMPLETE WORKFLOW

```bash
 nextflow run main.nf --matrix dataset/sub_pbmc_matrice.csv \
                      --transcriptome dataset/Homo_sapiens.GRCh38.cdna.all.fa \
                      --features gene_name \
                      --gtf dataset/GRCh38-2020-A-genes.gtf \
                      --sim_celltypes true \
                      --cell_types_annotation dataset/sub_pbmc_cell_type.csv \
                      --build_model true \
                      --fastq_model dataset/sub_pbmc_reads.fq \
                      --ref_genome dataset/GRCh38-2020-A-genome.fa \
                      --pcr_cycles 2 \
                      --pcr_dup_rate 0.7 \
                      --pcr_error_rate 0.00003
```

## $\textcolor{#FF7F00}{Output}$

After execution, results will be available in the specified `--outdir`. This includes simulated Nanopore reads `simulated.fastq.gz`, along with log file and QC report.

```bash
QC_report.html                    # final QC report
pipeline_info                     # Pipeline execution trace, timeline and Dag
simulated.fastq.gz                # Simulated reads including sequencing errors
template.fa.gz                    # Simulated template
```

#### Cleaning Up

To clean up temporary files generated by Nextflow:

```bash
nextflow clean -f
```

## $\textcolor{#FF7F00}{Workflow}$

![Workflow Schema](images/workflow.png)

## $\textcolor{#FF7F00}{Acknowledgements}$

- We would like to express our gratitude to [Youyupei](https://github.com/youyupei) for the development of [SLSim](https://github.com/youyupei/SLSim), which has been helpful to the `AsaruSim` workflow.
- Additionally, our thanks go to the teams behind [Badread](https://github.com/rrwick/Badread), [SPARSim](https://gitlab.com/sysbiobig/sparsim) and [Trans-NanoSim](https://github.com/bcgsc/NanoSim) whose tools are integral to the `AsaruSim` workflow.

## $\textcolor{#FF7F00}{Support\ and\ Contributions}$

For support, please open an issue in the repository's "Issues" section. Contributions via Pull Requests are welcome. Follow the contribution guidelines specified in `CONTRIBUTING.md`.

## $\textcolor{#FF7F00}{License}$

`AsaruSim` is distributed under a specific license. Check the `LICENSE` file in the GitHub repository for details.

## $\textcolor{#FF7F00}{Citation}$

If you use AsaruSim in your research, please cite this manuscript:<br>
> Ali Hamraoui, Laurent Jourdren, Morgane Thomas-Chollier, AsaruSim: a single-cell and spatial RNA-Seq Nanopore long-reads simulation workflow, Bioinformatics, 2025;, btaf087, https://doi.org/10.1093/bioinformatics/btaf087
