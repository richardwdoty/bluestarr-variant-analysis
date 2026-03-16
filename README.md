# BlueSTARR Variant Analysis

Workflow and analysis code for generating and analyzing regulatory variant predictions using **BlueSTARR**.

This repository contains the pipeline used to generate allele-specific regulatory predictions from genomic regions and perform downstream analyses of predicted variant effects. The workflow produces per-allele regulatory activity predictions, integrates population and genomic annotations, and supports analyses of observed versus unobserved variants in human populations.

The repository is organized around three main components:

* **Prediction generation** – running BlueSTARR across genomic regions with systematic single-nucleotide perturbations
* **Prediction processing** – aggregating predictions, annotating variants, and producing structured datasets
* **Analysis** – statistical analyses of predicted variant effects

The pipeline is designed for **high-performance computing environments** and assumes access to a Slurm scheduler for large-scale prediction generation.

---

## Repository Structure

```
bluestarr-variant-analysis/
├── config/        Configuration files
├── workflow/      Prediction generation and processing workflow
│   ├── workers/   Individual task scripts
│   ├── runners/   Scripts orchestrating workflow stages
│   └── slurm/     Slurm job templates and submission scripts
├── analysis/      Downstream analysis scripts
├── refs/          Reference metadata and resources
├── bin/           External utilities required by the workflow
├── examples/      Example inputs and configuration
└── docs/          Workflow and data documentation
```

---

## Workflow Overview

### 1. Prediction Generation

Genomic regions of interest are processed to generate allele-specific BlueSTARR predictions.

Steps include:

1. Parse and filter genomic regions
2. Retrieve genomic sequences
3. Generate systematic single-nucleotide perturbations
4. Run BlueSTARR predictions for each allele at each position

Predictions are produced as **log(RNA/DNA)** regulatory activity values.

Because each region is evaluated across multiple sequence windows and alleles, prediction generation is executed as a **large Slurm array job**.

---

### 2. Prediction Processing

Raw prediction files are aggregated and annotated.

Processing includes:

* Combining prediction outputs by chromosome
* Annotating variants using external datasets (e.g., gnomAD)
* Adding gene and TSS distance annotations
* Converting long-format predictions into condensed allele-level tables
* Filtering regions not used during BlueSTARR training
* Filtering by regulatory annotations when appropriate

The final processed dataset contains, for each genomic site:

* predicted regulatory activity for all four alleles
* allele type annotations (reference, observed SNV, unobserved)
* allele frequency information
* gene and TSS distance annotations

---

### 3. Analysis

Downstream analyses evaluate patterns in predicted regulatory effects.

Analyses include:

* rank-based comparisons of observed vs unobserved variants
* statistical testing of allele effect distributions
* mutation-class–controlled analyses
* stratified analyses across genomic and population annotations

These analyses assess whether variants observed in human populations preferentially occupy particular regulatory effect configurations relative to unobserved variants.

---

## Configuration

Workflow behavior is controlled through configuration files in the `config/` directory.

Paths to input data, reference datasets, and workflow outputs are specified in the configuration file. This allows the pipeline to run in different computing environments without modifying scripts directly.

---

## Requirements

This workflow assumes access to:

* Python 3
* BlueSTARR model checkpoints
* genomic reference data
* population variant data (e.g., gnomAD)
* Slurm workload manager (recommended for prediction generation)

Additional utilities required by the workflow may be stored in the `bin/` directory.

---

## Running the Workflow

Typical workflow stages are executed through runner scripts located in:

```
workflow/runners/
```

These scripts coordinate worker tasks and may generate or submit Slurm job arrays for large-scale computation.

Detailed usage instructions will be provided in the `docs/` directory.

---

## Data Storage

The workflow produces large intermediate datasets. These outputs are written to configurable run directories specified in the configuration file and are **not stored within the repository**.

---

## Citation

If you use this code, please cite the corresponding dissertation or publication describing the BlueSTARR analyses.

---

## License

This project is licensed under the MIT License.