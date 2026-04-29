# BlueSTARR Variant Analysis Workflow Handbook

This document describes the computational workflow used to generate and process BlueSTARR regulatory variant predictions. It serves as a technical reference for pipeline structure, inputs, outputs, and intermediate artifacts.

The workflow proceeds in sequential stages. Each stage produces files required by later stages.

---

# Stage 01 — Region Preparation

Script:

```
workflow/workers/01_prepare_regions.py
```

## Purpose

Prepare genomic regions for BlueSTARR prediction by:

* filtering regions based on size and TSS distance
* removing sex chromosomes
* generating padded prediction intervals
* writing region metadata for downstream annotation
* writing prediction intervals for sequence extraction

---

## Inputs

### Primary input

Configured in:

```
prediction_generation.input_file
```

Expected columns are defined in:

```
prediction_generation.columns
```

Required semantic fields:

| Field        | Description             |
| ------------ | ----------------------- |
| chrom        | Chromosome name         |
| start        | Region start coordinate |
| end          | Region end coordinate   |
| tss_distance | Distance to nearest TSS |

Optional:

| Field | Description         |
| ----- | ------------------- |
| gene  | Nearest gene symbol |

Input may be:

* headerless (BED-like)
* headered (tabular)

Column interpretation is controlled entirely by config.

---

## Filtering rules

Regions are removed if:

* region length < `min_size`
* midpoint distance from TSS < `min_distance`
* chromosome ends in `X` or `Y`

Midpoint distance is defined as:

```
abs(tss_distance) + (region_length / 2)
```

---

## Interval construction

Prediction intervals are constructed as follows:

If:

```
region_length ≥ region_size
```

then:

```
prediction_region =
(midpoint ± region_size/2) padded by 150 bp (left) and 149 bp (right)
```

Otherwise:

```
prediction_region =
(original region) padded by 150 bp (left) and 149 bp (right)
```

Padding ensures compatibility with the 300 bp BlueSTARR window.

---

## Outputs

Written to:

```
workflow.prediction_generation.intermediate_dir
```

### File 1

```
region_keys.txt
```

Columns:

| Column            | Description                |
| ----------------- | -------------------------- |
| source_region     | Original region interval   |
| prediction_region | Padded prediction interval |
| gene              | Optional nearest gene      |
| tss_distance      | Midpoint distance to TSS   |

Example:

```
chr1:10000-10200   chr1:9850-10349   MYC   5421
```

---

### File 2

```
regions_all.txt
```

Contains:

```
prediction_region
```

only, one per line, no header.

Example:

```
chr1:9850-10349
chr1:20400-20699
chr2:99800-100099
```

This file is used for sequence extraction and chunk generation.

---

## Configuration parameters

Controlled by:

```
prediction_generation
```

Relevant fields:

| Parameter        | Description                         |
| ---------------- | ----------------------------------- |
| input_file       | Region input file                   |
| has_header       | Whether input contains column names |
| columns          | Column role specification           |
| region_size      | Target prediction window            |
| min_size         | Minimum region size                 |
| min_distance     | Minimum midpoint TSS distance       |
| region_keys_file | Output metadata filename            |
| regions_all_file | Output interval list filename       |

---

## Downstream dependencies

Outputs are consumed by:

```
02_split_regions.py
```

which:

* retrieves genomic sequence
* splits intervals into chunk files
* prepares prediction jobs

---
