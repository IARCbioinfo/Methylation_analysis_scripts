# Methylation_analysis_scripts
R scripts for Illumina array analyses

## Pre-processing: *Methylation_pre-processing.R*

This script performs Illumina EPIC 850K array pre-processing and QC from idat files. 

## Prerequisites
This R script requires the following packages:
- minfi
- ENmix
- MASS
- broom

### Usage
```bash
Rscript Methylation_pre-processing.R -f input_folder -t targetfile.txt [options]
```

| **PARAMETER** | **DEFAULT** | **DESCRIPTION** |
|-----------|--------------:|-------------| 
*-f* | . | folder with idat files |
*-o* | out | output directory name |
*-t* |  . | target file with sample information |
*-c* | NULL | file with list of cross-reactive probes |
*-s* | FALSE | filter SNPs-associated probes |
*-m* | FALSE | filter multimodal probes |
*-h*    |  | Show help message and exit|

### Details
The script involves X steps
- **Raw data QC** using ...
- **Data normalization** using ...
- **Data filtering** using ...
- **Normalized data QC** using ...

### Output
- a .RData file with X objects: 

#### QC
- QC plots with...
