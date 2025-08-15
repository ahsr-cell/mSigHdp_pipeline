## Introduction

![mSigHdp pipeline overview](/assets/images/mSigHdp_overview.png)

**mSigHdp pipeline** is a bioinformatics pipeline for standardised mutational signature extraction using [mSigHdp](https://github.com/steverozen/mSigHdp).

There are three main processes of the pipeline: mSigHdp, SigProfilerPlotting, and SigProfilerAssignment. 

The pipeline executes all three processes (by default). mSigHdp runs first and the results of which (i.e., a matrix containing the *de novo* signatures) are subsequently fed to [SigProfilerPlotting](https://github.com/AlexandrovLab/SigProfilerPlotting) for signature spectra plotting and [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment) for a preliminary decomposition. Depending on user needs, SigProfilerPlotting and SigProfilerAssignment can be turned off. 


### mSigHdp
mSigHdp uses hierarchical Dirichlet processes to identify mutational signatures present within samples. 

The primary input of mSigHdp is a `mutational matrix` (specified via /path/to/mutation_matrix), with an expected format of one row per mutation type (e.g., the 96 SBS, A[C>A]A) and one column per sample. 

mSigHdp has two run modes, set by `analysis_type`: `analysis` and `testing`. 

**Testing** is intended for initial runs (aka "is this working" scenarios). It is run with minimal settings (1 Gibbs sampling chain using one thread, running 100 burn-in iterations, collecting 5 posterior samples off of each chain with 5 iterations between each, to allow for a short execution time. As this run is for testing purposes, these settings cannot be changed. 

**Analysis** is used for full-analysis runs, utilising 20 Gibbs sampling chains across 20 threads. These chains run 50,000 burn-in iterations (5,000 burn-in iterations with a 10x multiplier) and collect 250 posterior samples from each chain, with 100 iterations collected between each sample. Users can change these settings by specifying `--burnin_iterations`, `--burnin_multiplier`, `--posterior`, and `--posterior_iterations`.

mSigHdp can be run with **hierarchy** (setting `hierarchy = true`, specifying the hierarchy table via `--hierarchy_matrix /path/to/hierarchy_matrix`, and specifying the hierarchy parameter by `--hierarchy_parameter`) or **no hierarchy** (setting `hierarchy=false`).

The primary output of mSigHdp, found under the subdirectory `/deNovo_signatures/`) include the extracted *de novo* signature table (`mSigHdp_deNovoSignatures.txt`) and, if detected, the low confidence signature table (`mSigHdp_lowConfSignatures.txt`). Standard QC plots and matrices are also included within `/deNovo_signatures/`.

### SigProfilerPlotting
SigProfilerPlotting is used to plot the extracted *de novo* signature spectra. The pipeline will feed the extracted *de novo* signature table (i.e., `mSigHdp_deNovoSignatures.txt`) and the user-specified mutational context (done via `--mutational_context SBS96`, default is SBS96) to SigProfilerPlotting for plotting. The output of SigProfilerPlotting can be found under the directory `/Signature_Spectra/`, further broken down to subdirectories `/DeNovoSignatures/` and, if detected, `/LowConfidenceSignatures/`, both containing PDFs of their plots. 

Setting `plotting` to `false` will turn off this functionality. 

### SigProfilerAssignment
SigProfilerAssignment is used to decompose the extracted *de novo* signatures. The pipeline feeds the extracted *de novo* signature table and the previous, user-provided mutational matrix (specified via `mutational matrix`) to SigProfilerAssignment, which executes its `decompose_fit()` function. The output of SigProfilerAssignment is located under the directory `/SigProfilerDecomposition/`, containing subdirectories `/Activities/`, `/Signatures/`, and `/Solution_Stats/` and decomposition plots and mappings to COSMIC signatures (e.g., [COSMIC SBS](https://cancer.sanger.ac.uk/signatures/sbs/)).

Setting `decompose` to `false` will turn off this functionality. 

## Dependencies
* Nextflow >= 24.04.2 required
* Python, required packages: [SigProfilerPlotting](), [SigProfilerAssignment](),[argparse](https://docs.python.org/3/library/argparse.html), []()
* R, required packages: [mSigHdp](https://github.com/steverozen/mSigHdp), [hdpx](https://github.com/steverozen/hdpx), [tidyverse](https://www.tidyverse.org/), [argparse](https://cran.r-project.org/web/packages/argparse/index.html)

## Installation
Clone this repository via

 > git clone git@github.com:ahsr-cell/mSigHdp_pipeline.git

## Usage

### Input files

| Input      | Description |
| ----------- | ----------- |
| `mutation_matrix`      | Required input file, provided as a tab-delimited file (.tsv). The expected format is a matrix with one row per mutation type and one column per sample. See /docs/example_data/example_mutation_matrix.tsv for an example.        |
| `hierarchy`   | Required value, provided as a string. Options are <`true` or `false`>         |
| `hierarchy_matrix`   | Optional input file, provided if `hierarchy = true`. Expected . See /docs/example_data/example_hierarchy_matrix.tsv for an example.         |
| `hierarchy_parameter`   | Optional value, provided as a string if `hierarchy = true`. This should be formatted exactly as the column name specifying hierarchy in the input hierarchy_matrix. E.g., if a user provided /docs/example_data/example_hierarchy_matrix.tsv, `hierarchy_parameter` would be `hierarchy_parameter=sample_type`             |
| `analysis_type`   | Required value, provided as a string. Options are <`analysis` or `testing`>, default is `analysis`         |
| `mutational_context`   | Required value, provided as a string. Options are <`SBS96`, `SBS288`, `SBS1536`, `DBS78`, `ID83`>, default is `SBS96`         |
| `plotting`   | Required value, provided as a string. Options are <`true` or `false`>, default is `true`         |
| `decompose`   | Required value, provided as a string. Options are <`true` or `false`>, default is `true`         |
| `outdir`   | Required path, specifying the location of output files generated by pipeline         |

The pipeline can be run using:

nextflow run /path/to/mSigHdp_pipeline/main.nf \
     -profile <docker/singularity/.../institute> \
     -c /path/to/config_file \
     --mutational_matrix /path/to/mutation_matrix.tsv \
     --hierarchy <true/false> \
     --hierarchy_matrix /path/to/hierarchy_key.tsv \
     --hierarchy_parameter <variable_used_for_hierarchy_as_the_column-name_in_hierarchy_key> \
     --mutational_context <SBS96/SBS288/SBS1536/DBS78/ID83> \
     --analysis_type <analysis/testing> \
     --outdir /path/to/outdir/ \
     --plotting <true/false> \
     --decompose <true/false> \


### Sanger Users
Sangers can run the pipeline using the following wrapper script. Refer to 
```
#! /usr/bin/env bash

#BSUB -J USER_JOB_NAME
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -u USER@sanger.ac.uk
#BSUB -q week #set according to sample type and time requirements
#BSUB -n 20 #if testing (e.g., analysis_type=test), set to 1; otherwise, leave at 20
#BSUB -M50000 #set according to sample type and memory requirements
#BSUB -R "select[mem>50000] rusage[mem=50000]" #set according to sample type and memory requirements

#Change these options according to run requirements
mutational_matrix=/path/to/input/mutation_matrix.tsv

hierarchy=false #if running hierarchy (single-tier only), set to true

hierarchy_matrix=/path/to/input/hierarchy_key.tsv #if running no hierarchy, delete/comment this line out 
hierarchy_parameter=sample_type #set this equal to the column name of input hierarchy_key if running no hierarchy, delete/comment this line out

mutational_context=SBS96 #For SigProfilerPlotting, options include SBS96, SBS288, SBS1536, DBS78, ID83

analysis_type=analysis #'testing' for test run, 'analysis' for full analysis run

outdir=/path/to/desired/output/ #set to location of output files

plotting=true #set to false if you do not want to plot extracted signature spectra

decompose=true #set to false if you do not want to decompose extracted signatures 

module load cellgen/nextflow/25.04.4
module load ISG/singularity/3.11.4
config_file=/lustre/scratch125/casm/teams/team267/projects/Pipelines/mSigHdp_pipeline/conf/sanger_lsf.config
main_script=/lustre/scratch125/casm/teams/team267/projects/Pipelines/mSigHdp_pipeline/main.nf
profile=singularity

nextflow run ${main_script} \
     -profile ${profile} \
     -c ${config_file} \
     --mutational_matrix ${mutational_matrix} \
     --hierarchy ${hierarchy} \
     --hierarchy_matrix ${hierarchy_matrix} \
     --hierarchy_parameter ${hierarchy_parameter} \
     --mutational_context ${mutational_context} \
     --analysis_type ${analysis_type} \
     --outdir ${outdir} \
     --plotting ${plotting} \
     --decompose ${decompose} \
```

## Pipeline output

For more details about the output files and reports, please refer to the [`output.md`](output.md) output documentation.

## Credits

The mSigHdp pipeline was written by Andrew Ramos and Phuong Le. 

We thank the following people and teams for their assistance in the development of this pipeline:

Sarah Moody, Mimy Pham, Yichen Wang, and CASM IT

## Contributions and Support

Please feel free to contribute by either creating a pull request or create a new issue on this GitHub repo.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.









