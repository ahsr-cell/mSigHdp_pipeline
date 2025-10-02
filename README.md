## Introduction

![mSigHdp pipeline overview](/assets/images/mSigHdp_overview.png)

**mSigHdp pipeline** is a bioinformatics pipeline for standardised mutational signature extraction using [mSigHdp](https://github.com/steverozen/mSigHdp).

There are three main processes of the pipeline: mSigHdp, SigProfilerPlotting, and SigProfilerAssignment. 

The pipeline executes all three processes (by default). mSigHdp runs first, and the results it generates (i.e., a matrix containing the *de novo* signatures) are subsequently fed to [SigProfilerPlotting](https://github.com/AlexandrovLab/SigProfilerPlotting) for signature spectra plotting and [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment) for a preliminary decomposition. Depending on user needs, SigProfilerPlotting and SigProfilerAssignment can be turned off. The pipeline is designed to be compatible with the following mutation type classifications: SBS96, SBS288, SBS1536, DBS78, ID83. 


### mSigHdp
mSigHdp uses hierarchical Dirichlet processes to identify mutational signatures present within samples. 

The primary input of mSigHdp is a `mutational matrix` (specified via /path/to/mutation_matrix), with an expected format of one row per mutation type (e.g., the 96 SBS, A[C>A]A) and one column per sample. 

mSigHdp has two run modes (`analysis` and `testing`) set by `analysis_type`. 

**Testing** is intended for initial runs (aka "is this working" scenarios). It is run with minimal settings (1 Gibbs sampling chain using one thread, running 100 burn-in iterations, collecting 5 posterior samples off of each chain with 5 iterations between each, to allow for a short execution time. As this is for testing purposes, these run settings cannot be changed. 

**Analysis** is used for full-analysis runs, utilising 20 Gibbs sampling chains across 20 threads. These chains run 50,000 burn-in iterations (5,000 burn-in iterations with a 10x multiplier) and collect 250 posterior samples from each chain, with 100 iterations collected between each sample. These are standardised settings, optimised and conducted in the [Cancer Grand Challenges Mutographs project](https://www.cancergrandchallenges.org/mutographs). Users can change these settings by specifying `--burnin_iterations`, `--burnin_multiplier`, `--posterior`, and `--posterior_iterations`.

mSigHdp can be run with **hierarchy** (setting `hierarchy = true`, specifying the hierarchy table via `--hierarchy_matrix /path/to/hierarchy_matrix`, and specifying the hierarchy parameter by `--hierarchy_parameter`) or **no hierarchy** (setting `hierarchy=false`).

The primary output of mSigHdp, found under the subdirectory `/deNovo_signatures/`) include the extracted *de novo* signature table (`mSigHdp_deNovoSignatures.txt`) and, if detected, the low confidence signature table (`mSigHdp_lowConfSignatures.txt`). Standard QC plots and matrices are also included within `/deNovo_signatures/`.

### SigProfilerPlotting
SigProfilerPlotting is used to plot the extracted *de novo* signature spectra. The pipeline will feed the extracted *de novo* signature table (i.e., `mSigHdp_deNovoSignatures.txt`) and the user-specified mutational context (done via `--mutational_context SBS96`, default is SBS96) to SigProfilerPlotting for plotting. The output of SigProfilerPlotting can be found under the directory `/Signature_Spectra/`, further broken down to subdirectories `/DeNovoSignatures/` and, if detected, `/LowConfidenceSignatures/`, both containing PDFs of their plots. 

Setting `plotting` to `false` will turn this functionality off. 

### SigProfilerAssignment
SigProfilerAssignment is used to decompose the extracted *de novo* signatures. The pipeline feeds the extracted *de novo* signature table and the previous, user-provided mutational matrix (specified via `mutational matrix`) to SigProfilerAssignment, which executes its `decompose_fit()` function. The output of SigProfilerAssignment is located under the directory `/SigProfilerDecomposition/`, containing subdirectories `/Activities/`, `/Signatures/`, and `/Solution_Stats/` and decomposition plots and mappings to COSMIC signatures (e.g., [COSMIC SBS](https://cancer.sanger.ac.uk/signatures/sbs/)).

Setting `decompose` to `false` will turn this functionality off. 

## Dependencies
* Nextflow >= 24.04.2 required
* Python, required packages: [SigProfilerPlotting](https://github.com/AlexandrovLab/SigProfilerPlotting), [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment), [argparse](https://docs.python.org/3/library/argparse.html), []()
* R, required packages: [mSigHdp](https://github.com/steverozen/mSigHdp), [hdpx](https://github.com/steverozen/hdpx), [tidyverse](https://www.tidyverse.org/), [argparse](https://cran.r-project.org/web/packages/argparse/index.html)

## Installation
Clone this repository via
 > git clone git@github.com:ahsr-cell/mSigHdp_pipeline.git

## Usage

### Inputs

| Input      | Description |
| ----------- | ----------- |
| `mutation_matrix`      | Required input file, provided as a tab-delimited file (.tsv). The expected format is a matrix with one row per mutation type and one column per sample. Include the mutation types under a column labelled as 'MutationType'. It is highly recommended to generate mutation matrices via [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator). Please see [example_mutation_matrix.tsv](https://github.com/ahsr-cell/mSigHdp_pipeline/blob/main/docs/example_input_data/example_mutation_matrix.tsv) for an example.        |
| `hierarchy`   | Required value, provided as a string. Options are `true` or `false`         |
| `hierarchy_matrix`   | Optional input file, provided if `hierarchy = true`. The expected format is a matrix with one column per sample ID (matching the input mutation matrix) and one column specifying hierarchy groupings. See [example_hierarchy_matrix.tsv](https://github.com/ahsr-cell/mSigHdp_pipeline/blob/main/docs/example_input_data/example_hierarchy_matrix.tsv) for an example.         |
| `hierarchy_parameter`   | Optional value, provided as a string if `hierarchy = true`. This should be formatted exactly as the column name specifying hierarchy groupings in the input hierarchy_matrix. E.g., if a user provided [example_hierarchy_matrix.tsv](https://github.com/ahsr-cell/mSigHdp_pipeline/blob/main/docs/example_input_data/example_hierarchy_matrix.tsv), `hierarchy_parameter` would be set as `hierarchy_parameter=Grouping`             |
| `analysis_type`   | Required value, provided as a string. Options are `analysis` or `testing`, default is `analysis`         |
| `mutational_context`   | Required value, provided as a string. Options are `SBS96`, `SBS288`, `SBS1536`, `DBS78`, `ID83`, default is `SBS96`         |
| `plotting`   | Required value, provided as a string. Options are `true` or `false`, default is `true`         |
| `decompose`   | Required value, provided as a string. Options are `true` or `false`, default is `true`         |
| `outdir`   | Required path, specifying the location of output files generated by pipeline. See [`output.md`](docs/output.md) for the files/directories that will be contained within this directory.        |
| `burnin_iterations`[^*]   | Optional value, provided as a double. Used to change number of burn-in iterations conducted. Default is `5000`          |
| `burnin_multiplier`[^*]   | Optional value, provided as a double. Used to change multiplier of burn-in iterations conducted. Default is `10`         |
| `posterior`[^*]   | Optional value, provided as a double. Used to change number of posterior samples. Default is `250`         |
| `posterior_iterations`[^*]   | Optional value, provided as a double. Used to change number of iterations conducted between each posterior sample. Default is `100`         |
| `concentration_parameter`[^*]   | Optional value, provided as a double. Used to change number of iterations concentration parameter sampling. Default is `3`         |
| `chains`[^*]   | Optional value, provided as a double. Used to change number of posterior sampling chains and number of CPUs to use (for parallelisation). Default is `20`         |
| `clusters`[^*]   | Optional value, provided as a double. Used to change suggested initial number of raw clusters. Default is `16`         |
| `alpha`[^*]   |  Optional value, provided as a double. Used to change gamma.alpha, the shape parameter of the gamma distribution prior. Default is `1`         |
| `beta`[^*] [^**]   | Optional value, provided as a double. Used to change gamma.beta, the inverse scale (rate) parameter of the gamma distribution prior. Default is `20`         |
| `confidence`[^*]   | Optional value, provided as a double. Used to change `high.confidence.prop`, threshold to determine signatures with high confidence. Default is `0.9`         |
[^*]: mSigHdp run options that can be changed in analysis runs. Recommended to change only after conducting an initial analysis run with default settings and having validated reasonable cause for changing number of iterations.
[^**]: 20 recommended for SBS signature extraction, while 50 is recommended for DBS and ID signature extraction. See [mSigHdp manual](https://github.com/steverozen/mSigHdp/blob/v2.1.2-branch/mSigHdp_2.1.2.pdf) 

The pipeline can be run using:

```
nextflow run /path/to/mSigHdp_pipeline/main.nf \
     -profile <docker/singularity/.../institute> \
     -c /path/to/config_file \
     --mutational_matrix /path/to/mutation_matrix.tsv \
     --mutational_context <SBS96/SBS288/SBS1536/DBS78/ID83> \
     --analysis_type <analysis/testing> \
     --burnin_iterations <desired_number_of_burnin_iterations> \
     --burnin_multiplier <desired_burnin_multiplier> \
     --posterior <desired_number_of_posterior_samples> \
     --posterior_iterations <desired_number_of_posterior_iterations> \
     --chains <desired_number_of_posterior_sampling_chains> \
     --clusters <desired_number_of_initial_HDP_clusters> \
     --alpha <desired_gamma.alpha_value> \
     --beta <desired_gamma.beta_value> \
     --confidence <desired_confidence_threshold> \
     --outdir /path/to/outdir/ \
     --plotting <true/false> \
     --decompose <true/false> \
     --hierarchy <true/false> \
     --hierarchy_matrix /path/to/hierarchy_key.tsv \
     --hierarchy_parameter <column-name_of_hierarchyparameter_hierarchy_key> 
```

### Sanger Users
Sangers can run the pipeline using the following wrapper script:

```
#! /usr/bin/env bash

#BSUB -J USER_JOB_NAME
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -u USER@sanger.ac.uk
#BSUB -q week #set according time requirements for largest local running job
#BSUB -n 1 #Single CPU node used for pipeline task submitter
#BSUB -M5000 #set according memory requirements for for largest local running job
#BSUB -R "select[mem>5000] rusage[mem=5000]" #set according to sample type and memory requirements

# Run inputs and options, change accordingly
mutational_matrix=/path/to/mutation_matrix.tsv

mutational_context=SBS96 #options include SBS96, SBS288, SBS1536, DBS78, ID83 
analysis_type=analysis #options include testing/analysis
outdir=/path/to/desired/output_directory #set to location of output files
plotting=true #set to false if just signature extraction
decompose=true #set to false if just signature extraction

hierarchy=false #if running hierarchy, set to true
hierarchy_matrix=/path/to/input/hierarchy_matrix.tsv #if running no hierarchy, delete/comment this line out 
hierarchy_parameter=sample_type #if running no hierarchy delete/comment this line out

#mSigHdp analysis run options. If desired, change accordingly to your needs
burnin_iterations=5000 
burnin_multiplier=10
posterior=250
posterior_iterations=100
concentration_parameter=3
chains=20
clusters=20
alpha=1
beta=20
confidence=0.9

### Necessary for the pipeline - do not change
module load cellgen/nextflow/25.04.4
module load ISG/singularity/3.11.4
config_file=/lustre/scratch125/casm/teams/team267/projects/Pipelines/mSigHdp_pipeline/conf/sanger_lsf.config
main_script=/lustre/scratch125/casm/teams/team267/projects/Pipelines/mSigHdp_pipeline/main.nf
profile=singularity

nextflow run ${main_script} \
     -profile ${profile} \
     -c ${config_file} \
     --mutational_matrix ${mutational_matrix} \
     --mutational_context ${mutational_context} \
     --analysis_type ${analysis_type} \
     --burnin_iterations ${burnin_iterations} \
     --burnin_multiplier ${burnin_multiplier} \
     --posterior ${posterior} \
     --posterior_iterations ${posterior_iterations} \
     --chains ${chains} \
     --clusters ${clusters} \
     --alpha ${alpha} \
     --beta ${beta} \
     --confidence ${confidence} \
     --outdir ${outdir} \
     --plotting ${plotting} \
     --decompose ${decompose} \
     --hierarchy ${hierarchy} \
     --hierarchy_matrix ${hierarchy_matrix} \
     --hierarchy_parameter ${hierarchy_parameter} 
```

## Pipeline output

For more details about the output files and reports, please refer to the [`output.md`](docs/output.md) output documentation.

## Credits

The mSigHdp pipeline was written by Andrew Ramos and Phuong Le. 

We thank the following people and teams for their assistance in the development of this pipeline:

Sarah Moody, Mimy Pham, Yichen Wang, and CASM IT

## Contributions and Support

Please feel free to contribute by either creating a pull request or create a new issue on this GitHub repo.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
