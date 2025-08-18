# mSigHdp_pipeline: Output

## Introduction

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.


## Pipeline output overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [deNovo_signatures](#deNovo_signatures) - mSigHdp run
- [Signature_Spectra](#Signature_Spectra) - SigProfilerPlotting spectra
- [SigProfilerDecomposition](#SigProfilerDecomposition) - SigProfilerAssignment decomposition
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### mSigHdp output
- `deNovo_signatures/`
  - `diagnostics.components.in.which.gibbs.samples.pdf`
  - `diagnostics.likelihood.pdf`
  - `diagnostics.numcluster.pdf`
  - `diagnostics.signatures.pdf`
  - `extracted.signatures.csv`
  - `extracted.signatures.post.samp.number.csv`
  - `hdp.retval.Rdata`
  - `inferred.exposure.count.pdf`
  - `inferred.exposure.proportion.pdf` 
  - `inferred.exposures.csv`
  - `low.confidence.signatures.csv`
  - `low.confidence.signatures.post.samp.number.csv`
  - `mSigHdp_deNovoSignatures.txt`
  - `mSigHdp_deNovoSigs_sigPADecomp.txt`
  - `mSigHdp_lowConfSignatures.txt`

### SigProfilerPlotting output
- `Signature_Spectra/`
  - `DeNovoSignatures/`
    - `pkl/`
      - `SBS96.pkl`
    - `SBS_96_plots_deNovoSignatures.pdf`
  - `LowConfidenceSignatures/`
    - `pkl/`
      - `SBS96.pkl`
    - `SBS_96_plots_LowConfidenceSignatures.pdf`

### SigProfilerAssignment output
- `SigProfilerDecomposition/`
  - `Decompose_Solution/`
    - `Activities/`
      - `Decompose_Solution_Activities.txt`
      - `Decompose_Solution_Activity_Plots.pdf`
      - `Decompose_Solution_TMB_plot.pdf`
      - `Decomposed_MutationType_Probabilities.txt`
    - `Signatures/`
      - `Decompose_Solution_Activities.txt`
      - `SBS_96_plots_Decompose_Solution.pdf`
    - `Solution_Stats/`
      - `Cosmic_SBS96_Decomposition_Log.txt`
      - `Decompose_Solution_Samples_Stats.txt`
      - `Decompose_Solution_Signature_Assignment_log.txt`
    - `De_Novo_map_to_COSMIC_SBS96.csv`
    - `SBS96_Decomposition_Plots.pdf`     
  - `pkl/`
    - `CosmicTemplates/`
      - `COSMIC_v3.4_SBS_GRCH38.json`
    - `SBS96.pkl`
  - `JOB_METADATA_SPA.txt/`