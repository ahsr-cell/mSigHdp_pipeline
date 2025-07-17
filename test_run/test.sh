#! /usr/bin/env/bash

 nextflow run /Users/ar39/Documents/mSigHDP_HDP/7_Pipelines/0_mSigHdp/mSigHdp_pipeline/main.nf \
     --sample_matrix /Users/ar39/Documents/mSigHDP_HDP/0_Background/2_Example_Data/1_Test/EAC_Manuscript_sample_key.csv \
     --mutational_context SBS96 \
     --analysis_type testing \
     --burnin_iterations 100 \
     --burnin_multiplier 10 \
     --posterior 10 \
     --posterior_iterations 10 \
     --mutational_matrix /Users/ar39/Documents/mSigHDP_HDP/0_Background/2_Example_Data/1_Test/EAC_GEC_Manuscript_Subset_v1.SBS96.txt \
     --outdir /Users/ar39/Documents/mSigHDP_HDP/7_Pipelines/0_mSigHdp/0_test/ -resume