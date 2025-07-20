#!/bin/bash -ue
Rscript --vanilla /lustre/scratch125/casm/teams/team267/projects/Pipelines/mSigHdp_pipeline/bin/mSigHdp.R -s EAC_Manuscript_sample_key.csv -c SBS96 -a analysis -b 100 -x 10 -o 10 -i 10 EAC_GEC_Manuscript_Subset_v1.SBS96.txt
