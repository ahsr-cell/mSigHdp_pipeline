#! /usr/bin/env/bash

#! /bin/bash
#BSUB -J mSigHdp_pipeline_test1
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -u ar39@sanger.ac.uk
#BSUB -q basement
#BSUB -n 20
#BSUB -M10000
#BSUB -R "select[mem>10000] rusage[mem=10000]"

module load cellgen/nextflow/24.10.2
module load ISG/singularity/3.11.4

config_file = /lustre/scratch125/casm/teams/team267/projects/Pipelines/mSigHdp_pipeline/conf/sanger_lsf.config
main_script = /lustre/scratch125/casm/teams/team267/projects/Pipelines/mSigHdp_pipeline/main.nf
mutational_matrix = /lustre/scratch125/casm/teams/team267/users/ar39/0_Projects/0_mSigHdp/0_memorytests/0_data/1_Testing/EAC_GEC_Manuscript_Subset_v1.SBS96.txt
sample_matrix = /lustre/scratch125/casm/teams/team267/users/ar39/0_Projects/0_mSigHdp/0_memorytests/0_data/1_Testing/EAC_Manuscript_sample_key.csv
mutational_context = SBS96
outdir = /lustre/scratch125/casm/teams/team267/projects/Pipelines/1_tests/

nextflow run ${main_script} \
     --sample_matrix ${sample_matrix} \
     --mutational_context ${mutational_context} \
     --analysis_type analysis \
     --burnin_iterations 100 \
     --burnin_multiplier 10 \
     --posterior 10 \
     --posterior_iterations 10 \
     --mutational_matrix ${mutational_matrix} \
     --outdir ${outdir} \
     -resume \
     -profile singularity \
     -c ${config_file}
