#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { mSigHdp } from './workflows/mSigHdp.nf'
include { SigProfilerPlotting_ExtractedSigs } from './workflows/SigProfilerPlotting_ExtractedSigs.nf'
include { SigProfilerPlotting_LowConfidenceSigs } from './workflows/SigProfilerPlotting_LowConfidenceSigs.nf'
include { SigProfilerAssignment_ExtractedSigs } from './workflows/SigProfilerAssignment_ExtractedSigs.nf'
include { SigProfilerAssignment_LowConfidenceSigs } from './workflows/SigProfilerAssignment_LowConfidenceSigs.nf'

//
// WORKFLOW: Run main analysis pipeline depending on user inputs
//

workflow {

    main:

    //
    // WORKFLOW: Full suite of analysis: mSigHdp, SigProfilerPlotting, and SigProfilerAssignment
    //
    if (params.plotting == true){
        if (params.decompose == true){
            mSigHdp(
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )

            SigProfilerPlotting_ExtractedSigs(
                mSigHdp.out.deNovo_extractedsigs,
                params.mutational_context
            )

            SigProfilerAssignment_ExtractedSigs(
                mSigHdp.out.deNovo_extractedsigs,
                params.mutational_matrix
            )
            if (mSigHdp.out.deNovo_lowconfsigs != null){
                SigProfilerPlotting_LowConfidenceSigs(
                    mSigHdp.out.deNovo_lowconfsigs,
                    params.mutational_context
                )
                SigProfilerAssignment_LowConfidenceSigs(
                    mSigHdp.out.deNovo_lowconfsigs,
                    params.mutational_matrix
                )
            }
        }
        else {
            mSigHdp(
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )
            SigProfilerPlotting_ExtractedSigs(
                mSigHdp.out.deNovo_extractedsigs,
                params.mutational_context
            )
            if (mSigHdp.out.deNovo_lowconfsigs != null){
                SigProfilerPlotting_LowConfidenceSigs(
                    mSigHdp.out.deNovo_lowconfsigs,
                    params.mutational_context
                )
            }
        }
    }
    else {
        if (params.decompose == true) {
            mSigHdp(
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )
            SigProfilerAssignment_ExtractedSigs(
                mSigHdp.out.deNovo_extractedsigs,
                params.mutational_matrix
            )
            if (mSigHdp.out.deNovo_lowconfsigs != null){
                SigProfilerAssignment_LowConfidenceSigs(
                    mSigHdp.out.deNovo_lowconfsigs,
                    params.mutational_matrix
                )
            }
        }
        mSigHdp(
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )
    } 
}