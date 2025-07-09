#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { mSigHdp } from './workflows/mSigHdp.nf'
include { SigProfilerPlotting } from './workflows/SigProfilerPlotting.nf'
include { SigProfilerAssignment } from './workflows/SigProfilerAssignment.nf'

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
                params.mutational_matrix, 
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin.iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations
            )

            SigProfilerPlotting(
                mSigHdp.out,
                params.mutational_context,
                
            )
            SigProfilerAssignment(
                mSigHdp.out,
                params.mutational_context,
            )
        }
        else {
            mSigHdp(
                params.mutational_matrix, 
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin.iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations
            )
            SigProfilerPlotting(
                mSigHdp.out,
                params.mutational_context
            )
        }
    }
    else {
        if (params.decompose == true) {
            mSigHdp(
                params.mutational_matrix, 
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin.iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations
            )
            SigProfilerAssignment(
                mSigHdp.out,
                params.mutational_context
            )
        }
        mSigHdp(
                params.mutational_matrix, 
                params.sample_matrix, 
                params.mutational_context, 
                params.analysis_type,
                params.burnin.iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations
            )
    } 
}