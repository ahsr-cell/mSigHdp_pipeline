#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { mSigHdp_hierarchy } from './modules/mSigHdp_hierarchy.nf'
include { mSigHdp_flat } from './modules/mSigHdp_flat.nf'
include { SigProfilerPlotting as SigPlt_Extracted } from './modules/SigProfilerPlotting.nf'
include { SigProfilerPlotting as SigPlt_LowConfidence } from './modules/SigProfilerPlotting.nf'
include { SigProfilerAssignment as SigPA_Extracted } from './modules/SigProfilerAssignment.nf'
//include { SigProfilerAssignment as SigPA_LowConfidence } from './workflows/SigProfilerAssignment.nf'

//
// WORKFLOW: Run main analysis pipeline depending on user inputs
//
////def file = new File(params.mutational_matrix)
//file.eachLine { line ->
//    def values = line.split(",")  // Use ',' as delimiter
//}
//println values
//
////def matrix = file.readLines().collect { line ->
//    line.split(',') as List
////}

// Print the matrix
////matrix.each { row ->
////    println row
////}

workflow {

    main:
    //
    // WORKFLOW: Full suite of analysis: mSigHdp, SigProfilerPlotting, and SigProfilerAssignment
    //
    if (params.hierarchy == true) {
        mSigHdp_hierarchy(
                params.hierarchy_matrix,
                params.hierarchy_parameter,
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )
            if (params.plotting == true) {
                if (params.decompose == true) {
                    SigPlt_Extracted(
                    mSigHdp_hierarchy.out.deNovo_extractedsigs,
                    params.mutational_context,
                    sig_type = "DeNovoSignatures"
                    )
                    if (mSigHdp.out.deNovo_lowconfsigs != null) {
                        SigPlt_LowConfidence(
                        mSigHdp_hierarchy.out.deNovo_lowconfsigs,
                        params.mutational_context,
                        sig_type = "LowConfidenceSignatures"
                        )
                    }
                    SigPA_Extracted(
                        mSigHdp_hierarchy.out.deNovo_extsigs_sigPA,
                        params.mutational_matrix
                    )
                } else {
                    SigPlt_Extracted(
                    mSigHdp_hierarchy.out.deNovo_extractedsigs,
                    params.mutational_context,
                    sig_type = "DeNovoSignatures"
                    )
                    if (mSigHdp_hierarchy.out.deNovo_lowconfsigs != null) {
                        SigPlt_LowConfidence(
                            mSigHdp_hierarchy.out.deNovo_lowconfsigs,
                            params.mutational_context,
                            sig_type = "LowConfidenceSignatures"
                        )
                    }
                }
            } else if (params.decompose == true) {
                SigPA_Extracted(
                        mSigHdp_hierarchy.out.deNovo_extsigs_sigPA,
                        params.mutational_matrix
                    )
            }
    } else {
        mSigHdp_flat(
                params.mutational_context, 
                params.analysis_type,
                params.burnin_iterations,
                params.burnin_multiplier,
                params.posterior,
                params.posterior_iterations,
                params.mutational_matrix
            )
            if (params.plotting == true) {
                if (params.decompose == true) {
                    SigPlt_Extracted(
                    mSigHdp_flat.out.deNovo_extractedsigs,
                    params.mutational_context,
                    sig_type = "DeNovoSignatures"
                    )
                    if (mSigHdp_flat.out.deNovo_lowconfsigs != null) {
                        SigPlt_LowConfidence(
                        mSigHdp_flat.out.deNovo_lowconfsigs,
                        params.mutational_context,
                        sig_type = "LowConfidenceSignatures"
                        )
                    }
                    SigPA_Extracted(
                        mSigHdp_flat.out.deNovo_extsigs_sigPA,
                        params.mutational_matrix
                    )
                } else {
                    SigPlt_Extracted(
                    mSigHdp_flat.out.deNovo_extractedsigs,
                    params.mutational_context,
                    sig_type = "DeNovoSignatures"
                    )
                    if (mSigHdp_flat.out.deNovo_lowconfsigs != null) {
                        SigPlt_LowConfidence(
                            mSigHdp_flat.out.deNovo_lowconfsigs,
                            params.mutational_context,
                            sig_type = "LowConfidenceSignatures"
                        )
                    }
                }
            } else if (params.decompose == true) {
                SigPA_Extracted(
                        mSigHdp_flat.out.deNovo_extsigs_sigPA,
                        params.mutational_matrix
                    )
            }
    }
}