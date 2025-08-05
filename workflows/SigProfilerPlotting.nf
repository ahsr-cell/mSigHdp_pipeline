process SigProfilerPlotting {

    publishDir "${params.outdir}/Signature_Spectra/", mode: "copy"

    input: 
    path deNovoSignatures_matrix
    val mutational_context
    val sig_type

    output:
    path "${sig_type}"

    script:
    """
    export MPLCONFIGDIR=${workDir}
    rm -rf Signature_Spectra
    mkdir Signature_Spectra

    SigProfilerPlotting.py --deNovoSignatures_matrix ${deNovoSignatures_matrix} --mutational_context ${mutational_context} --output_directory "${sig_type}"
    """
}