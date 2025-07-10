process SigProfilerPlotting {

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input: 
    path deNovoSignatures_matrix
    val mutational_context

    output:
    path "Signature_Spectra"

    script:
    """
    python3 SigProfilerPlotting.py --mutational_context "${mutational_context}" --deNovoSignatures_matrix "${deNovoSignatures_matrix}"
    """
}