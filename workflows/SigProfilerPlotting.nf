process SigProfilerPlotting {

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input: 
    path deNovo_signatures
    val mutational_context

    output:
    path "Signature_Spectra"

    script:
    """
    python3 SigProfilerPlotting.py --mutational_context "${params.mutational_context}" --deNovoSignatures_matrix "${mSigHdp.out}"
    """
}