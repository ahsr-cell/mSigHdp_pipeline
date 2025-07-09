process SigProfilerAssignment {

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input: 
    path deNovo_signatures
    path mutational_matrix

    output:
    path "SigProfilerDecomposition"

    script:
    """
    python3 SigProfilerAssignment.py --mutational_matrix "${params.mutational_matrix}" --deNovoSignatures_matrix "${mSigHdp.out}"
    """
}