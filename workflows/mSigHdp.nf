process mSigHdp {

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    path sample_matrix
    val mutational_context
    val analysis_type
    val burnin_iterations 
    val burnin_multiplier
    val posterior 
    val posterior_iterations
    path mutational_matrix

    output:
    path "deNovo_signatures", emit: deNovo_signaturesdir
    path "deNovo_signatures/mSigHdp_deNovoSignatures.txt", emit: deNovo_extractedsigs
    path "deNovo_signatures/low.confidence.signatures.csv", emit: deNovo_lowconfsigs, optional: true

    script:
    """
    Rscript --vanilla ${projectDir}/bin/mSigHdp.R -s ${sample_matrix} -c ${mutational_context} -a ${analysis_type} -b ${burnin_iterations} -x ${burnin_multiplier} -o ${posterior} -i ${posterior_iterations} ${mutational_matrix}
    """
}