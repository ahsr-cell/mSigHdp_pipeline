process mSigHdp_hierarchy {
    
    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    path hierarchy_matrix
    val hierarchy_parameter
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
    path "deNovo_signatures/mSigHdp_deNovoSigs_sigPADecomp.txt", emit: deNovo_extsigs_sigPA
    path "deNovo_signatures/mSigHdp_lowConfSignatures.txt", emit: deNovo_lowconfsigs, optional: true

    script:
    """
    Rscript --vanilla ${projectDir}/bin/mSigHdp_hierarchy.R -hierarchy ${hierarchy_matrix} -hp ${hierarchy_parameter} -c ${mutational_context} -a ${analysis_type} -b ${burnin_iterations} -x ${burnin_multiplier} -o ${posterior} -i ${posterior_iterations} ${mutational_matrix}
    """
}