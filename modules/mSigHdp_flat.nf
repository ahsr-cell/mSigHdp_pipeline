process mSigHdp_flat {
    
    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    tuple val(Sample_number), val(Mutation_burden), val(Memory_required)
    val mutational_context
    val analysis_type
    val burnin_iterations 
    val burnin_multiplier
    val posterior 
    val posterior_iterations
    val chains
    val clusters
    val alpha
    val beta
    val confidence
    path mutational_matrix

    output:
    path "deNovo_signatures", emit: deNovo_signaturesdir
    path "deNovo_signatures/mSigHdp_deNovoSignatures.txt", emit: deNovo_extractedsigs
    path "deNovo_signatures/mSigHdp_deNovoSigs_sigPADecomp.txt", emit: deNovo_extsigs_sigPA
    path "deNovo_signatures/mSigHdp_lowConfSignatures.txt", emit: deNovo_lowconfsigs, optional: true

    script:
    """
    Rscript --vanilla ${projectDir}/bin/mSigHdp_flat.R -c ${mutational_context} -a ${analysis_type} -b ${burnin_iterations} -x ${burnin_multiplier} -o ${posterior} -i ${posterior_iterations} -ch ${chains} -k ${clusters} -ga ${alpha} -gb ${beta} -h ${confidence} ${mutational_matrix}
    """
}