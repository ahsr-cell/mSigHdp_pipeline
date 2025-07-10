process mSigHdp {

    publishDir "${params.outdir}/deNovo_signatures", mode: "copy", overwrite: true

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
    path "deNovo_signatures"

    script:
    """
    Rscript --vanilla mSigHdp.R -s "${params.sample_matrix}" -c "${params.mutational_context}" -a "${params.analysis_type}" -b "${params.burnin_iterations}" -x "${params.burnin_multiplier}" -o "${params.posterior}" -i "${params.posterior_iterations}" "${params.mutational_matrix}"
    """
}