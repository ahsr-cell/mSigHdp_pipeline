process SigProfilerAssignment {

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input: 
    path deNovoSignatures_matrix
    path mutational_matrix

    output:
    path "SigProfilerDecomposition"

    script:
    """
    export MPLCONFIGDIR=${workDir}
    rm -rf SigProfilerDecomposition
    mkdir SigProfilerDecomposition
    
    SigProfilerAssignment_wrapper.py --mutational_matrix ${mutational_matrix} --deNovoSignatures_matrix ${deNovoSignatures_matrix} --output_directory SigProfilerDecomposition/
    """
}