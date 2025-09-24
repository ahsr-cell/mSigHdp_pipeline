process MutMatrix_resourcereqs_hierarchy {

    input:
    path hierarchy_matrix
    path mutational_matrix

    output:
    path "memory_requirements/memory_requirements.csv", emit: memory_reqs_matrix

    script:
    """
    Rscript --vanilla ${projectDir}/bin/MutMatrix_resourcereqs_hierarchy.R -hierarchy ${hierarchy_matrix} ${mutational_matrix}
    """
}