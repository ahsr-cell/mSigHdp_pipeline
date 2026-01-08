process MutMatrix_resourcereqs_hierarchy {

    input:
    path mutational_matrix
    path hierarchy_matrix
    val hierarchy_parameter
    val user_defmemory

    output:
    path "memory_requirements/memory_requirements.csv", emit: memory_reqs_matrix

    script:
    """
    Rscript --vanilla ${projectDir}/bin/MutMatrix_resourcereqs_hierarchy.R -hierarchy ${hierarchy_matrix} -hp ${hierarchy_parameter} -umem ${user_defmemory} ${mutational_matrix}
    """
}