process MULTIQC_BIOTYPE_COUNTS_QUANTIFICATION {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path counts          // merged gene count matrix (gene_id × samples TSV)
    path gtf             // genome annotation GTF (used to look up gene_biotype / gene_type)
    val  quant_label     // label for plot title/id, e.g. 'star_salmon', 'rsem', 'salmon', 'kallisto'
    val  biotype_attr    // GTF attribute to use for biotype, e.g. 'gene_biotype' or 'gene_type'
    val  min_count       // minimum per-sample count to include a gene (default: 5)

    output:
    path "*.tsv"       , emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'biotype_counts_quantification.py'

    stub:
    """
    touch biotype_counts_quantification_${quant_label}_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
