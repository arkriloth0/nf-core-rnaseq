process MULTIQC_CUSTOM_POSTDEDUP_CATEGORIES {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    path bams, stageAs: 'bams/*'
    path bais, stageAs: 'bams/*'
    path header

    output:
    path "*.tsv"       , emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Build TSV from header
    cp ${header} postdedup_categories_mqc.tsv

    # Column header
    printf 'Sample\\tUniquely mapped\\tMulti-mapped\\tUnmapped\\n' >> postdedup_categories_mqc.tsv

    for bam in bams/*.bam; do
        sample=\$(basename "\$bam" .bam | sed 's/\\.sorted\\.bam\$//; s/\\.umi_dedup//')

        # Primary mapped with MAPQ=255 (STAR uniquely mapped)
        unique=\$(samtools view -c -F 0x904 -q 255 "\$bam")

        # Total primary mapped
        total_primary=\$(samtools view -c -F 0x904 "\$bam")

        # Multi-mapped = primary mapped with MAPQ < 255
        multi=\$(( total_primary - unique ))

        # Unmapped reads
        unmapped=\$(samtools view -c -f 4 "\$bam")

        printf '%s\\t%s\\t%s\\t%s\\n' "\$sample" "\$unique" "\$multi" "\$unmapped" >> postdedup_categories_mqc.tsv
    done

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """

    stub:
    """
    touch postdedup_categories_mqc.tsv

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
END_VERSIONS
    """
}
