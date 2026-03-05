process MULTIQC_CUSTOM_DEDUP {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(dedup_log)
    path  header

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - << 'PYEOF'
import re, sys

log_file  = '${dedup_log}'
sample_id = '${meta.id}'
prefix    = '${prefix}'
header    = '${header}'


def parse_umitools_log(content):
    reads_in  = None
    reads_out = None
    for line in content.splitlines():
        m = re.search(r'Reads: Input Reads:\\s+(\\d+)', line)
        if m:
            reads_in = int(m.group(1))
        m = re.search(r'Number of reads out:\\s+(\\d+)', line)
        if m:
            reads_out = int(m.group(1))
    return reads_in, reads_out


def parse_umicollapse_log(content):
    reads_in  = None
    reads_out = None
    for line in content.splitlines():
        m = re.search(r'Total reads:\\s*([\\d,]+),\\s*reads output:\\s*([\\d,]+)', line, re.IGNORECASE)
        if m:
            reads_in  = int(m.group(1).replace(',', ''))
            reads_out = int(m.group(2).replace(',', ''))
            break
        m = re.search(r'Number of input reads\\s+([\\d,]+)', line, re.IGNORECASE)
        if m and reads_in is None:
            reads_in = int(m.group(1).replace(',', ''))
        m = re.search(r'Number of reads after deduplicat\\w*\\s+([\\d,]+)', line, re.IGNORECASE)
        if m and reads_out is None:
            reads_out = int(m.group(1).replace(',', ''))
        m = re.search(r'Total input[^:]*:\\s*([\\d,]+)', line, re.IGNORECASE)
        if m and reads_in is None:
            reads_in = int(m.group(1).replace(',', ''))
        m = re.search(r'Total output[^:]*:\\s*([\\d,]+)', line, re.IGNORECASE)
        if m and reads_out is None:
            reads_out = int(m.group(1).replace(',', ''))
        m = re.search(r'Reads written[^:]*:\\s*([\\d,]+)\\s+of\\s+([\\d,]+)', line, re.IGNORECASE)
        if m:
            reads_out = int(m.group(1).replace(',', ''))
            reads_in  = int(m.group(2).replace(',', ''))
            break
    return reads_in, reads_out


with open(log_file) as fh:
    content = fh.read()

if '_UMICollapse.log' in log_file:
    reads_in, reads_out = parse_umicollapse_log(content)
else:
    reads_in, reads_out = parse_umitools_log(content)

if reads_in is None or reads_out is None:
    sys.exit(
        'ERROR: Could not parse reads_in (' + str(reads_in) + ') or reads_out ('
        + str(reads_out) + ') from ' + log_file
        + '. Please report this at https://github.com/nf-core/rnaseq/issues'
    )

dup_reads = reads_in - reads_out

with open(header) as fh:
    header_text = fh.read()

out_file = prefix + '.umi_dedup_genome_mqc.tsv'
with open(out_file, 'w') as fh:
    fh.write(header_text)
    fh.write('Sample\\tUnique reads\\tDuplicate reads\\n')
    fh.write(sample_id + '\\t' + str(reads_out) + '\\t' + str(dup_reads) + '\\n')
PYEOF

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_dedup_genome_mqc.tsv

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """
}
