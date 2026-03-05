process MULTIQC_CUSTOM_DEDUP {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path dedup_logs
    path header

    output:
    path "*.tsv"       , emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 - << 'PYEOF'
import glob, os, re, sys

header    = '${header}'
log_files = glob.glob('*_dedup*.log') + glob.glob('*.dedup*.log')
if not log_files:
    log_files = [f for f in os.listdir('.') if f.endswith('.log')]


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


def sample_id_from_filename(fname):
    name = os.path.basename(fname)
    name = re.sub(r'\\.umi_dedup\\..*', '', name)
    name = re.sub(r'\\.log\$', '', name)
    return name


with open(header) as fh:
    header_text = fh.read()

rows = []
for log_file in sorted(log_files):
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

    sample_id = sample_id_from_filename(log_file)
    dup_reads = reads_in - reads_out
    rows.append(sample_id + '\\t' + str(reads_out) + '\\t' + str(dup_reads))

with open('umi_dedup_genome_mqc.tsv', 'w') as fh:
    fh.write(header_text)
    fh.write('Sample\\tUnique reads\\tDuplicate reads\\n')
    for row in rows:
        fh.write(row + '\\n')
PYEOF

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """

    stub:
    """
    touch umi_dedup_genome_mqc.tsv

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """
}
