process MULTIQC_CUSTOM_DEDUP {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path pre_flagstats, stageAs: 'pre/*'
    path post_flagstats, stageAs: 'post/*'
    path header

    output:
    path "*.tsv"       , emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 - << 'PYEOF'
import os, re, sys


def parse_primary_mapped(flagstat_file):
    with open(flagstat_file) as fh:
        for line in fh:
            if 'primary mapped' in line:
                m = re.match(r'(\\d+)', line)
                if m:
                    return int(m.group(1))
    return None


def sample_id_from_flagstat(fname):
    name = os.path.basename(fname)
    name = re.sub(r'\\.flagstat', '', name)
    name = re.sub(r'\\.sorted\\.bam', '', name)
    name = re.sub(r'\\.umi_dedup', '', name)
    return name


pre_dir  = 'pre'
post_dir = 'post'

pre_files  = {sample_id_from_flagstat(f): os.path.join(pre_dir, f)
              for f in os.listdir(pre_dir) if f.endswith('.flagstat')}
post_files = {sample_id_from_flagstat(f): os.path.join(post_dir, f)
              for f in os.listdir(post_dir) if f.endswith('.flagstat')}

samples = sorted(set(pre_files) & set(post_files))
if not samples:
    sys.exit(
        'ERROR: No matching samples found between pre-dedup ('
        + str(sorted(pre_files)) + ') and post-dedup ('
        + str(sorted(post_files)) + ') flagstat files.'
    )

with open('${header}') as fh:
    header_text = fh.read()

rows = []
for sample in samples:
    pre_count  = parse_primary_mapped(pre_files[sample])
    post_count = parse_primary_mapped(post_files[sample])
    if pre_count is None or post_count is None:
        sys.exit(
            'ERROR: Could not parse primary mapped reads for sample '
            + sample + ' (pre=' + str(pre_count) + ', post=' + str(post_count) + ')'
        )
    dup_reads = pre_count - post_count
    rows.append(sample + '\\t' + str(post_count) + '\\t' + str(dup_reads))

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
