process MULTIQC_CUSTOM_POSTDEDUP_CATEGORIES {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
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


def parse_flagstat(flagstat_file):
    counts = {}
    with open(flagstat_file) as fh:
        for line in fh:
            m = re.match(r'(\\d+)', line)
            if not m:
                continue
            val = int(m.group(1))
            if 'primary mapped' in line:
                counts['primary_mapped'] = val
            elif 'secondary' in line:
                counts['secondary'] = val
            elif 'supplementary' in line:
                counts['supplementary'] = val
    return counts


def sample_id_from_flagstat(fname):
    name = os.path.basename(fname)
    name = re.sub(r'\\.flagstat', '', name)
    name = re.sub(r'\\.sorted\\.bam', '', name)
    name = re.sub(r'\\.umi_dedup', '', name)
    return name


post_dir = 'post'

post_files = {sample_id_from_flagstat(f): os.path.join(post_dir, f)
              for f in os.listdir(post_dir) if f.endswith('.flagstat')}

samples = sorted(post_files)
if not samples:
    sys.exit('ERROR: No flagstat files found in post-dedup directory.')

with open('${header}') as fh:
    header_text = fh.read()

columns = [
    'Primary mapped',
    'Secondary',
    'Supplementary',
]

rows = []
for sample in samples:
    counts = parse_flagstat(post_files[sample])
    values = [
        counts.get('primary_mapped', 0),
        counts.get('secondary', 0),
        counts.get('supplementary', 0),
    ]
    rows.append(sample + '\\t' + '\\t'.join(str(v) for v in values))

with open('star_dedup_categories_mqc.tsv', 'w') as fh:
    fh.write(header_text)
    fh.write('Sample\\t' + '\\t'.join(columns) + '\\n')
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
    touch star_dedup_categories_mqc.tsv

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """
}
