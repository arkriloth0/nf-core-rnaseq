process MULTIQC_CUSTOM_STAR_DEDUP {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path star_logs, stageAs: 'star/*'
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


def parse_star_log(log_file):
    categories = {}
    with open(log_file) as fh:
        for line in fh:
            line = line.strip()
            if 'Uniquely mapped reads number' in line:
                categories['uniquely_mapped'] = int(line.split('|')[-1].strip())
            elif 'Number of reads mapped to multiple loci' in line:
                categories['multi_mapped'] = int(line.split('|')[-1].strip())
            elif 'Number of reads mapped to too many loci' in line:
                categories['too_many_loci'] = int(line.split('|')[-1].strip())
            elif 'Number of reads unmapped: too many mismatches' in line:
                categories['unmapped_mismatches'] = int(line.split('|')[-1].strip())
            elif 'Number of reads unmapped: too short' in line:
                categories['unmapped_tooshort'] = int(line.split('|')[-1].strip())
            elif 'Number of reads unmapped: other' in line:
                categories['unmapped_other'] = int(line.split('|')[-1].strip())
    return categories


def parse_primary_mapped(flagstat_file):
    with open(flagstat_file) as fh:
        for line in fh:
            if 'primary mapped' in line:
                m = re.match(r'(\\d+)', line)
                if m:
                    return int(m.group(1))
    return None


def sample_id_from_star_log(fname):
    name = os.path.basename(fname)
    name = re.sub(r'\\.Log\\.final\\.out', '', name)
    return name


def sample_id_from_flagstat(fname):
    name = os.path.basename(fname)
    name = re.sub(r'\\.flagstat', '', name)
    name = re.sub(r'\\.sorted\\.bam', '', name)
    name = re.sub(r'\\.umi_dedup', '', name)
    return name


star_dir = 'star'
pre_dir  = 'pre'
post_dir = 'post'

star_files = {sample_id_from_star_log(f): os.path.join(star_dir, f)
              for f in os.listdir(star_dir) if f.endswith('.Log.final.out')}
pre_files  = {sample_id_from_flagstat(f): os.path.join(pre_dir, f)
              for f in os.listdir(pre_dir) if f.endswith('.flagstat')}
post_files = {sample_id_from_flagstat(f): os.path.join(post_dir, f)
              for f in os.listdir(post_dir) if f.endswith('.flagstat')}

samples = sorted(set(star_files) & set(pre_files) & set(post_files))
if not samples:
    sys.exit(
        'ERROR: No matching samples found between STAR logs ('
        + str(sorted(star_files)) + '), pre-dedup ('
        + str(sorted(pre_files)) + ') and post-dedup ('
        + str(sorted(post_files)) + ') flagstat files.'
    )

with open('${header}') as fh:
    header_text = fh.read()

columns = [
    'Uniquely mapped (unique)',
    'Uniquely mapped (duplicates)',
    'Mapped to multiple loci',
    'Mapped to too many loci',
    'Unmapped (too many mismatches)',
    'Unmapped (too short)',
    'Unmapped (other)',
]

rows = []
for sample in samples:
    cats = parse_star_log(star_files[sample])
    pre_count  = parse_primary_mapped(pre_files[sample])
    post_count = parse_primary_mapped(post_files[sample])
    if pre_count is None or post_count is None:
        sys.exit(
            'ERROR: Could not parse primary mapped reads for sample '
            + sample + ' (pre=' + str(pre_count) + ', post=' + str(post_count) + ')'
        )

    dup_reads = pre_count - post_count
    unique_mapped = cats.get('uniquely_mapped', 0) - dup_reads

    values = [
        unique_mapped,
        dup_reads,
        cats.get('multi_mapped', 0),
        cats.get('too_many_loci', 0),
        cats.get('unmapped_mismatches', 0),
        cats.get('unmapped_tooshort', 0),
        cats.get('unmapped_other', 0),
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
