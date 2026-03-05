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
    """
    python3 - << 'PYEOF'
import re, sys
from collections import defaultdict

counts_file  = '${counts}'
gtf_file     = '${gtf}'
quant_label  = '${quant_label}'
biotype_attr = '${biotype_attr}'
min_count    = int('${min_count}')

# ---------------------------------------------------------------------------
# 1. Parse GTF to build gene_id -> biotype mapping
# ---------------------------------------------------------------------------

gene_biotype = {}
attr_re = re.compile(r'%s\\s+"([^"]+)"' % re.escape(biotype_attr))

for feature_type in ('gene', 'transcript'):
    with open(gtf_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\\n').split('\\t')
            if len(fields) < 9 or fields[2] != feature_type:
                continue
            attrs = fields[8]
            gid_m = re.search(r'gene_id\\s+"([^"]+)"', attrs)
            bio_m = attr_re.search(attrs)
            if gid_m and bio_m:
                gene_biotype.setdefault(gid_m.group(1), bio_m.group(1))
    if gene_biotype:
        break

if not gene_biotype:
    sys.exit(
        'ERROR: No gene biotype annotations found in ' + gtf_file
        + ' using attribute "' + biotype_attr + '"'
        + '. Check --featurecounts_group_type / --gencode settings.'
    )

# ---------------------------------------------------------------------------
# 2. Parse count matrix
# ---------------------------------------------------------------------------

with open(counts_file) as fh:
    header_line = fh.readline().rstrip('\\n').split('\\t')

# Detect RSEM format: second column is 'transcript_id(s)'
if len(header_line) > 1 and header_line[1].lower().startswith('transcript_id'):
    gene_col     = 0
    sample_start = 2
else:
    gene_col     = 0
    sample_start = 1

sample_names = header_line[sample_start:]
per_sample_biotype = {s: defaultdict(int) for s in sample_names}

with open(counts_file) as fh:
    fh.readline()  # skip header
    for line in fh:
        parts   = line.rstrip('\\n').split('\\t')
        gene_id = parts[gene_col]
        biotype = gene_biotype.get(gene_id, 'other')
        counts  = parts[sample_start:]
        for sample, cnt in zip(sample_names, counts):
            try:
                if float(cnt) >= min_count:
                    per_sample_biotype[sample][biotype] += 1
            except ValueError:
                pass

# ---------------------------------------------------------------------------
# 3. Collect ordered biotypes (by total gene count across all samples)
# ---------------------------------------------------------------------------

biotype_totals = defaultdict(int)
for sample in sample_names:
    for bt, n in per_sample_biotype[sample].items():
        biotype_totals[bt] += n

biotypes_ordered = sorted(biotype_totals, key=lambda b: -biotype_totals[b])

# ---------------------------------------------------------------------------
# 4. Write MultiQC custom_content file
# ---------------------------------------------------------------------------

safe_label = quant_label.replace('-', '_')
out_file   = 'biotype_counts_quantification_' + safe_label + '_mqc.tsv'

tool_names = {
    'star_salmon': 'STAR + Salmon',
    'rsem':        'RSEM',
    'salmon':      'Salmon',
    'kallisto':    'Kallisto',
}
tool_name = tool_names.get(quant_label, quant_label)

lines = []
lines.append("# id: 'biotype_counts_quantification_" + safe_label + "'")
lines.append("# section_name: 'Detected Genes by Biotype (" + tool_name + ")'")
lines.append("# description: >")
lines.append("#     Number of genes with &ge;" + str(min_count) + " counts in each sample, broken down by gene biotype")
lines.append("#     as annotated in the supplied GTF (attribute <code>" + biotype_attr + "</code>).")
lines.append("#     Counts are drawn from the <b>" + tool_name + "</b> quantification output.")
lines.append("#     Only genes detected in at least one sample are shown.")
lines.append("#     <br><br>")
lines.append("#     <b>How this differs from featureCounts Biotypes:</b> The featureCounts biotype")
lines.append("#     section counts <em>reads</em> overlapping genomic features on the <em>genome BAM</em>,")
lines.append("#     making it a measure of library composition (e.g. rRNA contamination).")
lines.append("#     This section counts <em>distinct genes</em> passing a minimum expression")
lines.append("#     threshold in the final quantification output, making it a measure of")
lines.append("#     transcriptome complexity and gene detection sensitivity per sample.")
lines.append("#     For <code>--with_umi</code> samples the quantification is post-deduplication.")
lines.append("# plot_type: 'bargraph'")
lines.append("# anchor: 'biotype_counts_quantification_" + safe_label + "'")
lines.append("# pconfig:")
lines.append("#     id: 'biotype_counts_quantification_" + safe_label + "_plot'")
lines.append("#     title: 'Detected Genes by Biotype: " + tool_name + "'")
lines.append("#     cpswitch_counts_label: 'Number of Genes'")
lines.append("#     hide_zero_cats: true")
lines.append("Sample\\t" + "\\t".join(biotypes_ordered))
for sample in sample_names:
    row = [sample] + [str(per_sample_biotype[sample].get(bt, 0)) for bt in biotypes_ordered]
    lines.append("\\t".join(row))

with open(out_file, 'w') as fh:
    fh.write("\\n".join(lines) + "\\n")
PYEOF

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """

    stub:
    """
    touch biotype_counts_quantification_${quant_label}_mqc.tsv

    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """
}
