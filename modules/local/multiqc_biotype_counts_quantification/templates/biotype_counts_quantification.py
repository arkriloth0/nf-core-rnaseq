#!/usr/bin/env python3
"""
Build a MultiQC custom bargraph showing, per sample, the number of genes with
>= min_count raw counts broken down by gene biotype.

Inputs
------
counts      : merged gene count matrix (TSV).  Expected formats:
              - RSEM  : gene_id <tab> transcript_id(s) <tab> sample1 ... sampleN
              - Salmon/Kallisto (tximeta) : gene_id <tab> sample1 ... sampleN
gtf         : genome annotation GTF (gene lines used to extract biotype)
quant_label : short label inserted into section id / title (e.g. 'star_salmon')
biotype_attr: GTF attribute name for biotype ('gene_biotype' or 'gene_type')
min_count   : integer minimum per-sample count threshold (default 5)

Output
------
biotype_counts_quantification_<quant_label>_mqc.tsv
    MultiQC custom_content bargraph (counts/% toggle via cpswitch).
"""

import re
import sys
from collections import defaultdict

counts_file  = "$counts"
gtf_file     = "$gtf"
quant_label  = "$quant_label"
biotype_attr = "$biotype_attr"
min_count    = int("$min_count")

# ---------------------------------------------------------------------------
# 1. Parse GTF to build gene_id -> biotype mapping
# ---------------------------------------------------------------------------

gene_biotype = {}
attr_re = re.compile(r'%s\s+"([^"]+)"' % re.escape(biotype_attr))

with open(gtf_file) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        attrs = fields[8]
        gid_m = re.search(r'gene_id\s+"([^"]+)"', attrs)
        bio_m = attr_re.search(attrs)
        if gid_m and bio_m:
            gene_biotype[gid_m.group(1)] = bio_m.group(1)

if not gene_biotype:
    sys.exit(
        f"ERROR: No gene biotype annotations found in {gtf_file} "
        f"using attribute '{biotype_attr}'. "
        f"Check --featurecounts_group_type / --gencode settings."
    )

# ---------------------------------------------------------------------------
# 2. Parse count matrix
# ---------------------------------------------------------------------------

with open(counts_file) as fh:
    header_line = fh.readline().rstrip("\n").split("\t")

# Detect RSEM format: second column is 'transcript_id(s)'
if header_line[1].lower().startswith("transcript_id"):
    gene_col      = 0
    sample_start  = 2
else:
    gene_col      = 0
    sample_start  = 1

sample_names = header_line[sample_start:]

# per_sample_biotype[sample][biotype] = count of genes >= min_count
per_sample_biotype = {s: defaultdict(int) for s in sample_names}

with open(counts_file) as fh:
    fh.readline()  # skip header
    for line in fh:
        parts   = line.rstrip("\n").split("\t")
        gene_id = parts[gene_col]
        biotype = gene_biotype.get(gene_id, "other")
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

safe_label  = quant_label.replace("-", "_")
out_file    = f"biotype_counts_quantification_{safe_label}_mqc.tsv"

# Human-readable tool name for headers
tool_names  = {
    "star_salmon" : "STAR + Salmon",
    "rsem"        : "RSEM",
    "salmon"      : "Salmon",
    "kallisto"    : "Kallisto",
}
tool_name = tool_names.get(quant_label, quant_label)

with open(out_file, "w") as fh:
    fh.write(f"# id: 'biotype_counts_quantification_{safe_label}'\n")
    fh.write(f"# section_name: 'Detected Genes by Biotype ({tool_name})'\n")
    fh.write(
        f"# description: >\n"
        f"#     Number of genes with &ge;{min_count} counts in each sample, broken down by gene biotype\n"
        f"#     as annotated in the supplied GTF (attribute <code>{biotype_attr}</code>).\n"
        f"#     Counts are drawn from the <b>{tool_name}</b> quantification output.\n"
        f"#     Only genes detected in at least one sample are shown.\n"
        f"#     <br><br>\n"
        f"#     <b>How this differs from featureCounts Biotypes:</b> The featureCounts biotype\n"
        f"#     section counts <em>reads</em> overlapping genomic features on the <em>genome BAM</em>,\n"
        f"#     making it a measure of library composition (e.g. rRNA contamination).\n"
        f"#     This section counts <em>distinct genes</em> passing a minimum expression\n"
        f"#     threshold in the final quantification output, making it a measure of\n"
        f"#     transcriptome complexity and gene detection sensitivity per sample.\n"
        f"#     For <code>--with_umi</code> samples the quantification is post-deduplication.\n"
        f"# plot_type: 'bargraph'\n"
        f"# anchor: 'biotype_counts_quantification_{safe_label}'\n"
        f"# pconfig:\n"
        f"#     id: 'biotype_counts_quantification_{safe_label}_plot'\n"
        f"#     title: 'Detected Genes by Biotype: {tool_name}'\n"
        f"#     cpswitch_counts_label: 'Number of Genes'\n"
        f"#     hide_zero_cats: true\n"
    )

    # Header row: Sample <tab> biotype1 <tab> biotype2 ...
    fh.write("Sample\t" + "\t".join(biotypes_ordered) + "\n")

    for sample in sample_names:
        row = [sample] + [str(per_sample_biotype[sample].get(bt, 0)) for bt in biotypes_ordered]
        fh.write("\t".join(row) + "\n")

# ---------------------------------------------------------------------------
# 5. Write versions.yml
# ---------------------------------------------------------------------------

import subprocess
ver = subprocess.check_output(["python", "--version"]).decode().strip().replace("Python ", "")
with open("versions.yml", "w") as fh:
    fh.write('"${task.process}":\n')
    fh.write(f"    python: {ver}\n")
