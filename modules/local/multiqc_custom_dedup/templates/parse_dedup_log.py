#!/usr/bin/env python3
"""
Parse UMI deduplication log (umitools or UMICollapse) for MultiQC custom bargraph.

Produces a per-sample TSV appended to a MultiQC custom_content header that renders
as a stacked bargraph showing unique vs duplicate reads in the transcriptome BAM.
"""

import re
import sys

log_file = "$dedup_log"
sample_id = "$meta.id"
prefix    = "$prefix"
header    = "$header"


def parse_umitools_log(content):
    """
    UMI-tools dedup log lines of interest:
        INFO ... dedup - Reads: Input Reads: 1000000
        INFO ... dedup - Number of reads out: 800000
    """
    reads_in  = None
    reads_out = None
    for line in content.splitlines():
        m = re.search(r"Reads: Input Reads:\s+(\d+)", line)
        if m:
            reads_in = int(m.group(1))
        m = re.search(r"Number of reads out:\s+(\d+)", line)
        if m:
            reads_out = int(m.group(1))
    return reads_in, reads_out


def parse_umicollapse_log(content):
    """
    UMICollapse writes its summary statistics to stdout (tee'd to the log).
    Multiple output formats are tried to handle version differences:
        Total reads: N, reads output: M
        Total input: N / Total output: M
        Reads written (unique): M of N
    """
    reads_in  = None
    reads_out = None
    for line in content.splitlines():
        # Format: "Total reads: 1000000, reads output: 800000"
        m = re.search(r"Total reads:\s*([\d,]+),\s*reads output:\s*([\d,]+)", line, re.IGNORECASE)
        if m:
            reads_in  = int(m.group(1).replace(",", ""))
            reads_out = int(m.group(2).replace(",", ""))
            break

        # Format: "Total input: 1000000" / "Total output: 800000"
        m = re.search(r"Total input[^:]*:\s*([\d,]+)", line, re.IGNORECASE)
        if m and reads_in is None:
            reads_in = int(m.group(1).replace(",", ""))
        m = re.search(r"Total output[^:]*:\s*([\d,]+)", line, re.IGNORECASE)
        if m and reads_out is None:
            reads_out = int(m.group(1).replace(",", ""))

        # Format: "Reads written (unique): 800000 of 1000000"
        m = re.search(r"Reads written[^:]*:\s*([\d,]+)\s+of\s+([\d,]+)", line, re.IGNORECASE)
        if m:
            reads_out = int(m.group(1).replace(",", ""))
            reads_in  = int(m.group(2).replace(",", ""))
            break

    return reads_in, reads_out


# --- parse log ---------------------------------------------------------------

with open(log_file) as fh:
    content = fh.read()

if "_UMICollapse.log" in log_file:
    reads_in, reads_out = parse_umicollapse_log(content)
else:
    reads_in, reads_out = parse_umitools_log(content)

if reads_in is None or reads_out is None:
    sys.exit(
        f"ERROR: Could not parse reads_in ({reads_in}) or reads_out ({reads_out}) "
        f"from {log_file}. Please report this at https://github.com/nf-core/rnaseq/issues"
    )

dup_reads = reads_in - reads_out

# --- write output ------------------------------------------------------------

data_file = prefix + ".umi_dedup_transcriptome_mqc.tsv"

with open(header) as fh:
    header_text = fh.read()

with open(data_file, "w") as fh:
    fh.write(header_text)
    fh.write("Sample\tUnique reads\tDuplicate reads\n")
    fh.write(f"{sample_id}\t{reads_out}\t{dup_reads}\n")

# --- versions ----------------------------------------------------------------

with open("versions.yml", "w") as fh:
    fh.write('"${task.process}":\n')
    import subprocess
    ver = subprocess.check_output(["python", "--version"]).decode().strip().replace("Python ", "")
    fh.write(f"    python: {ver}\n")
