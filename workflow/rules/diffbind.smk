import pandas as pd
import os

OUTDIR = config["output_dir"]
samples_df = pd.read_csv(config["samples"], sep="\t").set_index("name", drop=False)
SAMPLES = samples_df.index.tolist()

# 1. Count ALL samples together to form the consensus peakset
rule diffbind_count:
    input:
        samples = config["samples"],
        peaks = expand(os.path.join(OUTDIR, "6.peak_calls", "{sample}_peaks.narrowPeak"), sample=SAMPLES),
        bams = expand(os.path.join(OUTDIR, "4.filtered", "{sample}.final.bam"), sample=SAMPLES)
    output:
        rds = os.path.join(OUTDIR, "8.diffbind", "{contrast}_counted.rds")
    params:
        outdir = OUTDIR, # We pass this to R so it can find the BAMs correctly!
        column = lambda wildcards: config["diffbind_contrasts"][wildcards.contrast]["column"],
        target = lambda wildcards: config["diffbind_contrasts"][wildcards.contrast]["target"],
        reference = lambda wildcards: config["diffbind_contrasts"][wildcards.contrast]["reference"]
    threads: 4
    resources:
        mem_mb = 64000,
        runtime = 360
    log:
        os.path.join(OUTDIR, "logs", "diffbind", "{contrast}_count.log")
    script:
        "/hpc/local/Rocky8/prekovic/theo/ATAC-seq_pipeline_WIP/scripts/run_diffbind_count.R"

# 2. Run pairwise comparisons dynamically based on config.yaml
rule diffbind_analyze:
    input:
        rds = os.path.join(OUTDIR, "8.diffbind", "{contrast}_counted.rds")
    output:
        diff = os.path.join(OUTDIR, "8.diffbind", "{contrast}_DiffPeaks.csv"),
        bed_up = os.path.join(OUTDIR, "8.diffbind", "{contrast}_UP.bed"),
        bed_down = os.path.join(OUTDIR, "8.diffbind", "{contrast}_DOWN.bed"),
        volc = os.path.join(OUTDIR, "8.diffbind", "{contrast}_Volcano.pdf"),
        pca = os.path.join(OUTDIR, "8.diffbind", "{contrast}_PCA.pdf")
    params:
        column = lambda wildcards: config["diffbind_contrasts"][wildcards.contrast]["column"],
        target = lambda wildcards: config["diffbind_contrasts"][wildcards.contrast]["target"],
        reference = lambda wildcards: config["diffbind_contrasts"][wildcards.contrast]["reference"]
    container:
        "docker://theo46/genomics-r-env:v1"
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        runtime = 60
    log:
        os.path.join(OUTDIR, "logs", "diffbind", "{contrast}_analyze.log")
    script:
        "/hpc/local/Rocky8/prekovic/theo/ATAC-seq_pipeline_WIP/scripts/run_diffbind_analyze.R"