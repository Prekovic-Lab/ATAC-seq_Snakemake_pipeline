import pandas as pd
import os
import re

wildcard_constraints:
    condition="[^/]+",
    sample="[^/]+"

# Load the samples file
# We use 'dtype=str' to prevent replicates like "1" being read as integers
samples_df = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("name", drop=False, verify_integrity=True)

# --- Define Constraints (Global) ---

# Create lists
PE_SAMPLES_LIST = pd.Index(samples_df.loc[
    (samples_df["seq_type"] == "paired") & 
    (samples_df["file_2"].notna()) & 
    (samples_df["file_2"] != "-")
]["name"]).unique().tolist()

SE_SAMPLES_LIST = [s for s in samples_df["name"] if s not in PE_SAMPLES_LIST]

# Create Regex Strings for Constraints
# If the list is empty, use a placeholder that will never match a real sample name.
#PE_SAMPLES = "|".join(PE_SAMPLES_LIST) if PE_SAMPLES_LIST else "NO_PE_SAMPLES_FOUND"                           <-- trying to see if + works now
PE_SAMPLES = "|".join(re.escape(s) for s in PE_SAMPLES_LIST) if PE_SAMPLES_LIST else "NO_PE_SAMPLES_FOUND"
#SE_SAMPLES = "|".join(SE_SAMPLES_LIST) if SE_SAMPLES_LIST else "NO_SE_SAMPLES_FOUND"                           <-- trying to see if + works now
SE_SAMPLES = "|".join(re.escape(s) for s in SE_SAMPLES_LIST) if SE_SAMPLES_LIST else "NO_SE_SAMPLES_FOUND"


# --- Peak Type Constraints ---
# Get lists based on the 'peak_type' column
NARROW_SAMPLES_LIST = samples_df[samples_df["peak_type"] == "narrow"]["name"].unique().tolist()
BROAD_SAMPLES_LIST = samples_df[samples_df["peak_type"] == "broad"]["name"].unique().tolist()

# Create Regex strings for wildcard_constraints
# If list is empty, use a placeholder that matches nothing
#NARROW_SAMPLES = "|".join(NARROW_SAMPLES_LIST) if NARROW_SAMPLES_LIST else "NO_NARROW_SAMPLES"                           <-- trying to see if + works now
NARROW_SAMPLES = "|".join(re.escape(s) for s in NARROW_SAMPLES_LIST) if NARROW_SAMPLES_LIST else "NO_NARROW_SAMPLES"
#BROAD_SAMPLES = "|".join(BROAD_SAMPLES_LIST) if BROAD_SAMPLES_LIST else "NO_BROAD_SAMPLES"                               <-- trying to see if + works now
BROAD_SAMPLES = "|".join(re.escape(s) for s in BROAD_SAMPLES_LIST) if BROAD_SAMPLES_LIST else "NO_BROAD_SAMPLES"

# --- Helper Functions ---

# Helper function to get FASTQ paths based on the sample wildcard
def get_fastqs(wildcards):
    row = samples_df.loc[wildcards.sample]
    r1 = os.path.join(config["fastq_dir"], row["file_1"])
    if row["seq_type"] == "single" or pd.isna(row["file_2"]) or row["file_2"] == "-":
        return [r1]
    else:
        r2 = os.path.join(config["fastq_dir"], row["file_2"])
        return [r1, r2]

# define the output targets for trimming step
def get_trimmed_targets(wildcards):
    targets = []
    for sample, row in samples_df.iterrows():
        # Trim Galore Naming Convention
        if row["seq_type"] == "single" or pd.isna(row["file_2"]) or row["file_2"] == "-":
            # Single end: {sample}_trimmed.fq.gz
            targets.append(os.path.join(config["output_dir"], "1.trimmed_fastq", f"{sample}_trimmed.fq.gz"))
        else:
            # Paired end: {sample}_val_1.fq.gz and {sample}_val_2.fq.gz
            targets.append(os.path.join(config["output_dir"], "1.trimmed_fastq", f"{sample}_val_1.fq.gz"))
            targets.append(os.path.join(config["output_dir"], "1.trimmed_fastq", f"{sample}_val_2.fq.gz"))
    return targets

# --------------------------- heatmap.smk -----------------------------
# Helper function to get extendReads flag for merged bigwig
def get_merged_extend_reads_flag(wildcards):
    # 1. Identify a representative sample for the current condition.
    # Assuming you have a global pandas DataFrame 'samples' with a 'condition' column.
    # (If you don't use pandas, retrieve the list of samples for this condition from your config).
    
    # Get all samples belonging to this condition
    associated_samples = samples_df[samples_df['condition'] == wildcards.condition].index.tolist()
    
    if not associated_samples:
        raise ValueError(f"No samples found for condition: {wildcards.condition}")

    # 2. Pick the first sample to check if it is PE or SE.
    # (This assumes all samples within one condition are either all PE or all SE).
    representative_sample = associated_samples[0]

    if representative_sample in PE_SAMPLES_LIST:
        # PE: Calculate length from mates
        return "--extendReads" 
    else:
        # SE: Extend to fixed average fragment length
        frag_len = config.get("bigwig_fragment_length", 200)
        return f"--extendReads {frag_len}"

# Add this to rules/common.smk
def get_bams_by_condition(wildcards):
    """Helper to get all filtered BAMs for a specific condition."""
    # Assuming samples_df is defined globally in your Snakefile or common.smk
    sel_samples = samples_df[samples_df["condition"] == wildcards.condition]["name"]
    return [os.path.join(config["output_dir"], "4.filtered", f"{s}.shifted.bam") for s in sel_samples]


################################################################### UNTESTED ##################################################################

# def get_final_output(wildcards):
#     targets = []
    
#     # Iterate over dataframe rows safely
#     for idx, row in samples_df.iterrows():
#         sample = row["name"]
#         p_type = row.get("peak_type", "narrow") # defaults to narrow if missing
        
#         # Core sample-level outputs
#         targets.extend([
#             # Keep unshifted BAM for peaks/insert size/FRiP
#             os.path.join(config["output_dir"], "4.filtered", f"{sample}.final.bam.bai"),
#             # Request Shifted BAM for BigWigs/TSS/Heatmaps
#             os.path.join(config["output_dir"], "4.filtered", f"{sample}.shifted.bam.bai"),
#             os.path.join(config["output_dir"], "logs", "filter", f"{sample}.flagstat"),
#             os.path.join(config["output_dir"], "5.bigwig", f"{sample}.bw"),
#             os.path.join(config["output_dir"], "qc", "tss", f"{sample}.tss_enrichment.pdf"),
#             os.path.join(config["output_dir"], "6.peak_calls", f"{sample}_peaks.{p_type}Peak"),
#             # Move FRiP into this core loop since it applies to all samples
#             os.path.join(config["output_dir"], "qc", "frip", f"{sample}.frip.mqc.tsv")
#         ])

#     # PE only insert sizes (Unshifted)
#     for sample in PE_SAMPLES_LIST:
#        targets.append(os.path.join(config["output_dir"], "qc", "insert_size", f"{sample}.isize.pdf"))

#     # Consensus Targets (MSPC)
#     unique_groups = samples_df[["condition", "peak_type"]].drop_duplicates()
#     for _, row in unique_groups.iterrows():
#         targets.append(
#             os.path.join(config["output_dir"], "7.consensus", row["condition"], row["peak_type"], "ConsensusPeaks.bed")
#         )

#     # Merged BAMs and BigWigs by Condition
#     for condition in samples_df["condition"].unique():
#         targets.extend([
#             os.path.join(config["output_dir"], "4.filtered", f"{condition}.merged.bam.bai"),
#             os.path.join(config["output_dir"], "5.bigwig", "merged", f"{condition}.merged.bw")
#         ])

#     # Heatmaps
#     if "heatmap_comparisons" in config:
#         for comparison_name in config["heatmap_comparisons"]:
#             targets.append(os.path.join(config["output_dir"], "8.heatmap", f"{comparison_name}_Heatmap.pdf"))

#     # DiffBind 
#     if "diffbind_contrasts" in config and config["diffbind_contrasts"]:
#         # 1. Global outputs that are always created once per diffbind run
#         targets.extend([
#             os.path.join(config["output_dir"], "8.diffbind", "Global_PCA.pdf"),
#             os.path.join(config["output_dir"], "8.diffbind", "DiffBind_object.rds")
#         ])
        
#         # 2. Contrast-specific outputs (CSV and Volcano plot)
#         for contrast in config["diffbind_contrasts"].keys():
#             targets.extend([
#                 os.path.join(config["output_dir"], "8.diffbind", f"{contrast}_DiffPeaks.csv"),
#                 os.path.join(config["output_dir"], "8.diffbind", f"{contrast}_Volcano.pdf")
#             ])

#     # MultiQC Report
#     targets.append(os.path.join(config["output_dir"], "qc", "multiqc_report.html"))

#     return targets

def get_final_output(wildcards):
    targets = []
    
    # Iterate over dataframe rows safely
    for idx, row in samples_df.iterrows():
        sample = row["name"]
        p_type = row.get("peak_type", "narrow") # defaults to narrow if missing
        
        # Core sample-level outputs
        targets.extend([
            # Keep unshifted BAM for peaks/insert size/FRiP
            os.path.join(config["output_dir"], "4.filtered", f"{sample}.final.bam.bai"),
            # Request Shifted BAM for BigWigs/TSS/Heatmaps
            os.path.join(config["output_dir"], "4.filtered", f"{sample}.shifted.bam.bai"),
            os.path.join(config["output_dir"], "logs", "filter", f"{sample}.flagstat"),
            os.path.join(config["output_dir"], "5.bigwig", f"{sample}.bw"),
            os.path.join(config["output_dir"], "qc", "tss", f"{sample}.tss_enrichment.pdf"),
            os.path.join(config["output_dir"], "6.peak_calls", f"{sample}_peaks.{p_type}Peak"),
            # Move FRiP into this core loop since it applies to all samples
            os.path.join(config["output_dir"], "qc", "frip", f"{sample}.frip.mqc.tsv")
        ])

    # PE only insert sizes (Unshifted)
    for sample in PE_SAMPLES_LIST:
       targets.append(os.path.join(config["output_dir"], "qc", "insert_size", f"{sample}.isize.pdf"))

    # Consensus Targets (MSPC)
    unique_groups = samples_df[["condition", "peak_type"]].drop_duplicates()
    for _, row in unique_groups.iterrows():
        targets.append(
            os.path.join(config["output_dir"], "7.consensus", row["condition"], row["peak_type"], "ConsensusPeaks.bed")
        )

    # Merged BAMs and BigWigs by Condition
    for condition in samples_df["condition"].unique():
        targets.extend([
            os.path.join(config["output_dir"], "4.filtered", f"{condition}.merged.bam.bai"),
            os.path.join(config["output_dir"], "5.bigwig", "merged", f"{condition}.merged.bw")
        ])

    # Heatmaps
    if "heatmap_comparisons" in config:
        for comparison_name in config["heatmap_comparisons"]:
            targets.append(os.path.join(config["output_dir"], "8.heatmap", f"{comparison_name}_Heatmap.pdf"))

    # DiffBind Dynamically Requested Outputs
    if "diffbind_contrasts" in config and config["diffbind_contrasts"]:
        contrasts = list(config["diffbind_contrasts"].keys())
        out_dir = config["output_dir"]
        
        # This safely creates the paths for all 5 file types across all contrasts!
        targets.extend(expand(
            os.path.join(out_dir, "8.diffbind", "{contrast}_{ext}"),
            contrast=contrasts,
            ext=["DiffPeaks.csv", "UP.bed", "DOWN.bed", "Volcano.pdf", "PCA.pdf"]
        ))

    # MultiQC Report
    targets.append(os.path.join(config["output_dir"], "qc", "multiqc_report.html"))

    return targets