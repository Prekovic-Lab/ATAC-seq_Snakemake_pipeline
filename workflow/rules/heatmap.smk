# workflow/rules/heatmap.smk
import os
import pandas as pd

# Helper Functions
def get_comparison_conditions(comparison_name):
    """Retrieve the list of conditions for a specific comparison from config."""
    try:
        return config["heatmap_comparisons"][comparison_name]["conditions"]
    except KeyError:
        raise ValueError(f"Comparison '{comparison_name}' not found in config['heatmap_comparisons']")

def get_plot_type(comparison_name):
    """Returns 'merged' or 'replicates' based on config (default: replicates)"""
    return config["heatmap_comparisons"][comparison_name].get("plot_type", "replicates")

def get_kmeans_k(wildcards):
    """Get the specific K value for this comparison."""
    return config["heatmap_comparisons"][wildcards.comparison]["k"]

def get_peaks_for_comparison(wildcards):
    """Get Consensus Peaks ONLY for conditions in this comparison."""
    conditions = get_comparison_conditions(wildcards.comparison)
    paths = []
    for cond in conditions:
        # Get peak type for this condition from samples_df
        # We grab the first occurrence of the condition to determine the peak_type
        p_type = samples_df[samples_df["condition"] == cond]["peak_type"].iloc[0]
        paths.append(
             os.path.join(config["output_dir"], "7.consensus", cond, p_type, "ConsensusPeaks.bed")
        )
    return paths

def get_bigwigs_for_comparison(wildcards):
    """Get BigWigs based on whether we want Replicates or Merged tracks."""
    conditions = get_comparison_conditions(wildcards.comparison)
    plot_type = get_plot_type(wildcards.comparison)
    
    paths = []
    
    if plot_type == "merged":
        # Return one bigwig per condition
        for cond in conditions:
            paths.append(os.path.join(config["output_dir"], "5.bigwig", "merged", f"{cond}.merged.bw"))
    else:
        # Return all replicates
        sub_df = samples_df[samples_df["condition"].isin(conditions)]
        sub_df = sub_df.sort_values(["condition", "replicate"])
        for sample in sub_df["name"]:
            paths.append(os.path.join(config["output_dir"], "5.bigwig", f"{sample}.bw"))
            
    return paths

def get_labels_for_comparison(wildcards):
    """Get labels based on plot type for the heatmap columns."""
    conditions = get_comparison_conditions(wildcards.comparison)
    plot_type = get_plot_type(wildcards.comparison)
    
    if plot_type == "merged":
        return conditions # e.g., ["ATCC", "EML4"]
    else:
        # Return sample names e.g., ["R1_ATCC", "R2_ATCC", ...]
        sub_df = samples_df[samples_df["condition"].isin(conditions)]
        sub_df = sub_df.sort_values(["condition", "replicate"])
        return sub_df["name"].tolist()

def get_bams_by_condition(wildcards):
    """Helper to get all BAMs for a specific condition (used for merging)."""
    sel_samples = samples_df[samples_df["condition"] == wildcards.condition]["name"]
    return [os.path.join(config["output_dir"], "4.filtered", f"{s}.final.bam") for s in sel_samples]

# Heatmap Rules

# Create a Specific Consensus Set for this Comparison
rule create_comparison_consensus:
    input:
        peaks = get_peaks_for_comparison
    output:
        master = os.path.join(config["output_dir"], "8.heatmap", "{comparison}", "Merged_Peaks.bed")
    log:
        os.path.join(config["output_dir"], "logs", "heatmap", "{comparison}_merge_peaks.log")
    container:
        "docker://quay.io/biocontainers/bedtools:2.31.0--h468198e_0"
    shell:
        """
        cat {input.peaks} | sort -k1,1 -k2,2n | bedtools merge -i - > {output.master} 2> {log}
        """

# Compute Matrix for this Comparison
rule compute_matrix_comparison:
    input:
        regions = os.path.join(config["output_dir"], "8.heatmap", "{comparison}", "Merged_Peaks.bed"),
        bigwigs = get_bigwigs_for_comparison
    output:
        matrix = os.path.join(config["output_dir"], "8.heatmap", "{comparison}", "matrix.gz")
    log:
        os.path.join(config["output_dir"], "logs", "heatmap", "{comparison}_compute_matrix.log")
    threads: 16
    resources:
        mem_mb = 64000,
        runtime = 120
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
    params:
        labels = get_labels_for_comparison,
        before = 5000,
        after = 5000
    shell:
        """
        computeMatrix reference-point \
            --referencePoint center \
            --regionsFileName {input.regions} \
            --scoreFileName {input.bigwigs} \
            --outFileName {output.matrix} \
            --samplesLabel {params.labels} \
            --beforeRegionStartLength {params.before} \
            --afterRegionStartLength {params.after} \
            --skipZeros \
            --numberOfProcessors {threads} \
            2> {log}
        """

# Plot Heatmap for this Comparison
rule plot_heatmap_comparison:
    input:
        matrix = os.path.join(config["output_dir"], "8.heatmap", "{comparison}", "matrix.gz")
    output:
        plot = os.path.join(config["output_dir"], "8.heatmap", "{comparison}_Heatmap.pdf"),
        matrix_out = os.path.join(config["output_dir"], "8.heatmap", "{comparison}_matrix_values.tab")
    log:
        os.path.join(config["output_dir"], "logs", "heatmap", "{comparison}_plot.log")
    resources:
        mem_mb = 64000,
        runtime = 60
    threads: 1
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
    params:
        k = get_kmeans_k,
        color = "Reds"
    shell:
        """
        plotHeatmap \
            -m {input.matrix} \
            -out {output.plot} \
            --outFileSortedRegions {output.matrix_out} \
            --kmeans {params.k} \
            --colorMap {params.color} \
            --whatToShow 'heatmap and colorbar' \
            --interpolationMethod 'bilinear' \
            2> {log}
        """