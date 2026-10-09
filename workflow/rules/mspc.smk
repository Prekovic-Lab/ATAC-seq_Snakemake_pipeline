# --- Target Generation ---

# We only want one consensus file per Condition + PeakType combination
CONSENSUS_TARGETS = []
unique_groups = samples_df[["condition", "peak_type"]].drop_duplicates()

for _, row in unique_groups.iterrows():
    CONSENSUS_TARGETS.append(
        os.path.join(config["output_dir"], "7.consensus", row["condition"], row["peak_type"], "ConsensusPeaks.bed")
    )

# --- Helper Functions for MSPC ---

def get_mspc_inputs(wildcards):
    """
    Finds all replicates for a specific condition and peak type.
    Example: If condition=ATCC and peak_type=narrow, returns [R1_ATCC...narrowPeak, R2_ATCC...narrowPeak]
    """
    # Filter the samples_df for the specific condition and peak_type
    sub_df = samples_df[
        (samples_df["condition"] == wildcards.condition) & 
        (samples_df["peak_type"] == wildcards.peak_type)
    ]
    
    # Determine the extension (narrowPeak or broadPeak)
    ext = "narrowPeak" if wildcards.peak_type == "narrow" else "broadPeak"
    
    # Return the list of expected file paths from the previous MACS3 step
    return [
        os.path.join(config["output_dir"], "6.peak_calls", f"{row.name}_peaks.{ext}")
        for _, row in sub_df.iterrows()
    ]

rule mspc_consensus:
    input:
        peaks = get_mspc_inputs
    output:
        consensus = os.path.join(config["output_dir"], "7.consensus", "{condition}", "{peak_type}", "ConsensusPeaks.bed")
    log:
        os.path.join(config["output_dir"], "logs", "mspc", "{condition}_{peak_type}.log")
    params:
        outdir = lambda w: os.path.join(config["output_dir"], "7.consensus", w.condition, w.peak_type),
        config_file = "/usr/local/bin/mspc_config.json"
    container: 
        "docker://theo46/mspc:v2" 
    resources:
        runtime = config.get('mspc_runtime', 10)
    run:
        if len(input.peaks) == 1:
            shell("cp {input.peaks[0]} {output.consensus} > {log} 2>&1")
        else:
            shell(
                """
                mspc -i {input.peaks} \
                --output {params.outdir} \
                -r bio \
                -w 1e-4 \
                -s 1e-8 \
                -p {params.config_file} \
                --excludeHeader > {log} 2>&1
                """
            )