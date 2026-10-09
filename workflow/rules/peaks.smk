# --- Helper Functions ---
def get_control_bam(wildcards):
    try:
        control_name = samples_df.loc[wildcards.sample, "control"]
        if control_name != "-" and pd.notna(control_name):
            return os.path.join(config["output_dir"], "4.filtered", f"{control_name}.final.bam")
    except KeyError:
        pass
    return []

def get_control_flag(wildcards, input):
    # Only add -c if control input exists
    if hasattr(input, "control") and input.control:
        return f"-c {input.control}"
    return ""

def get_macs3_params(wildcards):
    seq_type = samples_df.loc[wildcards.sample, "seq_type"]
    if seq_type == "paired":
        return "-f BAMPE --keep-dup all"
    else:
        # Shift for Tn5 cut site (Single End only)
        return "-f BAM --nomodel --shift -75 --extsize 150 --keep-dup all"

# --- Rules ---

rule macs3_call_narrow_peaks:
    input:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.final.bam"),
        control = get_control_bam
    output:
        peaks = os.path.join(config["output_dir"], "6.peak_calls", "{sample}_peaks.narrowPeak"),
        xls = os.path.join(config["output_dir"], "6.peak_calls", "{sample}_peaks.xls"),
        summits = os.path.join(config["output_dir"], "6.peak_calls", "{sample}_summits.bed")
    wildcard_constraints:
        sample = NARROW_SAMPLES
    log:
        os.path.join(config["output_dir"], "logs", "peaks", "{sample}.macs3.log")
    benchmark:
        os.path.join(config["output_dir"], "benchmarks", "peaks", "{sample}_macs3.tsv")
    params:
        outdir = lambda w: os.path.join(config["output_dir"], "6.peak_calls"),
        name = lambda w: w.sample,
        genome_size = config.get("macs_gsize", "hs"),
        control_cmd = get_control_flag,
        mode_opts = get_macs3_params
    container:
        "docker://quay.io/biocontainers/macs3:3.0.3--py312h4711d71_0"
    shell:
        """
        # -n {params.name} results in output: {params.outdir}/{params.name}_peaks.narrowPeak
        macs3 callpeak \
            -t {input.bam} \
            {params.control_cmd} \
            -g {params.genome_size} \
            -n {params.name} \
            {params.mode_opts} \
            -q 0.01 --call-summits \
            --outdir {params.outdir} \
            2> {log}
        """

rule macs3_call_broad_peaks:
    input:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.final.bam"),
        control = get_control_bam
    output:
        peaks = os.path.join(config["output_dir"], "6.peak_calls", "{sample}_peaks.broadPeak"),
        gapped = os.path.join(config["output_dir"], "6.peak_calls", "{sample}_peaks.gappedPeak"),
        xls = os.path.join(config["output_dir"], "6.peak_calls", "{sample}_peaks.xls")
    wildcard_constraints:
        sample = BROAD_SAMPLES
    log:
        os.path.join(config["output_dir"], "logs", "peaks", "{sample}.macs3.log")
    benchmark:
        os.path.join(config["output_dir"], "benchmarks", "peaks", "{sample}_macs3.tsv")
    params:
        outdir = lambda w: os.path.join(config["output_dir"], "6.peak_calls"),
        name = lambda w: w.sample,
        genome_size = config.get("macs_gsize", "hs"),
        control_cmd = get_control_flag,
        mode_opts = get_macs3_params
    container:
        "docker://quay.io/biocontainers/macs3:3.0.3--py312h4711d71_0"
    shell:
        """
        # -n {params.name} results in output: {params.outdir}/{params.name}_peaks.broadPeak
        # Note: --broad-cutoff 0.1 is standard, adjust if needed
        macs3 callpeak \
            -t {input.bam} \
            {params.control_cmd} \
            -g {params.genome_size} \
            -n {params.name} \
            {params.mode_opts} \
            --broad --broad-cutoff 0.1 \
            --outdir {params.outdir} \
            2> {log}
        """


