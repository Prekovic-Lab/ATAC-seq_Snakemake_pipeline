# Function to determine the extendReads flag based on sample type
def get_extend_reads_flag(wildcards):
    if wildcards.sample in PE_SAMPLES_LIST:
        # PE: Calculate length from mates (exact)
        return "--extendReads" 
    else:
        # SE: Extend to fixed average fragment length
        # Default to 200bp if not specified in config
        frag_len = config.get("bigwig_fragment_length", 200)
        return f"--extendReads {frag_len}"

# BAM to BigWig Conversion (deepTools)
rule bam_to_bigwig:
    input:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.shifted.bam"),
        bai = os.path.join(config["output_dir"], "4.filtered", "{sample}.shifted.bam.bai")
    output:
        bw = os.path.join(config["output_dir"], "5.bigwig", "{sample}.bw")
    log:
        os.path.join(config["output_dir"], "logs", "qc", "{sample}_bam_to_bigwig.log")
    benchmark:
        os.path.join(config["output_dir"], "benchmarks", "qc", "{sample}_bam_to_bigwig.tsv")
    resources:
        mem_mb = config.get('bigwig_mem_mb', 16000),
        runtime = config.get('bigwig_runtime', 240)
    params:
        binsize = 10,
        effective_genome_size = config.get("effective_genome_size"),
        ignoreForNormalization = config.get("bigwig_ignore_for_normalization", "chrM"),
        extend_reads_cmd = get_extend_reads_flag
    threads: 8
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
    shell:
        "bamCoverage --bam {input.bam} \
            --outFileName {output.bw} \
            --binSize {params.binsize} \
            --normalizeUsing RPGC \
            --ignoreDuplicates \
            --effectiveGenomeSize {params.effective_genome_size} \
            --ignoreForNormalization {params.ignoreForNormalization} \
            --minMappingQuality 1 \
            {params.extend_reads_cmd} \
            --numberOfProcessors {threads} \
            2> {log}"


# Helper to get all BAMs for a specific condition
def get_bams_by_condition(wildcards):
    # Filter for the specific condition requested
    # wildcards.condition comes from the output filename
    sel_samples = samples_df[samples_df["condition"] == wildcards.condition]["name"]
    return [os.path.join(config["output_dir"], "4.filtered", f"{s}.final.bam") for s in sel_samples]

# Merge Replicate BAMs
############################ Make output temp() once it's confirmed #######################
# rule merge_condition_bams:
#     input:
#         get_bams_by_condition
#     output:
#         bam = os.path.join(config["output_dir"], "4.filtered", "{condition}.merged.bam")
#     log:
#         os.path.join(config["output_dir"], "logs", "merge", "{condition}_merge.log")
#     resources:
#         mem_mb = 4000,
#         runtime = 120,
#         tmpdir = config.get("tmp_dir", "$TMPDIR") 
#     container:
#         "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
#     threads: 8
#     shell:
#         """
#         samtools merge -@ {threads} -o {output.bam} {input} 2> {log}
#         """
################################################ TROUBLESHOOTING RULE ABOVE ####################################################
################################################ IF IT WORKS, DELETE THE ONE ABOVE #############################################
rule merge_condition_bams:
    input:
        bams = get_bams_by_condition
    output:
        bam = os.path.join(config["output_dir"], "4.filtered", "{condition}.merged.bam")
    log:
        os.path.join(config["output_dir"], "logs", "merge", "{condition}_merge.log")
    resources:
        mem_mb = 4000,
        runtime = 120,
        tmpdir = config.get("tmp_dir", "$TMPDIR") 
    container:
        "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    threads: 8
    resources:
        mem_mb = 4000,
        runtime = 120,
        tmpdir = config.get("tmp_dir", "$TMPDIR") 
    run:
        # If there is only 1 file, just copy it to the merged name
        if len(input.bams) == 1:
            shell("cp {input.bams[0]} {output.bam} 2> {log}")
        # If there are multiple files, use samtools merge
        else:
            shell("samtools merge -@ {threads} -o {output.bam} {input.bams} 2> {log}")

# Create Merged BigWig (Same settings as your individual rule)
rule merged_bam_to_bigwig:
    input:
        bam = os.path.join(config["output_dir"], "4.filtered", "{condition}.merged.bam"),
        bai = os.path.join(config["output_dir"], "4.filtered", "{condition}.merged.bam.bai")
    output:
        bw = os.path.join(config["output_dir"], "5.bigwig", "merged", "{condition}.merged.bw")
    log:
        os.path.join(config["output_dir"], "logs", "qc", "{condition}_merged_bw.log")
    resources:
        mem_mb = config.get('bigwig_mem_mb', 16000),
        runtime = config.get('bigwig_runtime', 240)
    threads: 8
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
    params:
        binsize = 10,
        effective_genome_size = config.get("effective_genome_size"),
        ignoreForNormalization = config.get("bigwig_ignore_for_normalization", "chrM"),
        extend_reads_cmd = get_merged_extend_reads_flag
    shell:
        """
        bamCoverage --bam {input.bam} \
            --outFileName {output.bw} \
            --binSize {params.binsize} \
            --normalizeUsing RPGC \
            --ignoreDuplicates \
            --effectiveGenomeSize {params.effective_genome_size} \
            --ignoreForNormalization {params.ignoreForNormalization} \
            --minMappingQuality 1 \
            {params.extend_reads_cmd} \
            --numberOfProcessors {threads} \
            2> {log}
        """


############################################### NEW RULES FOR MERGING BAMS BY CONDITION ###############################################
############################################### THEO HAS TO CHECK THAT ################################################################

# rule bigwig_merged:
#     input:
#         bam = os.path.join(config["output_dir"], "4.filtered", "merged", "{condition}.merged.bam"),
#         bai = os.path.join(config["output_dir"], "4.filtered", "merged", "{condition}.merged.bam.bai")
#     output:
#         bw = os.path.join(config["output_dir"], "5.bigwig", "merged", "{condition}.merged.bw")
#     log:
#         os.path.join(config["output_dir"], "logs", "bigwig", "merged_{condition}.log")
#     threads: 16
#     container:
#         "docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
#     params:
#         eff_genome_size = config["effective_genome_size"],
#         ignore = " ".join(config["bigwig_ignore_for_normalization"])
#     shell:
#         """
#         bamCoverage -b {input.bam} \
#             -o {output.bw} \
#             --binSize 10 \
#             --normalizeUsing RPGC \
#             --effectiveGenomeSize {params.eff_genome_size} \
#             --ignoreForNormalization {params.ignore} \
#             --numberOfProcessors {threads} \
#             --extendReads 2> {log}
#         """