# Insert Size Distribution (Picard)
# Checks if we have the nucleosome banding pattern (147bp, 300bp, etc.)
rule picard_insert_size:
    input:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.final.bam")
    output:
        txt = os.path.join(config["output_dir"], "qc", "insert_size", "{sample}.isize.txt"),
        pdf = os.path.join(config["output_dir"], "qc", "insert_size", "{sample}.isize.pdf")
    log:
        os.path.join(config["output_dir"], "logs", "qc", "{sample}_isize.log")
    benchmark:
        os.path.join(config["output_dir"], "benchmarks", "qc", "{sample}_isize.tsv")
    wildcard_constraints:
        sample = PE_SAMPLES
    resources:
        mem_mb = 8000,
        runtime = 30
    container:
        # Official Picard container
        "docker://broadinstitute/picard:2.27.5"
    shell:
        """
        java -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
            I={input.bam} \
            O={output.txt} \
            H={output.pdf} \
            M=0.5 \
            2> {log}
        """
# Prepare Matrix for TSS Plot (deepTools)
# Calculates coverage signal around Transcription Start Sites
rule deeptools_compute_matrix:
    input:
        bw = os.path.join(config["output_dir"], "5.bigwig", "{sample}.bw"),
        gtf = config["ref_gtf"]
    output:
        matrix = temp(os.path.join(config["output_dir"], "qc", "tss", "{sample}.matrix.gz"))
    log:
        os.path.join(config["output_dir"], "logs", "qc", "{sample}_compute_matrix.log")
    benchmark:
        os.path.join(config["output_dir"], "benchmarks", "qc", "{sample}_compute_matrix.tsv")
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 80
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
    shell:
        """
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 2000 \
            -R {input.gtf} \
            -S {input.bw} \
            --numberOfProcessors {threads} \
            --missingDataAsZero \
            -o {output.matrix} \
            2> {log}
        """

# Plot TSS Enrichment (deepTools)
# Creates the actual Heatmap/Profile image
rule deeptools_plot_tss:
    input:
        matrix = os.path.join(config["output_dir"], "qc", "tss", "{sample}.matrix.gz")
    output:
        img = os.path.join(config["output_dir"], "qc", "tss", "{sample}.tss_enrichment.pdf")
    log:
        os.path.join(config["output_dir"], "logs", "qc", "{sample}_plot_tss.log")
    benchmark:
        os.path.join(config["output_dir"], "benchmarks", "qc", "{sample}_plot_tss.tsv")
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
    shell:
        """
        plotHeatmap \
            -m {input.matrix} \
            -out {output.img} \
            --colorMap Blues \
            --whatToShow "plot, heatmap and colorbar" \
            --refPointLabel "TSS" \
            2> {log}
        """

# --- FastQC on Raw Data (PE) ---
rule fastqc_raw_pe:
    input:
        r1 = lambda w: os.path.join(config["fastq_dir"], samples_df.loc[w.sample, "file_1"]),
        r2 = lambda w: os.path.join(config["fastq_dir"], samples_df.loc[w.sample, "file_2"])
    output:
        zips = [
            os.path.join(config["output_dir"], "qc", "fastqc", "{sample}_R1_fastqc.zip"),
            os.path.join(config["output_dir"], "qc", "fastqc", "{sample}_R2_fastqc.zip")
        ],
        htmls = [
            os.path.join(config["output_dir"], "qc", "fastqc", "{sample}_R1_fastqc.html"),
            os.path.join(config["output_dir"], "qc", "fastqc", "{sample}_R2_fastqc.html")
        ]
    wildcard_constraints:
        sample = PE_SAMPLES
    threads: 4
    resources:
        mem_mb = 4000
    container:
        "docker://quay.io/biocontainers/fastqc:0.11.9--0"
    shell:
        """
        # Create a temp directory for safe naming
        mkdir -p {config[output_dir]}/qc/fastqc/temp_{wildcards.sample}
        
        # Symlink inputs to standard names (R1.fastq.gz, R2.fastq.gz)
        # This forces FastQC to output standard filenames.
        ln -s {input.r1} {config[output_dir]}/qc/fastqc/temp_{wildcards.sample}/{wildcards.sample}_R1.fastq.gz
        ln -s {input.r2} {config[output_dir]}/qc/fastqc/temp_{wildcards.sample}/{wildcards.sample}_R2.fastq.gz
        
        # Run FastQC
        fastqc {config[output_dir]}/qc/fastqc/temp_{wildcards.sample}/*.fastq.gz \
            --outdir {config[output_dir]}/qc/fastqc/ \
            --threads {threads}
            
        # Cleanup temp links
        rm -rf {config[output_dir]}/qc/fastqc/temp_{wildcards.sample}
        """

# --- FastQC on Raw Data (SE) ---
rule fastqc_raw_se:
    input:
        r1 = lambda w: os.path.join(config["fastq_dir"], samples_df.loc[w.sample, "file_1"])
    output:
        zip = os.path.join(config["output_dir"], "qc", "fastqc", "{sample}_R1_fastqc.zip"),
        html = os.path.join(config["output_dir"], "qc", "fastqc", "{sample}_R1_fastqc.html")
    wildcard_constraints:
        sample = SE_SAMPLES
    threads: 2
    resources:
        mem_mb = 4000
    container:
        "docker://quay.io/biocontainers/fastqc:0.11.9--0"
    shell:
        """
        mkdir -p {config[output_dir]}/qc/fastqc/temp_{wildcards.sample}
        
        ln -s {input.r1} {config[output_dir]}/qc/fastqc/temp_{wildcards.sample}/{wildcards.sample}_R1.fastq.gz
        
        fastqc {config[output_dir]}/qc/fastqc/temp_{wildcards.sample}/*.fastq.gz \
            --outdir {config[output_dir]}/qc/fastqc/ \
            --threads {threads}
            
        rm -rf {config[output_dir]}/qc/fastqc/temp_{wildcards.sample}
        """

# --- FastQC on Filtered BAMs ---
rule fastqc_filtered:
    input:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.final.bam")
    output:
        zip = os.path.join(config["output_dir"], "qc", "fastqc_filtered", "{sample}.final_fastqc.zip"),
        html = os.path.join(config["output_dir"], "qc", "fastqc_filtered", "{sample}.final_fastqc.html")
    log:
        os.path.join(config["output_dir"], "logs", "qc", "{sample}_fastqc_filtered.log")
    threads: 2
    resources:
        mem_mb = 4000
    container:
        "docker://quay.io/biocontainers/fastqc:0.11.9--0"
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.zip})
        
        fastqc {input.bam} \
            --outdir $(dirname {output.zip}) \
            --threads {threads} \
            2> {log}
        """

# Helper to gather Raw FastQC files based on sample type
def get_raw_fastqc_input(wildcards):
    files = []
    for sample in samples_df["name"]:
        # Check if PE or SE based on your dataframe logic
        is_paired = (samples_df.loc[sample, "seq_type"] == "paired")
        
        files.append(os.path.join(config["output_dir"], "qc", "fastqc", f"{sample}_R1_fastqc.zip"))
        if is_paired:
            files.append(os.path.join(config["output_dir"], "qc", "fastqc", f"{sample}_R2_fastqc.zip"))
    return files

def get_peak_xls_for_multiqc(wildcards):
    files = []
    for sample in samples_df["name"]:
        # The .xls file is generated by BOTH narrow and broad modes.
        # Since get_final_output ensures the correct rule runs, 
        # the .xls file will exist.
        files.append(os.path.join(config["output_dir"], "6.peak_calls", f"{sample}_peaks.xls"))
    return files

# # multiqc rule
# rule multiqc:
#     input:
#         # --- QC Metrics ---
#         raw_fastqc = get_raw_fastqc_input,
#         trim_logs = expand(os.path.join(config["output_dir"], "logs", "trimming", "{sample}.log"), sample=samples_df["name"]),
#         align_logs = expand(os.path.join(config["output_dir"], "logs", "align", "{sample}.log"), sample=samples_df["name"]),
#         dedup_metrics = expand(os.path.join(config["output_dir"], "3.dedup", "{sample}.dedup_metrics.txt"), sample=samples_df["name"]),
#         filter_stats = expand(os.path.join(config["output_dir"], "logs", "filter", "{sample}.flagstat"), sample=samples_df["name"]),
#         insert_size = expand(os.path.join(config["output_dir"], "qc", "insert_size", "{sample}.isize.txt"), sample=samples_df["name"]),
#         filtered_fastqc = expand(os.path.join(config["output_dir"], "qc", "fastqc_filtered", "{sample}.final_fastqc.zip"), sample=samples_df["name"]),
#         peak_xls = get_peak_xls_for_multiqc,

#         # --- Benchmarks (for Resource Usage) ---
#         # We explicitly list them so Snakemake creates them before running MultiQC
#         bench_trim = expand(os.path.join(config["output_dir"], "benchmarks", "trimming", "{sample}.tsv"), sample=samples_df["name"]),
#         bench_align = expand(os.path.join(config["output_dir"], "benchmarks", "align", "{sample}.tsv"), sample=samples_df["name"]),
#         bench_dedup = expand(os.path.join(config["output_dir"], "benchmarks", "dedup", "{sample}.tsv"), sample=samples_df["name"]),
#         bench_filter = expand(os.path.join(config["output_dir"], "benchmarks", "filter", "{sample}.tsv"), sample=samples_df["name"]),
#         bench_peakcalling = expand(os.path.join(config["output_dir"], "benchmarks", "peaks", "{sample}_macs3.tsv"), sample=samples_df["name"])
#     output:
#         report = os.path.join(config["output_dir"], "qc", "multiqc_report.html")
#     log:
#         os.path.join(config["output_dir"], "logs", "qc", "multiqc.log")
#     container:
#         "docker://quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"
#     shell:
#         """
#         # We pass the root output directory. 
#         # MultiQC will recursively find logs AND benchmarks (if they are .tsv or .txt).
#         multiqc {config[output_dir]} \
#             --outdir $(dirname {output.report}) \
#             --filename multiqc_report.html \
#             --force \
#             --comment "ATAC-seq Pipeline Report" \
#             2> {log}
#         """

rule frip_score:
    input:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.final.bam"),
        # Grab the narrowPeak (you can adjust this if you mix broad/narrow)
        peaks = os.path.join(config["output_dir"], "6.peak_calls", "{sample}_peaks.narrowPeak")
    output:
        txt = os.path.join(config["output_dir"], "qc", "frip", "{sample}.frip.mqc.tsv")
    log:
        os.path.join(config["output_dir"], "logs", "qc", "{sample}_frip.log")
    threads: 2
    container:
        "docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
    shell:
        """
        # Count total reads
        TOTAL=$(samtools view -@ {threads} -c {input.bam})
        
        # Count reads overlapping peaks
        IN_PEAKS=$(samtools view -@ {threads} -c -L {input.peaks} {input.bam})
        
        # Calculate fraction using awk (outputs a format MultiQC loves)
        awk -v t=$TOTAL -v p=$IN_PEAKS -v s={wildcards.sample} 'BEGIN {{
            frip = p/t;
            print "Sample\\tTotal_Reads\\tReads_In_Peaks\\tFRiP";
            print s"\\t"t"\\t"p"\\t"frip;
        }}' > {output.txt} 2> {log}
        """

rule multiqc:
    input:
        # --- Trigger Upstream Rules via their Main Outputs ---
        # 1. FastQC requires raw data
        raw_fastqc = get_raw_fastqc_input,
        
        # 2. Dedup metrics require the Dedup rule (which requires Align rule)
        dedup_metrics = expand(os.path.join(config["output_dir"], "3.dedup", "{sample}.dedup_metrics.txt"), sample=samples_df["name"]),
        
        # 3. Filter stats require the Filter rule
        filter_stats = expand(os.path.join(config["output_dir"], "logs", "filter", "{sample}.flagstat"), sample=samples_df["name"]),
        
        # 4. Filtered FastQC requires the Filter rule
        filtered_fastqc = expand(os.path.join(config["output_dir"], "qc", "fastqc_filtered", "{sample}.final_fastqc.zip"), sample=samples_df["name"]),
        
        # 5. Peak Calling requires MACS/Genrich
        peak_xls = get_peak_xls_for_multiqc,
        
        # 6. Insert sizes
        insert_size = expand(os.path.join(config["output_dir"], "qc", "insert_size", "{sample}.isize.txt"), sample=samples_df["name"])

        # REMOVE explicit calls to logs and benchmarks here!
        # The rules generating the files above will automatically generate the logs/benchmarks 
        # as side effects before this rule runs.

    output:
        report = os.path.join(config["output_dir"], "qc", "multiqc_report.html")
    log:
        os.path.join(config["output_dir"], "logs", "qc", "multiqc.log")
    container:
        "docker://quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"
    shell:
        """
        # MultiQC scans the directory structure. 
        # Since the inputs above force the pipeline to run alignment/trimming/dedup,
        # the log files WILL exist on the disk by the time this command runs.
        
        multiqc {config[output_dir]} \
            --outdir $(dirname {output.report}) \
            --filename multiqc_report.html \
            --force \
            --comment "ATAC-seq Pipeline Report" \
            2> {log}
        """