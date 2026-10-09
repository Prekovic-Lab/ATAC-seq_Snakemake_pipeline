# workflow/rules/align.smk

# ALIGN PAIRED-END (Outputs Temp SAM)
rule align_pe:
    input:
        r1 = os.path.join(config["output_dir"], "1.trimmed_fastq", "{sample}_val_1.fq.gz"),
        r2 = os.path.join(config["output_dir"], "1.trimmed_fastq", "{sample}_val_2.fq.gz")
    output:
        # We mark this as temp() so Snakemake deletes it after sorting is done
        sam = temp(os.path.join(config["output_dir"], "2.mapped", "{sample}.unsorted.sam"))
    log:
        config["output_dir"] + "/logs/align/{sample}_bowtie2.log"
    benchmark:
        config["output_dir"] + "/benchmarks/align/{sample}_bowtie2.tsv"
    wildcard_constraints:
        sample = PE_SAMPLES
    threads: 12
    resources:
        mem_mb = config.get("align_mem_mb", 32000),
        runtime = config.get("align_runtime", 360)
    container:
        # Official Bowtie2 container
        "docker://quay.io/biocontainers/bowtie2:2.5.4--he96a11b_7"
    params:
        idx_prefix = config["bowtie2_index"],
        extra = "--very-sensitive -X 2000"
    shell:
        """
        bowtie2 \
            -p {threads} \
            {params.extra} \
            --rg-id {wildcards.sample} \
            --rg "SM:{wildcards.sample}" \
            --rg "PL:ILLUMINA" \
            -x {params.idx_prefix} \
            -1 {input.r1} \
            -2 {input.r2} \
            -S {output.sam} \
            2> {log}
        """

################################################ HAS YET TO BE TESTED ################################################
# ALIGN SINGLE-END (Outputs Temp SAM)
rule align_se:
    input:
        r1 = os.path.join(config["output_dir"], "1.trimmed_fastq", "{sample}_trimmed.fq.gz")
    output:
        sam = temp(os.path.join(config["output_dir"], "2.mapped", "{sample}.unsorted.sam"))
    log:
        config["output_dir"] + "/logs/align/{sample}.log"
    benchmark:
        config["output_dir"] + "/benchmarks/align/{sample}.tsv"
    wildcard_constraints:
        sample = SE_SAMPLES
    threads: 12
    resources:
        mem_mb = config.get("align_mem_mb", 32000),
        runtime = config.get("align_runtime", 240)
    container:
        "docker://quay.io/biocontainers/bowtie2:2.5.4--he96a11b_7"
    params:
        idx_prefix = config["bowtie2_index"],
        extra = "--very-sensitive"
    shell:
        """
        bowtie2 \
            -p {threads} \
            {params.extra} \
            --rg-id {wildcards.sample} \
            --rg "SM:{wildcards.sample}" \
            --rg "PL:ILLUMINA" \
            -x {params.idx_prefix} \
            -U {input.r1} \
            -S {output.sam} \
            2> {log}
        """
################################################ END OF UNTESTED CODE ################################################

# SORT SAM TO BAM
rule samtools_sort:
    input:
        sam = os.path.join(config["output_dir"], "2.mapped", "{sample}.unsorted.sam")
    output:
        bam = temp(os.path.join(config["output_dir"], "2.mapped", "{sample}.sorted.bam"))
    log:
        config["output_dir"] + "/logs/align/{sample}_sort.log"
    benchmark:
        config["output_dir"] + "/benchmarks/align/{sample}_sort.tsv"
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 60
    container:
        # Official Samtools container
        "docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
    shell:
        """
        samtools sort -@ {threads} -O BAM -o {output.bam} {input.sam} 2> {log}
        """

# INDEX BAM
rule samtools_index:
    input:
        "{filepath}.bam"
    output:
        "{filepath}.bam.bai"
    container:
        "docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
    shell:
        "samtools index {input}"