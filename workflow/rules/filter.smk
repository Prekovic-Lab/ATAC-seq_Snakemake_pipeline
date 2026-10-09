rule filter_bam:
    input:
        bam = os.path.join(config["output_dir"], "3.dedup", "{sample}.dedup.bam"),
        # We don't strictly need the input index for streaming, but it's good practice
        blacklist = config["blacklist_bed"]
    output:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.final.bam"),
        stats = os.path.join(config["output_dir"], "logs", "filter", "{sample}.flagstat")
    log:
        os.path.join(config["output_dir"], "logs", "filter", "{sample}.log")
    benchmark:
        os.path.join(config["output_dir"], "benchmarks", "filter", "{sample}.tsv")
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 100
    params:
        samtools_flags = lambda w: "-f 2" if w.sample in PE_SAMPLES_LIST else "-F 4"
    container:
        # Check it to change maybe to a quay.io image?
        "docker://danhumassmed/samtools-bedtools:1.0.2"
    shell:
        """
        # Filter using Samtools (View) | Bedtools (Intersect)
        samtools view -h -u -q 30 {params.samtools_flags} -F 2048 -e 'rname != "chrM"' {input.bam} \
        | bedtools intersect -v -abam stdin -b {input.blacklist} \
        | samtools sort -O BAM -o {output.bam} - 2> {log}
        
        # Calculate stats (QC)
        samtools flagstat {output.bam} > {output.stats}
        """

############################################### NEW RULES FOR MERGING BAMS BY CONDITION ###############################################
############################################### THEO HAS TO CHECK THAT ################################################################

# workflow/rules/filter.smk

rule tn5_shift_bam:
    input:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.final.bam"),
        # alignmentSieve requires the BAM to be indexed
        bai = os.path.join(config["output_dir"], "4.filtered", "{sample}.final.bam.bai")
    output:
        # We mark this as temp() so it deletes itself after sorting
        bam = temp(os.path.join(config["output_dir"], "4.filtered", "{sample}.shifted.unsorted.bam"))
    log:
        os.path.join(config["output_dir"], "logs", "filter", "{sample}_tn5_shift.log")
    threads: 8
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0"
    shell:
        """
        # --ATACshift shifts +4bp on the positive strand and -5bp on the negative strand
        alignmentSieve \
            --numberOfProcessors {threads} \
            --ATACshift \
            -b {input.bam} \
            -o {output.bam} \
            2> {log}
        """

rule sort_shifted_bam:
    input:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.shifted.unsorted.bam")
    output:
        bam = os.path.join(config["output_dir"], "4.filtered", "{sample}.shifted.bam")
    log:
        os.path.join(config["output_dir"], "logs", "filter", "{sample}_sort_shifted.log")
    threads: 4
    container:
        "docker://danhumassmed/samtools-bedtools:1.0.2"
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        """
        samtools sort -@ {threads} -O BAM -o {output.bam} {input.bam} 2> {log}
        """

# rule merge_bams_by_condition:
#     input:
#         bams = get_bams_by_condition
#     output:
#         bam = os.path.join(config["output_dir"], "4.filtered", "merged", "{condition}.merged.bam")
#     log:
#         os.path.join(config["output_dir"], "logs", "filter", "merge_{condition}.log")
#     threads: 8
#     container:
#         "docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
#     shell:
#         """
#         # If there is only one file, copy it. Otherwise, merge.
#         if [ $(echo {input.bams} | wc -w) -eq 1 ]; then
#             cp {input.bams} {output.bam}
#         else
#             samtools merge -@ {threads} -o {output.bam} {input.bams} 2> {log}
#         fi
#         """

# rule index_merged_bams:
#     input:
#         bam = os.path.join(config["output_dir"], "4.filtered", "merged", "{condition}.merged.bam")
#     output:
#         bai = os.path.join(config["output_dir"], "4.filtered", "merged", "{condition}.merged.bam.bai")
#     threads: 4
#     container:
#         "docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
#     shell:
#         """
#         samtools index -@ {threads} {input.bam}
#         """