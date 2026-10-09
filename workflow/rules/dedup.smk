# Mark duplicates rule
rule mark_duplicates:
    input:
        bam = os.path.join(config["output_dir"], "2.mapped", "{sample}.sorted.bam")
    output:
        bam = temp(os.path.join(config["output_dir"], "3.dedup", "{sample}.dedup.bam")),
        metrics = os.path.join(config["output_dir"], "3.dedup", "{sample}.dedup_metrics.txt")
    log:
        config["output_dir"] + "/logs/dedup/{sample}.log"
    benchmark:
        config["output_dir"] + "/benchmarks/dedup/{sample}.tsv"
    params:
        tmpdir = config["tmp_dir"]
    resources:
        mem_mb = config.get('dedup_mem_mb', 30000),
        runtime = config.get('dedup_runtime', 240)
    threads: 8
    container:
          "docker://broadinstitute/gatk:4.6.1.0"
    shell:
        """
        # Create the specific temp directory first!
        mkdir -p {params.tmpdir}/{wildcards.sample}

        # Mark duplicates using GATK Spark for multi-threaded performance
        gatk --java-options "-Xmx24g" MarkDuplicatesSpark \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --remove-sequencing-duplicates false \
            --create-output-bam-index false \
            --spark-master local[{threads}] \
            --read-name-regex null \
            --tmp-dir {params.tmpdir}/{wildcards.sample} \
            2> {log}

        # Clean up temp files after job finishes to save space
        rm -rf {params.tmpdir}/{wildcards.sample}
        """