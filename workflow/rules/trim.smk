# workflow/rules/trim.smk

rule trim_pe:
    input:
        get_fastqs
    output:
        # Trim Galore automatically appends _val_1.fq.gz and _val_2.fq.gz
        # We must match that pattern here.
        r1   = os.path.join(config["output_dir"], "1.trimmed_fastq", "{sample}_val_1.fq.gz"),
        r2   = os.path.join(config["output_dir"], "1.trimmed_fastq", "{sample}_val_2.fq.gz"),
        # Trim Galore also produces report files
        rep1 = temp(os.path.join(config["output_dir"], "1.trimmed_fastq", "{sample}_1.fastq.gz_trimming_report.txt")),
        rep2 = temp(os.path.join(config["output_dir"], "1.trimmed_fastq", "{sample}_2.fastq.gz_trimming_report.txt"))   
    log:
        config["output_dir"] + "/logs/trimming/{sample}.log"
    benchmark:
        config["output_dir"] + "/benchmarks/trimming/{sample}.tsv"
    threads: 6
    resources:
        mem_mb  = config.get('trim_mem_mb', 8000),
        runtime = config.get('trim_runtime', 120),
        tmpdir  = config["tmp_dir"]
    wildcard_constraints:
        sample = PE_SAMPLES
    container:
        # This container includes TrimGalore, Cutadapt, and FastQC
        "docker://quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0"
    params:
        outdir = os.path.join(config["output_dir"], "1.trimmed_fastq")
    shell:
        """
        # 1. Run Trim Galore 
        trim_galore \
            --paired \
            --gzip \
            --cores {threads} \
            --output_dir {params.outdir} \
            --basename {wildcards.sample} \
            {input[0]} {input[1]} \
            2> {log}

        # 2. Rename Report Files
        # TrimGalore names reports based on INPUT filenames but Snakemake expects wildcard-based names 
        
        # Capture input filenames
        f1=$(basename {input[0]})
        f2=$(basename {input[1]})

        # Move the generated reports to the expected output paths
        mv {params.outdir}/${{f1}}_trimming_report.txt {output.rep1}
        mv {params.outdir}/${{f2}}_trimming_report.txt {output.rep2}
        """

################################################ HAS YET TO BE TESTED ################################################`
rule trim_se:
    input:
        get_fastqs
    output:
        r1   = temp(os.path.join(config["output_dir"], "1.trimmed_fastq", "{sample}_trimmed.fq.gz")),
        rep  = os.path.join(config["output_dir"], "1.trimmed_fastq", "{sample}.fastq.gz_trimming_report.txt")
    log:
        config["output_dir"] + "/logs/trimming/{sample}.log"
    benchmark:
        config["output_dir"] + "/benchmarks/trimming/{sample}.tsv"
    threads: 6
    resources:
        mem_mb  = config.get('trim_mem_mb', 8000), 
        runtime = config.get('trim_runtime', 120),
        tmpdir  = config["tmp_dir"]
    wildcard_constraints:
        sample = SE_SAMPLES
    container:
        "docker://quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0"
    params:
        outdir = os.path.join(config["output_dir"], "1.trimmed_fastq")
    shell:
        """
        # 1. Run Trim Galore
        trim_galore \
            --gzip \
            --cores {threads} \
            --output_dir {params.outdir} \
            --basename {wildcards.sample} \
            {input[0]} \
            2> {log}
        
        # 2. Rename Report File
        # TrimGalore names reports based on INPUT filenames but Snakemake expects wildcard-based names 

        # Capture input filename
        f1=$(basename {input[0]})

        # Move the generated report to the expected output path
        mv {params.outdir}/${{f1}}_trimming_report.txt {output.rep}
        """