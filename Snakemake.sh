rule all:
    input:
        "output/variants.vcf"

rule fastqc:
    input:
        "raw_data/{sample}.R1.fastq.gz",
        "raw_data/{sample}.R2.fastq.gz"
    output:
        "qc/{sample}_R1_fastqc.html",
        "qc/{sample}_R2_fastqc.html"
    shell:
        "fastqc -o qc {input}"

rule cutadapt:
    input:
        "raw_data/{sample}.R1.fastq.gz",
        "raw_data/{sample}.R2.fastq.gz"
    output:
        "clean_data/{sample}.R1.fastq.gz",
        "clean_data/{sample}.R2.fastq.gz"
    shell:
        "cutadapt -o clean_data/{sample}.R1.fastq.gz -p clean_data/{sample}.R2.fastq.gz {input}"

rule download_reference:
    output:
        "reference/chr13.fa.gz",
        "reference/chr17.fa.gz"
    shell:
        "wget -P reference ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.13.fa.gz"
        "wget -P reference ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.17.fa.gz"

rule extract_reference:
    input:
        "reference/chr13.fa.gz",
        "reference/chr17.fa.gz"
    output:
        "reference/chr13.fa",
        "reference/chr17.fa"
    shell:
        "gunzip -k {input}"

rule concatenate_reference:
    input:
        "reference/chr13.fa",
        "reference/chr17.fa"
    output:
        "reference/BRCA1_2.fa"
    shell:
        "cat {input} > {output}"

rule index_genome:
    input:
        "reference/BRCA1_2.fa"
    output:
        "reference/BRCA1_2.fa.fai"
    shell:
        "bwa index {input}"

rule map_reads:
    input:
        "clean_data/{sample}.R1.fastq.gz",
        "clean_data/{sample}.R2.fastq.gz"
    output:
        "mapped_reads/{sample}.sam"
    shell:
        "bwa mem -t 4 reference/BRCA1_2.fa {input} > {output}"

rule convert_sam_to_bam:
    input:
        "mapped_reads/{sample}.sam"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "samtools view -b -S -o {output} {input}"

rule remove_duplicates:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "mapped_reads/{sample}_dedup.bam"
    shell:
        "samtools rmdup {input} {output}"

rule call_variants:
    input:
        "mapped_reads/{sample}_dedup.bam"
    output:
        "variants/{sample}.bcf"
    shell:
        "samtools mpileup -g -f reference/BRCA1_2.fa {input} | bcftools call -vmO b -o {output}"

rule convert_bcf_to_vcf:
    input:
        "variants/{sample}.bcf"
    output:
        "variants/{sample}.vcf"
    shell:
        "bcftools view {input} > {output}"

rule filter_variants:
    input:
        "variants/{sample}.vcf"
    output:
        "filtered_variants/{sample}.vcf"
    shell:
        "java -jar SnpSift.jar filter '((QUAL > 20) & (AF > 0.1))' {input} > {output}"

rule annotate_variants:
    input:
        "filtered_variants/{sample}.vcf"
    output:
        "annotated_variants/{sample}.vcf"
    shell:
        "vep -i {input} -o {output}"

