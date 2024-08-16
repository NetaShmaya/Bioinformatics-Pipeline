#!/usr/bin/env nextflow

// Define parameters
params.fastq = ''
params.ref_fasta = ''
params.ref_index = ''
params.outdir = './results'
params.annovar_db = ''
params.annovar_exe = ''

// Process definition for FastQC
process FastQC {
    tag "FastQC on ${fastq_file}"

    publishDir params.outdir, mode: 'copy'

    input:
    path fastq_file

    output:
    path "*.zip", emit: zip
    path "*.html", emit: html

    script:
    """
    fastqc ${fastq_file}
    """
}

// Process definition for BwaMem
process BwaMem {
    tag "BWA MEM on ${fastq_file}"

    publishDir params.outdir, mode: 'copy'

    input:
    path fastq_file
    path ref_fasta
    path ref_index_files

    output:
    tuple path("*.bam"), path("*.bai"), emit: bam_and_bai

    script:
    """
    bwa mem ${ref_fasta} ${fastq_file} | samtools view -Sb - | samtools sort -o ${fastq_file.simpleName}.sorted.bam -
    samtools index ${fastq_file.simpleName}.sorted.bam
    """
}

// Process definition for GATK HaplotypeCaller
process HaplotypeCaller {
    tag "GATK HaplotypeCaller on ${bam_file}"

    publishDir params.outdir, mode: 'copy'

    container 'broadinstitute/gatk:latest'

    input:
    tuple path(bam_file), path(bai_file)
    path ref_fasta
    path ref_fasta_fai
    path ref_dict

    output:
    path "*.vcf", emit: vcf

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${bam_file} \
        -O ${bam_file.simpleName}.vcf
    """
}

// Process definition for ANNOVAR annotation
process ANNOVAR {
    tag "ANNOVAR annotation on ${vcf_file}"

    publishDir params.outdir, mode: 'copy'

    input:
    path vcf_file

    output:
    path "*.annovar_output", emit: annovar_output

    script:
    """
    ${params.annovar_exe} \
        --buildver hg19 \
        --out ${vcf_file.simpleName} \
        --remove \
        --protocol refGene,cytoBand,snp138 \
        --operation g,r,f \
        --nastring . \
        ${vcf_file} \
        ${params.annovar_db}
    """
}

// Main workflow
workflow {
    // Create channels for the input files
    fastq_file = file(params.fastq)
    ref_fasta = file(params.ref_fasta)
    ref_fasta_fai = file("${params.ref_fasta}.fai")
    ref_dict = file("${params.ref_fasta.tokenize('.').init().join('.')}.dict")
    ref_index_files = Channel.fromPath(params.ref_index)

    // Run FastQC
    FastQC(fastq_file)

    // Run BwaMem
    bam_and_bai = BwaMem(fastq_file, ref_fasta, ref_index_files.collect())

    // Run GATK HaplotypeCaller
    vcf = HaplotypeCaller(bam_and_bai, ref_fasta, ref_fasta_fai, ref_dict)

    // Run ANNOVAR
    ANNOVAR(vcf)
}

