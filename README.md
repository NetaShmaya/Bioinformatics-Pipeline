# Bioinformatics-Pipeline

This repository contains a Nextflow pipeline that processes FASTQ files, performs quality control, alignment, variant calling, and annotation, leading to an annotated VCF file.

## Pipeline Overview

1. **Quality Control (QC)**: FastQC is used to assess the quality of the raw FASTQ files.
2. **Alignment**: BWA is used to align the reads to the reference genome, and the resulting alignments are sorted and converted to BAM format.
3. **Indexing**: SAMtools is used to index the sorted BAM file.
4. **Variant Calling**: GATK HaplotypeCaller is used to call variants and generate a VCF file.
5. **Annotation**: ANNOVAR is used to annotate the VCF file with functional information.

## Prerequisites

Before running the pipeline, ensure you have the following installed:
- [Nextflow](https://www.nextflow.io/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [SAMtools](http://www.htslib.org/)
- [GATK](https://gatk.broadinstitute.org/)
- [ANNOVAR](http://annovar.openbioinformatics.org/)

## Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/Bioinformatics-Pipeline.git
   cd Bioinformatics-Pipeline
2. Update the config.nf file with the correct paths to your input files, reference genome, ANNOVAR database, and executables.
3. Ensure the run_pipeline.sh script has execution permissions:

   ```chmod +x run_pipeline.sh```

## Execution
To execute the pipeline, simply run:
  
   ```./run_pipeline.sh```

## Pipeline Details

1. **QC**:
   - Command: `fastqc $fastq`
   - Description: This step uses FastQC to assess the quality of the input FASTQ files.

2. **Alignment and BAM Conversion**:
   - Command: `bwa mem $REF $fastq_file | samtools sort -o $output.bam -`
   - Description: BWA is used to align the reads to the reference genome, and the alignments are then sorted and converted to BAM format using SAMtools.

3. **Indexing**:
   - Command: `samtools index $output.bam`
   - Description: This step indexes the sorted BAM file to allow for efficient querying during downstream analysis.

4. **Variant Calling**:
   - Command: `gatk HaplotypeCaller -R $REF -I $output.bam -O $variants.vcf`
   - Description: GATK HaplotypeCaller is used to identify variants in the aligned reads and output them in VCF format.

5. **Annotation**:
   - Command: `table_annovar.pl $vcf ${params.annovar}/humandb/ -buildver hg38 -out annotated -remove -protocol refGene,cytoBand,clinvar_20220320 -operation g,r,f -nastring . -vcfinput`
   - Description: ANNOVAR annotates the VCF file with functional information from various databases, including refGene, cytoBand, and ClinVar.

   
