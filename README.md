# End-to-End NGS Variant Calling Pipeline (GATK + HPC)

A high-throughput, reproducible NGS pipeline for variant discovery following GATK Best Practices. The workflow automates quality control, alignment, post-processing, variant calling, and variant quality recalibration.

---

## Overview

Accurate variant discovery requires multiple processing and QC steps. This pipeline transforms raw FASTQ files into high-confidence variant calls while ensuring reproducibility in HPC environments.

**Key Goals:**

- Perform QC on raw sequencing reads  
- Trim adapters and low-quality bases  
- Align reads to the reference genome  
- Process alignments for downstream analysis  
- Call and genotype variants  
- Apply variant quality recalibration  
- Generate alignment and variant metrics  

---

## Workflow

### 1. Quality Control
- FastQC analysis for read quality, GC content, and sequencing artifacts  

### 2. Adapter Trimming
- Cutadapt removes adapters and low-quality bases  
- Applies minimum read length and quality thresholds  

### 3. Alignment
- BWA-MEM aligns reads to GRCh38 reference genome  
- SAM → sorted BAM conversion  
- BAM indexing for efficient access  

### 4. Post-Alignment Processing
- GATK MarkDuplicatesSpark for duplicate marking  
- Base quality score recalibration (BQSR) using known variant sites  
- Generates recalibrated BAM files  

### 5. Alignment Metrics
- Samtools for alignment statistics  
- Computes read depth, coverage, and alignment summary metrics  

### 6. Variant Calling
- GATK HaplotypeCaller generates intermediate GVCFs  
- Joint genotyping using GenotypeGVCFs  

### 7. Variant Quality Recalibration
- Sites-only VCF created  
- VQSR applied to SNPs and INDELs using multiple training datasets  

### 8. Final Variant Set
- Filtered, high-confidence variants  
- Variant-level quality metrics for downstream analysis  

---

## Technologies

- Bash (Unix scripting)  
- SLURM (HPC job scheduling)  
- FastQC  
- Cutadapt  
- BWA-MEM  
- Samtools  
- GATK  

---

## Key Features

- Fully automated, end-to-end NGS workflow  
- GATK Best Practices implementation  
- HPC-compatible parallel execution  
- Multi-step quality control and validation  
- Checkpointing and error recovery  

---

## Skills Demonstrated

- NGS data processing and variant calling  
- High-performance computing (SLURM)  
- Alignment and post-processing of sequencing reads  
- Variant quality control and recalibration  
- Automation of large-scale genomic workflows  
- Handling of large-scale sequencing datasets  

---

## Notes

- Optimized for HPC environments with SLURM  
- Intermediate files managed for storage efficiency  
- Variant recalibration uses established reference datasets  
- Error handling ensures reproducibility and recovery of intermediate outputs  

---

## Author

**Gagan Heda**
