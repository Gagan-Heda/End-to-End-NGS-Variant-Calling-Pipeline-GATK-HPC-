#!/bin/bash
#SBATCH -M teach
#SBATCH -A hugen2072-2025s
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=GOH21@pitt.edu
#SBATCH -c 24
#SBATCH -t 96:00:00
set -ve

module load fastqc/0.11.9
module load cutadapt/2.10
module load fastqc/0.11.9
module load gatk/4.5.0.0
module load gcc/8.2.0 bwa/0.7.17 samtools/1.21

cd $SLURM_SCRATCH

#Step - 1 : FASTQC
fastqc $SLURM_SUBMIT_DIR/p5/p5_1.fastq.gz -t 48 --outdir=$SLURM_SCRATCH
fastqc $SLURM_SUBMIT_DIR/p5/p5_2.fastq.gz -t 48 --outdir=$SLURM_SCRATCH

#Step - 2 : CUTADAPTS
cutadapt -j 0 -m 10 -q 20 $SLURM_SUBMIT_DIR/p5/p5_1.fastq.gz $SLURM_SUBMIT_DIR/p5/p5_2.fastq.gz \
-a AGATCGGAAGAG -A AGATCGGAAGAG \
-o $SLURM_SCRATCH/p5_1_trimmed.fastq.gz -p $SLURM_SCRATCH/p5_2_trimmed.fastq.gz

#Step - 3 : FASTQC
fastqc $SLURM_SCRATCH/p5_1_trimmed.fastq.gz -t 48 --outdir=$SLURM_SUBMIT_DIR
fastqc $SLURM_SCRATCH/p5_2_trimmed.fastq.gz -t 48 --outdir=$SLURM_SUBMIT_DIR

#Step - 4 : ALIGNMENT
bwa mem -t 48 \
$SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
$SLURM_SCRATCH/p5_1_trimmed.fastq.gz \
$SLURM_SCRATCH/p5_2_trimmed.fastq.gz \
-R "@RG\tID:P5\tLB:P5\tSM:P5\tPL:ILLUMINA" | \
samtools view -bh | samtools sort > $SLURM_SUBMIT_DIR/p5_12.bam
samtools index $SLURM_SUBMIT_DIR/p5_12.bam

bwa mem -t 48 \
$SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
$SLURM_SCRATCH/p5_1_trimmed.fastq.gz \
$SLURM_SCRATCH/p5_2_trimmed.fastq.gz \
-R "@RG\tID:P5\tLB:P5\tSM:P5\tPL:ILLUMINA" | \
samtools view -bh | head -n 10 > $SLURM_SUBMIT_DIR/first10_alignments.txt
samtools index $SLURM_SUBMIT_DIR/p5_12.bam

run_on_exit(){
    if [ ! -f "$SLURM_SUBMIT_DIR/p5_12.bam" ]; then
        echo "Script failed before completion — restoring project.bam"
        cp p5_12.bam* $SLURM_SUBMIT_DIR/
    fi
}

trap run_on_exit EXIT

cp $SLURM_SUBMIT_DIR/p5_12.bam* . && rm $SLURM_SUBMIT_DIR/p5_12.bam*

#QUALITY CONTROL
#Step - 5 : MARK FOR DUPLICATION
gatk MarkDuplicatesSpark -I $SLURM_SCRATCH/p5_12.bam \
-O $SLURM_SUBMIT_DIR/p5_12_dupsmarked.bam

if [ $? -eq 0 ]; then
    echo "MarkDuplicatesSpark finished successfully."
else
    echo "MarkDuplicatesSpark failed."
    cp p5_12.bam* $SLURM_SUBMIT_DIR/
fi

run_on_exit(){
    if [ ! -f "$SLURM_SUBMIT_DIR/p5_12_dupsmarked.bam" ]; then    
        cp $SLURM_SCRATCH/p5_12_dupsmarked.bam* $SLURM_SUBMIT_DIR/
    fi
}

trap run_on_exit EXIT

cp $SLURM_SUBMIT_DIR/p5_12_dupsmarked.bam* $SLURM_SCRATCH/ && rm $SLURM_SUBMIT_DIR/p5_12_dupsmarked.bam*

#Step - 6 : BASE RECALIBRATION
gatk BaseRecalibrator -I $SLURM_SCRATCH/p5_12_dupsmarked.bam \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-known-sites $SLURM_SUBMIT_DIR/p5/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-known-sites $SLURM_SUBMIT_DIR/p5/dbsnp_146.hg38.vcf.gz \
-O $SLURM_SCRATCH/p5_BQSR.table

if [ $? -eq 0 ]; then
    echo "Base Recalibration finished successfully."
else
    echo "Recalibration failed."
    cp p5_12_dupsmarked.bam* $SLURM_SUBMIT_DIR/
fi

run_on_exit(){
    if [ ! -f "$SLURM_SUBMIT_DIR/p5_BQSR.table" ]; then
        cp $SLURM_SCRATCH/p5_BQSR.table  $SLURM_SUBMIT_DIR/ && cp $SLURM_SCRATCH/p5_12_dupsmarked.bam* $SLURM_SUBMIT_DIR/
    fi
}

trap run_on_exit EXIT
echo "Starting Apply BQSR"

#Step - 7 : APPLY CALIBRATION
gatk ApplyBQSR -R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-I $SLURM_SCRATCH/p5_12_dupsmarked.bam \
--bqsr-recal-file $SLURM_SCRATCH/p5_BQSR.table \
-O $SLURM_SUBMIT_DIR/p5_12_dupsmarked_cleaned.bam


if [ $? -eq 0 ]; then
    echo "Calibration applied successfully."
else
    echo "Calibration application failed."
    cp $SLURM_SCRATCH/p5_dupsmarked.bam* $SLURM_SUBMIT_DIR/ && cp $SLURM_SCRATCH/p5_BQSR* $SLURM_SUBMIT_DIR/ 
fi

#Step - 8 : ALIGNMENT AND DEPTH STATISTICS
samtools flagstat $SLURM_SUBMIT_DIR/p5_12_dupsmarked_cleaned.bam > $SLURM_SUBMIT_DIR/p5_alignment_statistics.out
samtools depth $SLURM_SUBMIT_DIR/p5_12_dupsmarked_cleaned.bam \
| awk '{if ($3 >= 0) {total_bases += $3}} END {print total_bases / NR}' \
> $SLURM_SUBMIT_DIR/p5_avg_read_depth.out

#GENOTYPING
#Step - 9 : HAPLOTYPE CALLER 
gatk HaplotypeCaller -R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-I $SLURM_SUBMIT_DIR/p5_12_dupsmarked_cleaned.bam \
-O p5_genotypes.g.vcf.gz \
-ERC GVCF \
-OVI \
--native-pair-hmm-threads 48


if [ $? -eq 0 ]; then
    echo "HaplotypeCaller ran successfully."
    cp $SLURM_SCRATCH/p5_genotypes.g.vcf.gz $SLURM_SUBMIT_DIR/ && cp $SLURM_SUBMIT_DIR/p5_12_dupsmarked_cleaned.bam* $SLURM_SCRATCH
else
    echo "HaplotypeCaller failed."
fi


run_on_exit(){
    if [ ! -f "$SLURM_SUBMIT_DIR/p5_12_dupsmarked_cleaned.bam*" ]; then
        echo "Script failed before completion — restoring project.bam"
        cp p5_12_dupsmarked_cleaned* $SLURM_SUBMIT_DIR/
    fi
}
trap run_on_exit EXIT

#Step - 10 : CONVERSION OF GVCF TO VCF FILE
gatk GenotypeGVCFs -R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V $SLURM_SUBMIT_DIR/p5_genotypes.g.vcf.gz \
--java-options "-XX:ParallelGCThreads=48" \
-O p5_genotypes.vcf.gz

if [ $? -eq 0 ]; then
    echo "Genotyping ran successfully."
    cp $SLURM_SCRATCH/p5_genotypes.vcf.gz $SLURM_SUBMIT_DIR/
else
    echo "Genotyping failed."
fi


#GENOTYPE QUALITY CONTROL

#Step - 11 : SITE INFORMATION
gatk MakeSitesOnlyVcf -I $SLURM_SUBMIT_DIR/p5_genotypes.vcf.gz -O p5_sites_only.vcf.gz

#Step - 12a : VARIANT RECALIBRATOR [INDEL]
gatk VariantRecalibrator -mode INDEL \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V p5_sites_only.vcf.gz \
-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
-resource:mills,known=false,training=true,truth=true,prior=12 \
$SLURM_SUBMIT_DIR/p5/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=2 \
$SLURM_SUBMIT_DIR/p5/dbsnp_146.hg38.vcf.gz \
--java-options "-XX:ParallelGCThreads=48" \
-O p5_indels.recal \
--tranches-file p5_indels.tranches \

if [ $? -eq 0 ]; then
    echo "Variant Recalibration -INDEL ran successfully."
    cp $SLURM_SCRATCH/p5_indels* $SLURM_SUBMIT_DIR/
else
    echo "Variant Recalibration -INDEL failed."
fi  

#Step - 12b : VARIANT RECALIBRATOR [SNP]
gatk VariantRecalibrator -mode SNP \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V p5_sites_only.vcf.gz \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
-resource:hapmap,known=false,training=true,truth=true,prior=15 \
$SLURM_SUBMIT_DIR/p5/hapmap_3.3.hg38.vcf.gz \
-resource:omni,known=false,training=true,truth=true,prior=12 \
$SLURM_SUBMIT_DIR/p5/1000G_omni2.5.hg38.vcf.gz \
-resource:1000G,known=false,training=true,truth=false,prior=10 \
$SLURM_SUBMIT_DIR/p5/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=7 \
$SLURM_SUBMIT_DIR/p5/dbsnp_146.hg38.vcf.gz \
-O p5_snps.recal \
--tranches-file p5_snps.tranches

#Step - 12c : APPLY THE RECALIBRATION STATISTICS[INDEL] 
gatk ApplyVQSR -mode INDEL \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V $SLURM_SUBMIT_DIR/p5_genotypes.vcf.gz \
--recal-file p5_indels.recal \
--tranches-file p5_indels.tranches \
--truth-sensitivity-filter-level 99.0 \
--create-output-variant-index true \
--java-options "-XX:ParallelGCThreads=48" \
-O p5_genotypes_indelqc.vcf.gz


#Step - 12d : APPLY THE RECALIBRATION STATISTICS[SNP]
gatk ApplyVQSR -mode SNP \
-R $SLURM_SUBMIT_DIR/p5/Homo_sapiens_assembly38.fasta \
-V p5_genotypes_indelqc.vcf.gz \
--recal-file p5_snps.recal \
--tranches-file p5_snps.tranches \
--truth-sensitivity-filter-level 99.0 \
--create-output-variant-index true \
-O $SLURM_SUBMIT_DIR/p5.vcf.gz

#Step - 13 : STATISTICS
gatk CollectVariantCallingMetrics -I $SLURM_SUBMIT_DIR/p5.vcf.gz \
--DBSNP $SLURM_SUBMIT_DIR/p5/dbsnp_146.hg38.vcf.gz \
-O $SLURM_SUBMIT_DIR/p5_genotype_metrics