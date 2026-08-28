#!/bin/bash
#SBATCH --job-name=handrun_DATASET
#SBATCH --partition=compute
#SBATCH --output= DATASET
#SBATCH --time=48:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=16

export PATH= 

REF=…./pigeon_ref_v39
GTF=$REF/gencode.v39.annotation.sorted.gtf
GENOME=$REF/human_GRCh38_no_alt_analysis_set.fasta
CAGE=$REF/refTSS_v3.3_human_coordinate.hg38.sorted.bed
POLYA=$REF/polyA.list.txt
INTROPOLIS=$REF/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv

MAS_ADAPTERS=/data/.../sc_isoseq_run/inputs/mas16_primers.fasta
PRIMERS=/data/..../sc_isoseq_run/inputs/10x_5kit_primers.fasta
BC5P=/data/.../sc_isoseq_run/inputs/737K_august_2016.txt.gz
CCS_BAM=/data/…./dataset_PacBio_10x5/original/m84014_240128_063621_s2.hifi_reads.bcM0003.bam

OUTDIR=/data/handrun_DATASET
mkdir -p $OUTDIR
cd $OUTDIR

echo "[$(date)] STEP 1: skera split"
skera split -j $SLURM_CPUS_PER_TASK $CCS_BAM $MAS_ADAPTERS segmented.bam
echo "skera exit: $?"

echo "[$(date)] STEP 2: lima"
lima -j $SLURM_CPUS_PER_TASK --per-read --isoseq segmented.bam $PRIMERS demuxed.bam
echo "lima exit: $?"
LIMA_OUT=$(ls demuxed.5p--3p.bam 2>/dev/null | head -1)
[ -f "$LIMA_OUT" ] || { echo "FATAL: lima output not found"; ls -la demuxed*; exit 1; }

echo "[$(date)] STEP 3: isoseq tag (5' design)"
isoseq tag -j $SLURM_CPUS_PER_TASK $LIMA_OUT flt.bam --design 16B-10U-10X-T
echo "tag exit: $?"

echo "[$(date)] STEP 4: isoseq refine"
isoseq refine -j $SLURM_CPUS_PER_TASK flt.bam $PRIMERS fltnc.bam --require-polya
echo "refine exit: $?"

echo "[$(date)] STEP 5: isoseq correct (5' barcode whitelist)"
isoseq correct -j $SLURM_CPUS_PER_TASK --barcodes $BC5P --filter failing fltnc.bam fltnc.corrected.bam
echo "correct exit: $?"

echo "[$(date)] STEP 6: sort by CB tag"
samtools sort -t CB -@ $SLURM_CPUS_PER_TASK fltnc.corrected.bam -o fltnc.corrected.sorted.bam
echo "sort exit: $?"

echo "[$(date)] STEP 7: isoseq groupdedup"
isoseq groupdedup -j $SLURM_CPUS_PER_TASK fltnc.corrected.sorted.bam dedup.bam
echo "groupdedup exit: $?"

echo "[$(date)] STEP 8: pbmm2 align"
pbmm2 align --preset ISOSEQ --sort -j $SLURM_CPUS_PER_TASK $GENOME dedup.bam mapped.bam
echo "align exit: $?"

echo "[$(date)] STEP 9: isoseq collapse"
isoseq collapse -j $SLURM_CPUS_PER_TASK --max-fuzzy-junction 4 mapped.bam collapsed.gff
echo "collapse exit: $?"

echo "[$(date)] STEP 10: pigeon prepare + classify"
pigeon prepare collapsed.gff
pigeon classify collapsed.sorted.gff $GTF $GENOME \
    --fl collapsed.abundance.txt \
    --cage-peak $CAGE --poly-a $POLYA --coverage $INTROPOLIS \
    -d $OUTDIR -o handrun_5p_classification
echo "classify exit: $?"

echo "[$(date)] STEP 11: pigeon filter"
pigeon filter handrun_5p_classification_classification.txt -i collapsed.sorted.gff -j 1
echo "filter exit: $?"

echo "[$(date)] STEP 12: pigeon make-seurat"
pigeon make-seurat --dedup dedup.fasta -g collapsed.group.txt -j $SLURM_CPUS_PER_TASK -d . \
    handrun_5p_classification_classification.filtered_lite_classification.txt
echo "make-seurat exit: $?"

echo "[$(date)] DONE. Category breakdown:"
F=$(ls *filtered_lite_classification.txt 2>/dev/null | head -1)
echo "TOTAL:      $(tail -n +2 $F | wc -l)  (SMRT Link result: 761,833)"
echo "FSM:        $(awk -F'\t' '$6=="full-splice_match"' $F | wc -l)"
echo "ISM:        $(awk -F'\t' '$6=="incomplete-splice_match"' $F | wc -l)"
echo "NIC:        $(awk -F'\t' '$6=="novel_in_catalog"' $F | wc -l)"
echo "NNC:        $(awk -F'\t' '$6=="novel_not_in_catalog"' $F | wc -l)"
echo "intergenic: $(awk -F'\t' '$6=="intergenic"' $F | wc -l)"