#!/bin/bash

# ========================================================================================
# --- Set up and run tapeworm RNAseq analysis
# --- for RNAseq workshop July 2017
# ========================================================================================

cd /gpfs/scratch/cxb585/tapeworm_RNAseq

mkdir tapeworm-rnaseq
cd tapeworm-rnaseq

# ----------------------------------------------------------------------------------------
# --- Download genomes
# ----------------------------------------------------------------------------------------

mkdir genomes
cd genomes
sh ../scripts/download_tapeworm_genomes.sh
sh ../scripts/download_tapeworm_annotations.sh

# ----------------------------------------------------------------------------------------
# --- Build index
# ----------------------------------------------------------------------------------------

module load bowtie/2.2.3

###export PATH=$PATH:/gpfs/cyberstar/ghp3/Bergey/bin/bowtie2-2.3.1/

mkdir reports

bowtie2-build genomes/Tsolium_Mexico_v1/Tsolium_Mexico_v1.fa Tsolium_Mexico_v1 \
    &> reports/bowtie2-build.log

mv Tsolium_Mexico_v1.* genomes/Tsolium_Mexico_v1/

# ----------------------------------------------------------------------------------------
# --- Link reads
# ----------------------------------------------------------------------------------------

for FQ in `ls $GRP/tapeworm_RNAseq/Taenia_RNAseq_2017-03/GP02272017/*.fastq.gz`; do
    IND=`basename $FQ | sed -e "s/_S.*//"`
    READ=`basename $FQ | sed -e "s/.*_\(R.\).*/\1/"`
    ln -s $FQ data/Taenia_PR${IND}_$READ.fastq.gz
done

# ----------------------------------------------------------------------------------------
# --- Trim reads
# ----------------------------------------------------------------------------------------

# Can be called with pbs/trim_reads.pb

# This module load actually does nothing.
module load Trimmomatic/0.32

TRIMMOMATIC=/usr/global/Trimmomatic/0.32/trimmomatic-0.32.jar

ADAPTERS=`dirname $TRIMMOMATIC`/adapters/TruSeq3-PE.fa

for THIS_IND in `ls data/*R1.fastq.gz`; do
    THIS_IND=`basename $THIS_IND | sed -e "s/_R1.*//"`

    echo "Trimming for individual $THIS_IND";

    java -jar $TRIMMOMATIC PE -phred33 \
        data/${THIS_IND}_R1.fastq.gz \
        data/${THIS_IND}_R2.fastq.gz \
        data/${THIS_IND}_TRIM_R1.fastq.gz \
        data/${THIS_IND}_trim_unpaired_R1.fastq.gz \
        data/${THIS_IND}_TRIM_R2.fastq.gz \
        data/${THIS_IND}_trim_unpaired_R2.fastq.gz \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done

# ----------------------------------------------------------------------------------------
# --- Align reads
# ----------------------------------------------------------------------------------------

# Bowtie module doesn't work and top version installed can't find bowtie2
# even if it's in the PATH
#export PATH=$PATH:/gpfs/cyberstar/ghp3/Bergey/bin/bowtie2-2.3.1/
module load bowtie/2.2.3
export PATH=$PATH:/gpfs/cyberstar/ghp3/Bergey/bin/tophat-2.1.1.Linux_x86_64

INDEX="genomes/Tsolium_Mexico_v1/Tsolium_Mexico_v1"
GTF="genomes/Tsolium_Mexico_v1/taenia_solium.PRJNA170813.WBPS8.canonical_geneset.gtf"

mkdir results/bam

for THIS_IND in `ls data/*_TRIM_R1.fastq.gz`; do

    THIS_IND=`basename $THIS_IND | sed -e "s/_TRIM_R1.*//"`

    READ_1="data/${THIS_IND}_TRIM_R1.fastq.gz"
    READ_2="data/${THIS_IND}_TRIM_R2.fastq.gz"

    IND_DIR=results/tophat/$THIS_IND
    mkdir -p $IND_DIR

    tophat2 \
        -o $IND_DIR \
        -g 1 \
        --library-type fr-firststrand \
        -G $GTF \
        -x 1 \
        -p 5 \
        $INDEX \
        $READ_1 $READ_2 &> $IND_DIR/run_tophat.log

    mv $IND_DIR/accepted_hits.bam results/bam/$THIS_IND.bam
done

# ----------------------------------------------------------------------------------------
# --- Select reads by strand and index all BAMs
# ----------------------------------------------------------------------------------------

module load samtools/1.2

for THIS_IND in `ls data/*_TRIM_R1.fastq.gz`; do

    THIS_IND=`basename $THIS_IND | sed -e "s/_TRIM_R1.*//"`
    IN_BAM=results/bam/$THIS_IND.bam

    OUT_MIN=results/bam/${THIS_IND}_min.bam
    OUT_PLUS=results/bam/${THIS_IND}_plus.bam

    samtools  view -f99  -hb $IN_BAM > ${OUT_MIN}_99.bam
    samtools  view -f147 -hb $IN_BAM > ${OUT_MIN}_147.bam
    samtools  view -f83  -hb $IN_BAM > ${OUT_PLUS}_83.bam
    samtools  view -f163 -hb $IN_BAM > ${OUT_PLUS}_163.bam

    samtools merge -h ${OUT_MIN}_99.bam  ${OUT_MIN}  ${OUT_MIN}_99.bam  ${OUT_MIN}_147.bam
    samtools merge -h ${OUT_PLUS}_83.bam ${OUT_PLUS} ${OUT_PLUS}_83.bam ${OUT_PLUS}_163.bam

    rm -f ${OUT_MIN}_99.bam ${OUT_MIN}_147.bam ${OUT_PLUS}_83.bam ${OUT_PLUS}_163.bam

    samtools index $IN_BAM
    samtools index $OUT_MIN
    samtools index $OUT_PLUS

done

# ----------------------------------------------------------------------------------------
# --- Index genome FASTA
# ----------------------------------------------------------------------------------------

samtools faidx genomes/Tsolium_Mexico_v1/Tsolium_Mexico_v1.fa
cut -f1,2 genomes/Tsolium_Mexico_v1/Tsolium_Mexico_v1.fa.fai > \
    genomes/Tsolium_Mexico_v1/Tsolium_Mexico_v1.fa.index

# ----------------------------------------------------------------------------------------
# --- Create bigwigs for both strands
# ----------------------------------------------------------------------------------------

module load python/2.7.11.ucs4
wget https://raw.githubusercontent.com/dnanexus/rseqc/master/rseqc/scripts/bam2wig.py

# Call to bam2wig function requires option, chrom_file, which isn't given in bam2wig.py
# def bamTowig(self, outfile, chrom_sizes, chrom_file,
#              skip_multi=True, strand_rule=None, WigSumFactor=None, q_cut=30)
# So fix it:
sed -e "s/\(\tobj\.bamTowig(.*chromSizes,\).*\(skip_multi=.*\)/\1 chrom_file=options.chromSize, \2/" bam2wig.py > bam2wig.fix.py

module load ucsc_tools/312

mkdir results/bwig

for THIS_IND in `ls data/*_TRIM_R1.fastq.gz`; do

    THIS_IND=`basename $THIS_IND | sed -e "s/_TRIM_R1.*//"`
    IN_BAM=results/bam/$THIS_IND.bam

    IN_MIN=results/bam/${THIS_IND}_min.bam
    IN_PLUS=results/bam/${THIS_IND}_plus.bam

    OUT_MIN_BW=${IN_MIN/bam$/bw}
    OUT_PLUS_BW=${IN_PLUS/bam$/bw}

    CHR_SIZES=genomes/Tsolium_Mexico_v1/Tsolium_Mexico_v1.fa.index

    python bam2wig.fix.py -i $IN_PLUS \
        -s $CHR_SIZES \
        -o results/bwig/${THIS_IND}_plus &> ${OUT_PLUS_BW}.log
    wigToBigWig -clip results/bwig/${THIS_IND}_plus.wig $CHR_SIZES ${OUT_PLUS_BW} \
        2>&1 >> ${OUT_PLUS_BW}.log
    rm -f results/bwig/${THIS_IND}_plus.wig

    python bam2wig.fix.py -i $IN_MIN \
        -s $CHR_SIZES \
        -o results/bwig/${THIS_IND}_min &> ${OUT_MIN_BW}.log
    wigToBigWig -clip results/bwig/${THIS_IND}_min.wig $CHR_SIZES ${OUT_MIN_BW} \
        2>&1 >> ${OUT_MIN_BW}.log
    rm -f results/bwig/${THIS_IND}_min.wig

done

# ----------------------------------------------------------------------------------------
# --- Run cufflinks to find novel transcripts
# ----------------------------------------------------------------------------------------

module load cufflinks/2.2.1

GTF="genomes/Tsolium_Mexico_v1/taenia_solium.PRJNA170813.WBPS8.canonical_geneset.gtf"

mkdir -p results/cufflinks

for THIS_IND in `ls data/*_TRIM_R1.fastq.gz`; do

    THIS_IND=`basename $THIS_IND | sed -e "s/_TRIM_R1.*//"`
    IN_BAM=results/bam/$THIS_IND.bam
    cufflinks -g $GTF -p 5 --library-type fr-firststrand \
        -o results/cufflinks/$THIS_IND $IN_BAM &> results/cufflinks/$THIS_IND.log

done

# ----------------------------------------------------------------------------------------
# --- Compare novel and known transcripts
# ----------------------------------------------------------------------------------------

module load cufflinks/2.2.1

GTF="genomes/Tsolium_Mexico_v1/taenia_solium.PRJNA170813.WBPS8.canonical_geneset.gtf"

mkdir -p results/cuffmerge
ls -1 results/cufflinks/*/transcripts.gtf > results/cuffmerge/assembly.txt

cuffmerge \
    -o results/cuffmerge \
    -g $GTF \
    --keep-tmp \
    -s genomes/Tsolium_Mexico_v1/Tsolium_Mexico_v1.fa \
    -p 5 \
    results/cuffmerge/assembly.txt &> results/cuffmerge/merged.gtf.log

# Grab just novel transcripts that match class_code u, x, o, s, or j.
echo 'class_code "u";' >  novel_transcript_codes.txt
echo 'class_code "x";' >> novel_transcript_codes.txt
echo 'class_code "o";' >> novel_transcript_codes.txt
echo 'class_code "s";' >> novel_transcript_codes.txt
echo 'class_code "j";' >> novel_transcript_codes.txt

grep -f novel_transcript_codes.txt results/cuffmerge/merged.gtf \
    > results/cuffmerge/novel_transcript.gtf

# ----------------------------------------------------------------------------------------
# --- Add gene name to novel transcripts
# ----------------------------------------------------------------------------------------

NOVEL=results/cuffmerge/novel_transcript.gtf

sed '/gene_name/! s/.*gene_id \("[^ ]*"\).*/\0 gene_name \1;/g' \
    $NOVEL > ${NOVEL/.gtf/_gn.gtf}

# ----------------------------------------------------------------------------------------
# --- Combine novel and known transcripts
# ----------------------------------------------------------------------------------------

mkdir -p results/new_annotation

NOVEL_GN=results/cuffmerge/novel_transcript_gn.gtf
GTF="genomes/Tsolium_Mexico_v1/taenia_solium.PRJNA170813.WBPS8.canonical_geneset.gtf"

NOVEL_AND_KNOWN=results/new_annotation/all_transcripts.gtf

cat $NOVEL_GN $GTF > $NOVEL_AND_KNOWN

# ----------------------------------------------------------------------------------------
# --- Tally up genes
# ----------------------------------------------------------------------------------------

export PATH=$PATH:$HOME/bin/subread-1.5.2-Linux-x86_64/bin/

mkdir -p results/counts

NOVEL_AND_KNOWN=results/new_annotation/all_transcripts.gtf

BAMS=`ls results/bam/Taenia_PR*.bam | grep -v "min" | grep -v "plus"`

featureCounts -p -s 2 -T 15 -t exon -g gene_id \
    -a $NOVEL_AND_KNOWN \
    -o results/counts/gene_counts.txt \
    ${BAMS[*]} &> results/counts/gene_counts.txt.log

cut -f 1,7- results/counts/gene_counts.txt | \
    awk 'NR > 1' | \
    awk '{ gsub("samples/bam/","",$0); print }' > results/counts/gene_counts_mini.txt
