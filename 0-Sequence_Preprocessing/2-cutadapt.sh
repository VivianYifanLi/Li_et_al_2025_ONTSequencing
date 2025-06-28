### Convert bam files to fastq. ###

for i in *.bam;
do
    samtools bam2fq $i > ${i%.bam}.fastq;
done

####################################################################

# Combine all fastq per barcode into a single file.
for i in barcode*/
do
	cat $i/*fastq.gz > ${i%/}-allpass.fastq.gz
done


# Make directory for storing the cutadapt output sequences.
mkdir /storage/home/yvl6147/scratch/Symb_microbiomes_all_run1/cutadapt_out


### bash script ###

conda activate qiime2 # The cutadapt tool is downloaded to the qiime2 environment in Roar Collab.

export INDIR=/storage/home/yvl6147/scratch/Symb_microbiomes_all_run4-pilot_try2/dorado_seqs_simplex/run4_simplex_demux
export OUTDIR=/storage/home/yvl6147/scratch/Symb_microbiomes_all_run4-pilot_try2/cutadapt_out

# Run cutadapt
cd $INDIR
for i in *.fastq;
do
    cutadapt \
    -b TTTCTGTTGGTGCTGATATTGCAGRGTTYGATYMTGGCTCAG \
    -b ACTTGCCTGTCGCTCTATCTTCCGGYTACCTTGTTACGACTT \
    -j 0 \
    ${INDIR}/$i > ${OUTDIR}/${i%.fastq}_noprimers.fastq \
    2> $OUTDIR/${i%.fastq}_cutadapt_stats.txt \
    --info-file=${OUTDIR}/${i%.fastq}_cutadapt_info.tsv;
done

####################################################################

### Export QC report of summary statistics. ###

multiqc .