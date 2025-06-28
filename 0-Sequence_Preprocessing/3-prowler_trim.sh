# In Roar Collab, the python module must first be loaded.
module load python/3.6.8

# This runs pretty fast and the process is easily monitored by outputs in the terminal, so it seems to be easier to run and keep track of as an interactive script.

########## prowler-int.sh ##########

# Change directory to the folder that contains all dorado basecalled seqeunces for a specific barcode.
cd /storage/home/yvl6147/scratch/Symb_microbiomes_all_run4-pilot_try2/cutadapt_out/

# Make a directory for the trimmed seqs to go into.
#mkdir /storage/home/yvl6147/scratch/Symb_microbiomes_all_run4-pilot_try2/prowler
mkdir /storage/home/yvl6147/scratch/Symb_microbiomes_all_run4-pilot_try2/prowler/prowler_out

# Define variables.
export PROWLER=/storage/home/yvl6147/scratch/MinION_manuscript_field/ProwlerTrimmer/TrimmerLarge.py
export INDIR=/storage/home/yvl6147/scratch/MinION_manuscript_field/cutadapt_out/
export OUTDIR=/storage/home/yvl6147/scratch/MinION_manuscript_field/prowler_out/

# Unzip and trim all fastq files in the directory.
for i in barcode10/*.fastq;
do
  	python3 $PROWLER -f ${i} -i $INDIR -o $OUTDIR -w 60 -l 200 -c "LT" -g "U0" -m "S" -q 20 
done

# Parameters:
## -f, 	--file,		filename:	The name of the file you want to trim, wihtout the folderpath
## -i, 	--infolder, 	inFolder:	The folderpath where your file to be trimmed is located (default = cwd)
## -o, 	--outfolder,	outFolder:	The folderpath where your want to save the trimmed file (default = cwd)
## -w, 	--windowsize,	windowSize:	Change the size of the trimming window (default= 100bp)
## -l, 	--minlen,	minLen:		Change the minimum acceptable numer of bases in a read (default=100)
## -m, 	--trimmode,	mode:		Select trimming algorithm: S for static  or D for dynamic (default=S)
## -q, 	--qscore,	Qcutoff:	Select the phred quality score trimming threshold (default=7)
## -d, 	--datamax,	maxDataMB:	Select a maximum data subsample in MB (default=0, entire file)
## -r, 	--outformat,	outMode:	Select output format of trimmed file (fastq or fasta) (default=.fastq)
## -c, 	--clip,		clipping:	Select L to clip leading Ns, T to trim trialing Ns and LT to trim both (default=LT)
## -g, 	--fragments,	fragments:	Select fragmentation mode (default=U0)

########## END ##########

# After trimming, before running Qiime2, we have to concatenate all the separate fastq files into one big fastq file. This will make it easier to run the Qiime2 pipeline.
# Since each barcode contains multiple fastq files, we must first merge them into one single fastq file per barcode.
cat barcode01/*fastq.gz > barcode01-alltrimmed.fastq.gz # eg. for barcode 01

# We can write a loop to do this iteratively for each barcode folder.
for i in barcode*/;
do
        cat $i/*fastq.gz > ${i%/}-alltrimmed.fastq.gz;
done

# Then create a manifest file for the new concatenated trimmed fastq sequences. 
# Then run Qiime2 to visualize per-base read qualities and distribution of read lengths.


# Run fastqc on all fastq files for quality check purposes.
module load fastqc
fastqc *.fastq -t 4
# Uses 4 threads to run this command. This makes it run faster, but using > 4 threads may cause result in the process getting killed prematurely, idk why.

# Combine all fastqc reports into one multiqc report.
multiqc .
# In Roar Collab, multiqc is downloaded in the qiime2 conda environment.