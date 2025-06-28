# Download dorado basecalling software.
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.3-linux-x64.tar.gz

# Uncompress the .tar.gz folder.
tar -xzvf dorado-0.5.3-linux-x64.tar.gz
# The dorado program is located within the bin folder inside the uncompressed folder.
# When using dorado on Roar Collab, I find that I have to specify the path to where it is located, if not the system won't be able to find it.

# Optional: Download basecalling model. See https://github.com/nanoporetech/dorado/blob/release-v0.5.3/README.md#dna-models for a list of models to choose from.
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.3.0

# Below is the workflow for basecalling individual pod5 files that have not already been demultiplexed in the sequencer.

    ## 1. run your dataset through simplex basecalling with barcoding enabled: dorado basecaller <model> <pod5> --kit-name <barcode-kit> > calls.bam 
    ## 2. demultiplex the basecalled bam file: dorado demux --no-classify --output-dir classify and split the dataset
    ## 3. then fetch the read ids per barcode from the corresponding .bam and put it in a read.txt file, eg. for barcode01: samtools view SQK-NBD114-24_barcode01.bam | cut -f1 | sort | uniq > BC01_read_names.txt
    ## 3. run dorado duplex <model> <pod5> --read-ids BC01_read_names.txt and this will run duplex basecalling only with the read ids from that barcode

###################################################################################################################
# STEP 1: SIMPLEX BASECALLING & BARCODING #
###################################################################################################################

#!/bin/bash
#SBATCH --job-name=all_0
#SBATCH --output=all_0.out
#SBATCH --time=24:00:00
#SBATCH --account=mum55_1gc5gb 
#SBATCH --partition=sla-prio
#SBATCH --gpus=1

# Print starting date and time.
date

# Change directory to the folder that contains all the pod5 files for one barcode.
cd /storage/home/yvl6147/scratch/Symb_microbiomes_all_run4-pilot_try2/Symb_microbiomes_all_run4-pilot/run4-pilot/20240601_1557_MC-114792_FAZ36568_87d358a2/pod5

# Create a folder for the output fastq files to go into.
mkdir /storage/home/yvl6147/scratch/Symb_microbiomes_all_run4-pilot_try2/dorado_seqs_simplex/

# Define variables.
export MODEL=/storage/home/yvl6147/scratch/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 # define super accurate basecalling model for dorado
# Alternatively, MODEL=path_to_downloaded_basecalling_model
export OUT=/storage/home/yvl6147/scratch/Symb_microbiomes_all_run4-pilot_try2/dorado_seqs_simplex/FAZ36568_87d358a2_7a4fa985_0.bam # specify output directory
export IN=/storage/home/yvl6147/scratch/Symb_microbiomes_all_run4-pilot_try2/Symb_microbiomes_all_run4-pilot/run4-pilot/20240601_1557_MC-114792_FAZ36568_87d358a2/pod5/FAZ36568_87d358a2_7a4fa985_0.pod5 # spe$

# Run basecalling with barcoding for the single pod5 file.
echo $IN # sanity check to make sure the correct pod5 file is being run.
/storage/home/yvl6147/scratch/dorado-0.5.3-linux-x64/bin/dorado basecaller $MODEL $IN --kit-name SQK-NBD114-24 > $OUT
echo $OUT # sanity check to make sure the correct bam file is being outputted.

# Print completion date and time.
date


###################################################################################################################

# We need to run the above script for each pod5 files. So we have to copy the script, edit each script, and run each script. To do automate this I ran the following code:
## Eg. To run this for the firt 100 pod5 files,
for i in {0..1}; do echo $i; cp basecalling_allsimplex.sh basecalling_allsimplex_${i}.sh; sed -i "s/_0/_${i}/g" basecalling_allsimplex_${i}.sh; bash basecalling_allsimplex_${i}.sh; done

for i in {1121..1433}; do echo $i; cp basecalling_allsimplex.sh basecalling_allsimplex_${i}.sh; sed -i "s/_0/_${i}/g" basecalling_allsimplex_${i}.sh; bash basecalling_allsimplex_${i}.sh; done

    #Line 1: Checks that for loop is iterating 0 - 99 times, this should coincide with the number of pod5files you have. For barcode02 I had 100 files.
    #Line 2: This loop will copy the script you made with a new name.
    #Line 3: This loop will go into each script and change all instances of "_0" to "_${i}" so each script runs its coinciding pod5 file.
    #Line 4: submits each job


###################################################################################################################
# STEP 2 & 3: DEMULTIPLEXING ALL BARCODES #
###################################################################################################################

## Merge all basecalled BAM files into single BAM file.
samtools merge finalBamFile.bam *.bam

## Make directory to store demuxed BAM files.
mkdir demux

## Demultiplex
/storage/home/yvl6147/scratch/dorado-0.5.3-linux-x64/bin/dorado demux --output-dir demux --no-classify finalBamFile.bam

## Inside the demuxed folder, fetch the read ids per barcode from the corresponding .bam and put it in a read.txt file.
for i in *.bam;
do
  samtools view $i | cut -f1 | sort | uniq > ${i%.bam}'_read_names.txt';
done


###################################################################################################################
# STEP 4: RUN DUPLEX BASECALLING ON READ IDS PER BARCODE #
###################################################################################################################

# Script to run for basecalling all BC01 sequences:

# Print starting date and time.
date

# Change directory to the folder that contains all the pod5 files before demux.
cd /storage/home/yvl6147/scratch/Symb_microbiomes_all_run2/run2/20240527_1458_MC-114792_FAZ02813_1f0d390f/pod5

# Create a folder for the output fastq files to go into.
mkdir /storage/home/yvl6147/scratch/Symb_microbiomes_all_run2/dorado_fastq_demuxed

# Create a folder for BC01 fastq files to go into.
mkdir /storage/home/yvl6147/scratch/Symb_microbiomes_all_run2/dorado_fastq_demuxed/barcode01

# Define variables.
export MODEL=/storage/home/yvl6147/scratch/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 # define super accurate basecalling model for dorado
# Alternatively, MODEL=path_to_downloaded_basecalling_model
export OUT=/storage/home/yvl6147/scratch/Symb_microbiomes_all_run2/dorado_fastq_demuxed/barcode01/FAZ02813_1f0d390f_4ce703e5_0.fastq # specify output directory
export IN=/storage/home/yvl6147/scratch/Symb_microbiomes_all_run2/run2/20240527_1458_MC-114792_FAZ02813_1f0d390f/pod5/FAZ02813_1f0d390f_4ce703e5_0.pod5 # specify input directory
export READLIST=/storage/home/yvl6147/scratch/Symb_microbiomes_all_run2/dorado_seqs_simplex/demux/SQK-NBD114-24_barcode01_read_names.txt # Specify where to find the list of all read names with barcode01.

# Run duplex basecalling for only the sequences with the same barcode.
echo $IN # sanity check to make sure the correct pod5 file is being run.
/storage/home/yvl6147/scratch/dorado-0.5.3-linux-x64/bin/dorado duplex --emit-fastq $MODEL $IN --read-ids $READLIST > $OUT 
echo $OUT # sanity check to make sure the correct bam file is being outputted.

# Print completion date and time.
date

########## End of script ##########

# Run the script for all pod5 files. This is only for reads with barcode01.
for i in {2..1155}; do echo $i; cp duplex_basecalling_BC01.sh duplex_basecalling_BC01_${i}.sh; sed -i "s/_0/_${i}/g" duplex_basecalling_BC01_${i}.sh; bash duplex_basecalling_BC01_${i}.sh; done

# Run the script for all pod5 files. This is only for reads with barcode02.
for i in {0..1}; do echo $i; cp duplex_basecalling_BC02.sh duplex_basecalling_BC02_${i}.sh; sed -i "s/_0/_${i}/g" duplex_basecalling_BC02_${i}.sh; bash duplex_basecalling_BC02_${i}.sh barcode02; done


mkdir seqs
cd seqs
for i in {1..9}; do mkdir barcode0${i}_04; done
for i in {10..24}; do mkdir barcode${i}_04; done
cd ..
for i in {1..9}; do mv *barcode0${i}* seqs/barcode0${i}_04; done
for i in {10..24}; do mv *barcode${i}* seqs/barcode${i}_04; done


# To start interactive gpu session.
salloc --account=mum55_1gc5gb --partition=sla-prio --gpus=1 -t 24:00:00