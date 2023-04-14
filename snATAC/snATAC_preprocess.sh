#!/bin/sh
#PBS -q home-epigen
#PBS -m ae
#PBS -M rlancione@ucsd.edu
#PBS -V
#PBS -A epigen-group
#PBS -e /home/rlancione/jeo/err/eye_proj/r2-snATAC_array_preproccessingV2.py.e
#PBS -o /home/rlancione/jeo/out/eye_proj/r2-snATAC_array_preproccessingV2.py.o
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=20

# Navigate to oasis scratch
cd /oasis/tscc/scratch/rlancione/
# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snATAC


#### Array setup. ***Edit***: SAMPLE_LIST
# Path to sample sheet
SAMPLE_LIST='/home/rlancione/projects/eye_proj/eye_scic2/''sample_array.csv'
# Parse array data
NAME=`tail -n+2 $SAMPLE_LIST|sed -n ${PBS_ARRAYID}p|cut -f 1 -d ','`
R1=`tail -n+2 $SAMPLE_LIST|sed -n ${PBS_ARRAYID}p|cut -f 2 -d ','`
R2=`tail -n+2 $SAMPLE_LIST|sed -n ${PBS_ARRAYID}p|cut -f 3 -d ','`
#R1P=`ls /projects/ps-epigen/sciATAC_output/CTRP5_MFRP_Eye/LL_117/ | grep R1.fastq`

### Project specific Vars. ***Edit***: Read paths, species, and desired output locations
# R1 Path
R1P='/projects/ps-epigen/sciATAC_output/CTRP5_MFRP_Eye/'$NAME"/"$R1
# R2 Path
R2P='/projects/ps-epigen/sciATAC_output/CTRP5_MFRP_Eye/'$NAME"/"$R2
# Species: options = "mouse" or "human
SPECIES="mouse"
# Oasis output path
OUT='/oasis/tscc/scratch/rlancione/eye_project/eye_sciATACr2/'
# Final output path
FOUT='/projects/ps-epigen/users/rlan/eye_data/eye_sciATACr2/'
####### only "need" to edit above this line #######


# Set vars conditional by species. No need to edit.
if [ "$SPECIES" == "mouse" ]
	then
	## Vars for mouse
	# Reference
	RD='/projects/ps-epigen/GENOME/mm10/bwa_index/mm10_no_alt_analysis_set_ENCODE.fasta'
	# Chrom sizes
	CS='/projects/ps-epigen/GENOME/mm10/mm10.chrom.sizes'
	# Blacklist
	BL='/projects/ps-epigen/GENOME/mm10/mm10.blacklist.bed'
	# Promoter Regions
	PR='/projects/ps-epigen/GENOME/mm10/gencode.vM17.protein_coding.tr.prom.1kb.bed'

elif [ "$SPECIES" == "human" ]
	then
	## Vars for human
	# Reference
	RD='/projects/ps-epigen/GENOME/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta'
	# Chrom sizes
	CS='/projects/ps-epigen/GENOME/hg38/hg38.chrom.sizes'
	# Blacklist
	BL='/projects/ps-epigen/GENOME/hg38/hg38.blacklist.bed'
	# Promoter Regions
	PR='/projects/ps-epigen/GENOME/hg38/gencode.v33.annotation.geneUp2k.bed'
fi


# Picard.jar
PD='/home/rlancione/miniconda3/envs/snATAC/share/picard-2.25.2-0/picard.jar'

# Sanity check
echo "****INPUT****"
echo NAME: $NAME
echo Read1 Path: $R1P
echo Read2 Path: $R2P
echo Chrom sizes: $CS


#### Further define/create OUT locations
# oasis atac out
OOUT="$OUT"atac_prep/
OAOUT="$OOUT"$NAME"/"
# oasis signac out
OSOUT="$OUT"signac_out/
# Final output atac
FAOUT="$FOUT"atac_prep/
# Final output signac
FSOUT="$FOUT"signac_out/

# Instead create downstream dirs out out and fout. /atac_prep, /atac_prep/$NAME, /siganc out
### Create OUTDIR & FOUTDIR
[ ! -d "$OOUT" ] &&  mkdir $OOUT
[ ! -d "$OAOUT" ] &&  mkdir $OAOUT
[ ! -d "$OSOUT" ] &&  mkdir $OSOUT
[ ! -d "$FAOUT" ] &&  mkdir $FAOUT
[ ! -d "$FSOUT" ] &&  mkdir $FSOUT


################ Run snATAC_pipe #################
python3 ~/scripts/ATAC/snATAC_pipelineV2.py -r1 $R1P -r2 $R2P -o $OAOUT -n $NAME -t 18    \
                -ref $RD --picard $PD --chrom-sizes $CS --blacklist $BL --promoter-file $PR
# Use --minimum-reads 0 to see full read distribution. Caveat is much longer runtime

##### Perform indexing on rmdup bam
samtools index $OAOUT"/"$NAME".filt.rmdup.bam"


####  Sinto Setup
# Deactivate snATAC and activate sinto_pipe env
conda deactivate
conda activate sinto_pipe

# Nav to atac OUT for sinto
cd $OAOUT
# Format input
INPUT="$NAME".filt.rmdup.bam

################# Run Sinto Pipe #################
#echo  Create frag file with sinto
sinto fragments -b $INPUT -p 20 -t BX -f $NAME'_fragments.bed'
#echo Begin sorting
sort -k 1,1 -k2,2n $NAME'_fragments.bed' > $NAME'_fragments_sorted.bed'
#echo modify with awk
cat $NAME'_fragments_sorted.bed'  | awk -v NAME=$NAME '{print $1"\t"$2"\t"$3"\t"NAME"_"$4"\t"$5}' > $NAME'_fragments_sorted_barcode_modified.bed'
#echo zip awk output
bgzip $NAME'_fragments_sorted_barcode_modified.bed'
#echo run tabix on zipped bed
tabix -p bed $NAME'_fragments_sorted_barcode_modified.bed.gz'

# Copy atac preprocessing & sinto results to Final atac OUTDIR
cp -r $OAOUT $FAOUT

## catch these errors
#[bgzip] can't create MM_258_fragments_sorted_barcode_modified.bed.gz: File exists
#[tabix] the index file exists. Please use '-f' to overwrite.

#### Perform signac pre-processing
conda deactivate
conda activate signac

##################### Run signac preprocessing #####################
# args = 1.experiment name 2.input dir 3.outdir 4. min_tss 5.min_reads
Rscript ~/scripts/ATAC/signac_pre.R $NAME $FAOUT $OSOUT #5 1000
#****** add args for min tss and uniq usable reads

# Copy signac output.
cp -r $OSOUT/* $FSOUT


echo Done!
