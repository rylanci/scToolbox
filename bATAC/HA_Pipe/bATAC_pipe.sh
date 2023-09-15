#!/bin/bash
#PBS -q hotel
#PBS -N bATAC_pipe
#PBS -l nodes=1:ppn=8
#PBS -l walltime=4:00:00
#PBS -M rlancione@ucsd.edu
#PBS -m ae
#PBS -e /home/rlancione/jeo/FFS/bATAC_pipe.e
#PBS -o /home/rlancione/jeo/FFS/bATAC_pipe.o

# activate cond env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bATAC1

#make a project directory and navigate to the directory and create the following directories
#e.g. mkdir fastqc rawdata mapping trimmed bw log macs2
#soft link ATAC fastq files into fastq_dir
#e.g. cp -s fastq_files ./fastq_dir

proj_dir="/projects/ps-epigen/users/rlan/ATACSeq/HA_Pipe/output/cow_control/output/"
fastq_dir="/projects/ps-epigen/users/rlan/ATACSeq/HA_Pipe/output/cow_control/fastqs/"
dsheet="/projects/ps-epigen/users/rlan/ATACSeq/HA_Pipe/output/cow_control/datasheet.csv"


# For testing
#PBS_ARRAYID=2
#Take sample name from sample sheet one row per job
sample=`tail -n+2 $dsheet |sed -n ${PBS_ARRAYID}p|cut -f 1 -d ','` 
genome=`tail -n+2 $dsheet |sed -n ${PBS_ARRAYID}p|cut -f 3 -d ','`
r1=`ls $fastq_dir | grep $sample | grep _R1`
r2=`ls $fastq_dir | grep $sample | grep _R2`

echo "Sample: $sample"
echo "R1P: ${fastq_dir}${r1}"
echo "R2P: ${fastq_dir}${r2}"
echo "Genome: $genome"


fastqc_dir=${proj_dir}/fastqc
mkdir $fastqc_dir
map_dir=${proj_dir}/mapping
mkdir $map_dir
trim_dir=${proj_dir}/trimmed
mkdir $trim_dir
bw_dir=${proj_dir}/bw
mkdir $bw_dir
log_dir=${proj_dir}/log
mkdir $log_dir
peak_dir=${proj_dir}/macs2
mkdir $peak_dir 

#mm10_bowtie2=/projects/ps-renlab/share/bowtie2_indexes/mm10
hg38_bowtie2="/home/rlancione/ps-epigen/users/rlan/ATACSeq/PEPATAC/genomes/GRCh38_noalt_as/GRCh38_noalt_as"
#btar_bowtie2="/projects/ps-epigen/GENOME/cattle/cattle"
#btar_bowtie2="/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/bosTau8mask_indexes/bosTau8"
btar_bowtie2="/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/ARS-UCD1.3_bowtie/ARS-UCD1.3"
#gal_bowtie2=""
#mm10_bed=/projects/ps-renlab/y2xie/projects/genome_ref/gencode.vM25.annotation.gtf
hg38_bed="/projects/ps-epigen/GENOME/hg38/gencode.v33.annotation.gtf"
#btar_bed="/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/btaurus_gtf/Bos_taurus.ARS-UCD1.2.108.gtf"
#btar_bed="/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/btaur_genomic_info/ncbi_dataset/data/GCF_002263795.2/genomic.gtf"
btar_bed="/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/ARS-UCD1.3.gtf"
#gal_bed=""
### Rat
rn7_bowtie2="/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/Rat_Samples/rn7_index/rn7"
rn7_bed="/projects/ps-epigen/users/rlan/ATACSeq/zemke_scripts/Rat_Samples/ncbiRefSeq.gtf"

#if [ -z ${genome+x} ]; then echo "genome is not specified. Default to mm10 (bc I want)"; genome="mm10"; fi
if [ $genome == "mm10" ]; then ref=${mm10_bowtie2}; macs2_genome="mm"; TSS_bed=${mm10_bed}; fi
if [ $genome == "human" ]; then ref=${hg38_bowtie2}; macs2_genome="hs"; TSS_bed=${hg38_bed}; fi
if [ $genome == "cow" ]; then ref=${btar_bowtie2}; macs2_genome="tr"; TSS_bed=${btar_bed}; fi
#if [ $genome == "chicken" ]; then ref=${Gal_bowtie2}; macs2_genome="hs"; TSS_bed=${hg38_bed}; fi
if [ $genome == "rat" ]; then ref=${rn7_bowtie2}; macs2_genome="tr"; TSS_bed=${rn7_bed}; fi 

#echo before comment
#: <<'END'

fastqc -t 16 -o ${fastqc_dir} ${fastq_dir}/${r1}
fastqc -t 16 -o ${fastqc_dir} ${fastq_dir}/${r2}
trim_galore -q 20 --paired ${fastq_dir}/${r1} ${fastq_dir}/${r2} -o ${trim_dir} --basename ${sample}

# mapping
# 2> redirects the bowtie stder to the log
#r1t=$sample"_R1_001_val_1.fq.gz"
#r2t=$sample"_R2_001_val_2.fq.gz"
r1t=$sample"_val_1.fq.gz"
r2t=$sample"_val_2.fq.gz"

#ls Murdoch_QY_1755-1766/trimmed/ | grep -e QY_1755_S1_R1.*gz$ 
bowtie2 -x ${ref} -1 ${trim_dir}/${r1t} -2 ${trim_dir}/${r2t} --no-unal -p 16 -S ${map_dir}/${sample}_${genome}.sam 2> ${map_dir}/${sample}.log
samtools sort -@ 16 -T ${map_dir} -o ${map_dir}/${sample}_${genome}_sorted.bam ${map_dir}/${sample}_${genome}.sam
#if [[ -f "${map_dir}/${sample}_${genome}_sorted.bam" ]]; then rm ${map_dir}/${sample}_${genome}.sam; fi

# deduplicate
java -Xmx8G -XX:ParallelGCThreads=3 -jar /home/rlancione/miniconda3/envs/snATAC/share/picard-2.25.2-0/picard.jar \
	MarkDuplicates I=${map_dir}/${sample}_${genome}_sorted.bam TMP_DIR=${map_dir} \
	METRICS_FILE=${log_dir}/${sample}_${genome}_dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=TRUE \
	O=${map_dir}/${sample}_${genome}_sorted_rmdup.bam REMOVE_DUPLICATES=true
#END
#echo "after comment"

# genreate bw, TSS enrichment
samtools index ${map_dir}/${sample}_${genome}_sorted_rmdup.bam
#bamCoverage -b ${map_dir}/${sample}_${genome}_sorted_rmdup.bam -o ${bw_dir}/${sample}_${genome}_sorted_rmdup.bw -p max --normalizeUsing RPKM

#computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -R ${TSS_bed} \
#	-S ${bw_dir}/${sample}_${genome}_sorted_rmdup.bw --skipZeros \
#	-o ${bw_dir}/${sample}_${genome}_sorted.matrix.gz -p 16

#plotProfile -m ${bw_dir}/${sample}_${genome}_sorted.matrix.gz -o ${bw_dir}/${sample}_${genome}_sorted.profile.png --refPointLabel "TSS"

### give a one-number TSS enrichment calcualetd by ATACseqQC
conda deactivate 
conda activate renv4
### argss are bam, genome, sample name
Rscript /projects/ps-epigen/users/rlan/ATACSeq/HA_Pipe/scripts/bATAC_tss.R ${map_dir}/${sample}_${genome}_sorted_rmdup.bam ${genome} $sample

# call peaks with macs2
#/home/y2xie/anaconda3/envs/macs2/bin/macs2 callpeak -t ${map_dir}/${s}_${genome}_sorted_rmdup.bam --outdir ${peak_dir} -n ${s}_${genome} -q 0.05 --nomodel --keep-dup all --shift -100 --extsize 200 -g ${macs2_genome}


# Remove mito reads & calculate %mito
conda deactivate 
conda activate bATAC1

# Calculate/aquire qc metrics and add to file
samtools idxstats ${map_dir}/${sample}_${genome}_sorted_rmdup.bam | cut -f 1 | \
	grep -v chrM | xargs samtools view -b ${map_dir}/${sample}_${genome}_sorted_rmdup.bam \
	> ${map_dir}/${sample}_${genome}_sorted_rmdup_nomt.bam
### mito percentage
samtools index ${map_dir}/${sample}_${genome}_sorted_rmdup_nomt.bam
bvalue=`samtools view ${map_dir}/${sample}_${genome}_sorted_rmdup.bam | wc -l`
avalue=`samtools view ${map_dir}/${sample}_${genome}_sorted_rmdup_nomt.bam | wc -l`
mito_percent=`echo "(1 - $avalue / $bvalue) * 100" | bc -l` 
### Dup rate
duprate=`head -8 ${log_dir}/${sample}_${genome}_dup.qc | tail -1 | cut -f 9`
### Mapping rate 
maprate=`tail -1 ${map_dir}/${sample}.log`

# append to QY_1763_S9_cow_sorted_rmdup.bam.tsse.score.txt
echo Mito Percent: $mito_percent >> ${map_dir}/${sample}_${genome}_sorted_rmdup.bam.tsse.score.txt
echo Duplication Rate: $duprate >> ${map_dir}/${sample}_${genome}_sorted_rmdup.bam.tsse.score.txt
echo Mapping Rate: $maprate >> ${map_dir}/${sample}_${genome}_sorted_rmdup.bam.tsse.score.txt
echo >> ${map_dir}/${sample}_${genome}_sorted_rmdup.bam.tsse.score.txt




