#!/usr/local/packages/mambaforge/bin/zsh

## Extracting Cell Barcode and Mapping Reads


FASTQ_DIR=./fastq
INDEX=bowtie2/mm10/mm10
usechrs=chr.tsv

wd=`pwd`

I1=`ls ${FASTQ_DIR}/* | grep I1
I2=`ls ${FASTQ_DIR}/* | grep I2
R1=`ls ${FASTQ_DIR}/* | grep R1
R2=`ls ${FASTQ_DIR}/* | grep R2
[[ ! -d log ]] && mkdir log
[[ ! -d .tmp ]] && mkdir .tmp
[[ ! -d sep ]] && mkdir sep

seqkit fx2tab ${I1} > .tmp/i7.tsv & 
        seqkit subseq -r 1:6 ${I2}  | seqkit fx2tab | awk '{print $3,$4}' > .tmp/t5_seq.txt &
        seqkit subseq -r -8:-1 ${I2}  | seqkit fx2tab | awk '{print $3,$4}' > .tmp/i5_seq.txt
paste .tmp/i7.tsv .tmp/i5_seq.txt .tmp/t5_seq.txt |\
        awk '{print $1,$2"\t"$3$5$7"\t"$4$6$8}' |\
        seqkit tab2fx > bc.fastq
umi_tools extract --log2stderr -I bc.fastq -S .tmp/bc_out_R1.fastq.gz \
        --bc-pattern=NNNNNNNNNNNNNNNNCCCCCC --read2-in=${R1} \
        --read2-out=.tmp/R1_bc.fastq.gz 2>log/umi_tools_R1.log.txt &
umi_tools extract --log2stderr -I bc.fastq -S .tmp/bc_out_R2.fastq.gz \
        --bc-pattern=NNNNNNNNNNNNNNNNCCCCCC --read2-in=${R2} \
        --read2-out=.tmp/R2_bc.fastq.gz 2>log/umi_tools_R2.log.txt

trim_galore --paired --clip_R2 5 -j 6 .tmp/R1_bc.fastq.gz .tmp/R2_bc.fastq.gz 2>log/trim.log.txt
mv R1_bc.fastq.gz_trimming_report.txt log
mv R2_bc.fastq.gz_trimming_report.txt log

gzip bc.fastq &
bowtie2 -p16 --reorder -x ${INDEX} \
        -1 R1_bc_val_1.fq.gz -2 R2_bc_val_2.fq.gz 2>log/map.log.txt |\
        samtools sort -m 8G > mapped.bam
samtools view -@12 -hb -f 2 mapped.bam  > pairedmapped.bam 
samtools index -@12 pairedmapped.bam

cat ${usechrs} |\
        while read line; do 
                samtools view -@12 -m 4G -b pairedmapped.bam ${line} > sep/${line}.bam
	done 


###