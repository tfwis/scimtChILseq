#!/usr/local/packages/mambaforge/bin/zsh

## sci-ChIL-seq preprocessing

### Setup

cellbcs=def_cellbc.tsv.gz
# please download from Zenodo

#### mouse

INDEX=bowtie2/mm10/mm10 # Bowtie2 index plese download from https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
usechrs=chr.tsv
genomefile=mm10.genome

wd=`pwd`
FASTQ_DIR=${wd}/fastq # please download FASTQ files from GEO and set under this directory

### 

cd ${wd}

if [[ ! -d `dirname ${INDEX}` ]]; then
  echo No Bowtie2 INDEX file; exit 0 
fi

if [[ ! -f ${usechrs} ]]; then
  echo No Chromosome file; exit 0 
fi

./chilmap.bash ${FASTQ_DIR} ${INDEX} ${usechrs}
Rscript ./chilrmdups.R . cellbcs
Rscript ./chilqcs.R .

rm -rf sep
rm -f mapped.bam
rm -r -f .tmp

#####