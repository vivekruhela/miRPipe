#! /bin/sh

# Adapter Trimming

date
HOME_DIR=`pwd`
SEQ_DIR=$HOME_DIR/data
mkdir -p $SEQ_DIR/trimmed_fastq

echo "********Starting Adapter Trimming***********"

cd $HOME_DIR/Tools/TrimGalore-0.6.1
for fq in $SEQ_DIR/*.fastq.gz $SEQ_DIR/*.fq.gz; do
	if [ ${file: -9} == ".fastq.gz" ]
	then
		basename=$(basename "$fq" .fastq.gz | cut -f1 -d '_')
	else
		basename=$(basename "$fq" .fq.gz | cut -f1 -d '_')
	fi
	
	if [ -z "$adaptor" ]; 
	then
		./trim_galore $fq --length 17 -o $SEQ_DIR/trimmed_fastq --gzip --fastqc
	else
		./trim_galore --length 17 -a $adaptor $fq -o $SEQ_DIR/trimmed_fastq --gzip --fastqc
	fi
    echo "-----------Adapter Trimming of $fq is complete-----------"
done
echo "********Adapter Trimming is complete***********"

