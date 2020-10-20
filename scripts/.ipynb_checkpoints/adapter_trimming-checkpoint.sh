#! /bin/sh

# Adapter Trimming

date
HOME_DIR=`pwd`
SEQ_DIR=$HOME_DIR/data
mkdir -p $SEQ_DIR/trimmed_fastq

echo "********Starting Adapter Trimming***********"

cd $HOME_DIR/Tools/TrimGalore-0.6.0
for fq in $SEQ_DIR/*.fastq.gz; do
	basename=$(basename "$fq" .fastq.gz | cut -f1 -d '_')
	if [ -z "$adaptor" ]; 
		then
			./trim_galore $fq -o $SEQ_DIR/trimmed_fastq --gzip
		else
			./trim_galore -a $adaptor $fq -o $SEQ_DIR/trimmed_fastq --gzip
	fi
    echo "-----------Adapter Trimming of $fq is complete-----------"
done
echo "********Adapter Trimming is complete***********"

