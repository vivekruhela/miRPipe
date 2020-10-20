#! /bin/sh


# length and mean quality filtering
HOME_DIR=`pwd`
SEQ_DIR=$HOME_DIR/data
mkdir -p $SEQ_DIR/moRNA
mkdir -p $SEQ_DIR/moRNA/filtered_fastq2
cd $HOME_DIR/Tools/bbmap
echo "********Length and Quality Filteration***********"
for fq in $SEQ_DIR/moRNA/filtered_fastq1/*.fq.gz;do
    basename=$(basename "$fq" .fq.gz)
    echo "Processing $basename"
    ./bbduk.sh -Xmx1g maxlength=30 \
                      minavgquality=30 \
                      minlength=18 \
                      in=$fq \
                      out=$SEQ_DIR/moRNA/filtered_fastq2/$basename.fq.gz
done

cd $HOME_DIR

# Remove reads having more than 2 bases less then quality 20
mkdir -p $SEQ_DIR/moRNA/filtered_fastq3
echo "********Filtering sequence having at most 2 base less than phred quality 20***********"
Rscript $HOME_DIR/scripts/qual_trim.R
