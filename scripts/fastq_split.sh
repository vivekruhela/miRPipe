#! /bin/sh

# length and mean quality filtering
HOME_DIR=`pwd`
SEQ_DIR=$HOME_DIR/data
mkdir -p $SEQ_DIR/fastq_17_24
mkdir -p $SEQ_DIR/fastq_24_31
mkdir -p $SEQ_DIR/fastq_31
mkdir -p $SEQ_DIR/filtered_out
mkdir -p $SEQ_DIR/filtered_out1
cd $HOME_DIR/Tools/bbmap
echo "********Length and Quality Filteration (17nt <= Length <= 24 nt and Quality >= 20)***********"
for fq in $SEQ_DIR/trimmed_fastq/*.fq.gz;do
    basename=$(basename "$fq" .fq.gz)
    echo "Processing $basename"
    ./bbduk.sh -Xmx1g minlength=17 \
                      maxlength=24 \
                      minavgquality=20 \
                      in=$fq \
                      out=$SEQ_DIR/fastq_17_24/$basename.fastq.gz \
                      outm=$SEQ_DIR/filtered_out/$basename.fq.gz
done
pigz -p 4 -k -d $SEQ_DIR/fastq_17_24/*.fastq.gz

echo "********Length and Quality Filteration (24nt <= Length <= 31 nt and Quality >= 20)***********"
for fq in $SEQ_DIR/filtered_out/*.fq.gz;do
    basename=$(basename "$fq" .fq.gz)
    echo "Processing $basename"
    ./bbduk.sh -Xmx1g minlength=24 \
                      maxlength=31 \
                      minavgquality=20 \
                      in=$fq \
                      out=$SEQ_DIR/fastq_24_31/$basename.fq.gz \
                      outm=$SEQ_DIR/filtered_out1/$basename.fq.gz
done

echo "********Length and Quality Filteration (31nt <= Length and Quality >= 20)***********"
for fq in $SEQ_DIR/filtered_out1/*.fq.gz;do
    basename=$(basename "$fq" .fq.gz)
    echo "Processing $basename"
    ./bbduk.sh -Xmx1g minlength=32 \
                      minavgquality=20 \
                      in=$fq \
                      out=$SEQ_DIR/fastq_31/$basename.fq.gz
done

rm -r $SEQ_DIR/filtered_out
rm -r $SEQ_DIR/filtered_out1