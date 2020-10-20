#! /bin/sh

# exit when any command fails
set -e

HOME_DIR=`pwd`
SEQ_DIR=$HOME_DIR/data
REF_DIR=$HOME_DIR/refs

if [ -d "$REF_DIR/hg19" ]; then 
     echo "Selecting hg19 as reference genome."
     BWT_INDEX=$REF_DIR/hg19/hg19; 
     wget --content-disposition http://www.pirnadb.org/download/downloadarchive/gff_gtf/pirnadb.v1_7_5.hg19.gff3.gz -P $REF_DIR/pirnadb
     pigz -p 4 -k -d $REF_DIR/pirnadb/pirnadb.v1_7_5.hg19.gff3.gz
     GFF_FILE=$REF_DIR/pirnadb/pirnadb.v1_7_5.hg19.gff3

else
    if [ -d "$REF_DIR/hg38" ]; then 
         echo "Selecting hg38 as reference genome."
         BWT_INDEX=$REF_DIR/hg38/hg38; 
         pigz -p 4 -k -d $REF_DIR/pirnadb/pirnadb.hg38.gff3.gz
         GFF_FILE=$HOME_DIR/refs/pirnadb/pirnadb.hg38.gff3
         
     else
         trap 'echo "Please provide piRNA coordinate file (.gff3 or .gff3.gz)"' EXIT
         

    fi

fi

# piRNA pipeline

mkdir -p $SEQ_DIR/piRNA
mkdir -p $SEQ_DIR/piRNA/pirna_counts

for fq in $SEQ_DIR/fastq_24_31/*.fq.gz; do 
     basename=$(basename "$fq" .fq.gz)
     echo "Working on $basename"
     pigz -p 4 -k -d $fq
     # Trimming reads having more than 10% bp with quality score less than 20
     cd $HOME_DIR/Tools/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64
     ./fastq_quality_filter -q 20 -p 90 -Q33 -z \
                            -i $SEQ_DIR/fastq_24_31/$basename.fq \
                            -o $SEQ_DIR/piRNA/$basename.fq.gz

     
     echo "****Quality Trimming of $basename is complete"
     
     # Starting Alignment
     echo "**** Getting alignment for $basename ****"
     cd $HOME_DIR/Tools/bowtie-1.2.3-linux-x86_64
     ./bowtie -p 4 $BWT_INDEX \
              -v 0 --norc -l 24 $SEQ_DIR/piRNA/$basename.fq.gz \
              -S $SEQ_DIR/piRNA/$basename"_bowtie".sam 
     samtools view -h -b $SEQ_DIR/piRNA/$basename"_bowtie".sam -o $SEQ_DIR/piRNA/$basename"_bowtie".bam
     samtools sort -@ 8 $SEQ_DIR/piRNA/$basename"_bowtie".bam > $SEQ_DIR/piRNA/$basename"_bowtie_sorted".bam 
     
     echo "****Sequence alignemnt of $basename is complete"
     # Deleting the temporary files
     rm $SEQ_DIR/piRNA/$basename"_bowtie".sam
     rm $SEQ_DIR/piRNA/$basename"_bowtie".bam
     rm $SEQ_DIR/piRNA/$basename.fq.gz

     # Starting read counts
     echo "**** Getting piRNA read counts from $basename ****"
     bedtools bamtobed -i $SEQ_DIR/piRNA/$basename"_bowtie_sorted".bam > $SEQ_DIR/piRNA/pirna_counts/$basename.bed
     bedtools annotate -counts -i $GFF_FILE -files $SEQ_DIR/piRNA/pirna_counts/$basename.bed > $SEQ_DIR/piRNA/pirna_counts/$basename"_pirna_counts".txt
     pigz -p 4 $SEQ_DIR/piRNA/pirna_counts/$basename.bed
     rm $SEQ_DIR/piRNA/$basename"_bowtie_sorted".bam    
     echo "****Mapped sequence quantification and annotation of $basename is complete"

done
