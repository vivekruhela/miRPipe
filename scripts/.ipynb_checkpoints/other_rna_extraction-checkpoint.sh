#! /bin/sh
date

HOME_DIR=`pwd`
SEQ_DIR=$HOME_DIR/data
REF_DIR=$HOME_DIR/refs

if [ -d "$REF_DIR/hg19" ]; then BWT_INDEX=$REF_DIR/hg19/hg19; fi

if [ -d "$REF_DIR/hg38" ]; then BWT_INDEX=$REF_DIR/hg38/hg38; fi

# piand sno-RNA pipeline

mkdir -p $SEQ_DIR/piRNA
mkdir -p $SEQ_DIR/piRNA/pirna_counts
mkdir -p $SEQ_DIR/snoRNA
mkdir -p $SEQ_DIR/snoRNA/snorna_counts

GFF_FILE="$HOME_DIR/refs/pirnadb/pirnadb.hg19.gff3"

if ! [ -f "$GFF_FILE" ]; then
    pigz -p 5 -d $HOME_DIR/refs/pirnadb/pirnadb.hg19.gff3.gz
fi

for fq in $SEQ_DIR/trimmed_fastq/*.fq.gz; do 
     basename=$(basename "$fq" .fq.gz)
     
     # Trimmimng reads shorter than 15
     cd $HOME_DIR/Tools/bbmap
     ./bbduk.sh -Xmx1g minlength=15 \
                       in=$fq \
                       out=$SEQ_DIR/piRNA/$basename"_1".fq.gz

     # Trimming reads having more than 10% bp with quality score less than 20
     cd $HOME_DIR/Tools/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64
     pigz -p 5 -d $SEQ_DIR/piRNA/$basename"_1".fq.gz
     ./fastq_quality_filter -q 20 -p 90 -Q33 -z \
                            -i $SEQ_DIR/piRNA/$basename"_1.fq" \
                            -o $SEQ_DIR/piRNA/$basename"_2".fq.gz

     
     # Starting Alignment
     cd $HOME_DIR/Tools/bowtie-1.2.2-linux-x86_64
     echo "**** Getting alignment for $basename ****"
     ./bowtie -p 8 $BWT_INDEX \
              --no-unal $SEQ_DIR/piRNA/$basename"_2".fq.gz \
              -S $SEQ_DIR/piRNA/$basename"_bowtie".sam 
     samtools view -h -b $SEQ_DIR/piRNA/$basename"_bowtie".sam -o $SEQ_DIR/piRNA/$basename"_bowtie".bam
     samtools sort -@ 8 $SEQ_DIR/piRNA/$basename"_bowtie".bam > $SEQ_DIR/piRNA/$basename"_bowtie_sorted".bam 
     
     # Deleting the temporary files
     rm $SEQ_DIR/piRNA/$basename"_bowtie".sam
     rm $SEQ_DIR/piRNA/$basename"_bowtie".bam
     rm $SEQ_DIR/piRNA/$basename"_1.fq"
     rm $SEQ_DIR/piRNA/$basename"_2".fq.gz

     # Starting read counts
     echo "**** Getting piRNA read counts from $basename ****"
     bedtools bamtobed -i $SEQ_DIR/piRNA/$basename"_bowtie_sorted".bam > $SEQ_DIR/piRNA/pirna_counts/$basename.bed
     bedtools annotate -counts -i $HOME_DIR/refs/pirnadb/pirnadb.hg19.gff3 -files $SEQ_DIR/piRNA/pirna_counts/$basename.bed > $SEQ_DIR/piRNA/pirna_counts/$basename"_pirna_counts".txt
     bedtools annotate -counts -i $HOME_DIR/refs/snodb/snoDB.gff3 -files $SEQ_DIR/piRNA/pirna_counts/$basename.bed > $SEQ_DIR/snoRNA/snorna_counts/$basename"_snorna_counts".txt
     pigz -p 5 $SEQ_DIR/piRNA/pirna_counts/$basename.bed
     rm $SEQ_DIR/piRNA/$basename"_bowtie_sorted".bam    

done

date