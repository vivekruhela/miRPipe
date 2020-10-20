#! /bin/sh

# Alignment using mirdeep star

date
HOME_DIR=`pwd`
SEQ_DIR=$HOME_DIR/data
Tools_DIR=$Tools_DIR

if [ -d "$Tools_DIR/MDS_command_line_v37/MDS_command_line/genome/hg19" ]; then 

	REF_DIR=$Tools_DIR/MDS_command_line_v37/MDS_command_line
	MD_INDEX=hg19

fi

if [ -d "$Tools_DIR/MDS_command_line_v38/MDS_command_line/genome/hg19" ]; then 

	REF_DIR=$Tools_DIR/MDS_command_line_v38/MDS_command_line
	MD_INDEX=hg19

fi

if [ -d "$Tools_DIR/MDS_command_line_v38/MDS_command_line/genome/hg38" ]; then 

	REF_DIR=$Tools_DIR/MDS_command_line_v38/MDS_command_line	
	MD_INDEX=hg38

fi

if [ -z "$adaptor" ]; 
	then
		echo "Default Adaptor sequene is choosen. Then using default adaptor sequence $adaptor"
	else
		echo "Adaptor sequence $adaptor will be used for adaptor trimming."
fi

cd $REF_DIR

for file in $SEQ_DIR/* ; do		
	if ! [ -d "$file" ]; then
		if ! [[ $file =~ \.gz$ ]]; # check file extension if fastq.gz or fq.gz 
		then
			echo "Aligning $file with human reference genome...."
			# Please add adaptor according to your data type.		
			java -jar -Xmx8g MD.jar -g $MD_INDEX -a $adaptor -t 18 -l 23 -p 20 -m 101 -r 5 -s -10 $file
			echo "$file alignment is complete"
		else
			echo "Input $file is in compressed format. Uncompressing the $file...."
			if [[ $file =~ \.fq.gz$ ]];
			then
				basename=$(basename "$file" .fq.gz | cut -f1 -d '_')			
			else
				if [[ $file =~ \.fastq.gz$ ]];
					then
						basename=$(basename "$file" .fastq.gz | cut -f1 -d '_')					
				fi
			fi
			
			pigz -p 5 -f -k -d $file -c > $SEQ_DIR/$basename.fastq
			echo "Uncompression is complete. Now Aligning $file with human reference genome...."
			# Please add adaptor according to your data type.
			java -jar -Xmx8g MD.jar -g $MD_INDEX -a $adaptor -t 18 -l 23 -p 20 -m 101 -r 5 -s -10 $SEQ_DIR/$basename.fastq
			echo "$file alignment is complete"
			rm $SEQ_DIR/$basename.fastq

		fi
	fi		
done

