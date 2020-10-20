
##### Running FastQC on the files in file_input #####


SEQ_DIR <- file.path(getwd(),'data')
# dir.create('FastQC_reports') 
if(!file.exists('FastQC_reports')) {  # If the folder does not exist, create a new one                                    
                                    dir.create(file.path(getwd(),'data','FastQC_reports'))
                                   } else 
                                   {   # If it existed, delete and replace with a new one  
                                     unlink(file.path(getwd(),'data','FastQC_reports'), recursive = TRUE)
                                     dir.create(file.path(getwd(),'data','FastQC_reports'))
                                   }
output_path <- file.path(getwd(),'data','FastQC_reports')
fastqc_output_path <- paste('--outdir=',output_path,sep="")


for (file in list.files(SEQ_DIR, pattern = glob2rx("*f*q*"))){
    file = file.path(SEQ_DIR,file)
    command_fastqc <- paste('fastqc', file, fastqc_output_path, sep=" ")
    system(command_fastqc)    
}

# RUN Multiqc to compile all fastqc reports
command_multiqc <- paste("multiqc",output_path,"-o",output_path, sep = " ")
system(command_multiqc)