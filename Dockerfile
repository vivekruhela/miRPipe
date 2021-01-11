FROM jupyter/scipy-notebook

USER root

MAINTAINER Vivek Ruhela <vivekr@iiitd.ac.in>

# Copy the application folder inside the container
ADD . /MiRPipe_Docker

# Set the default directory where CMD will execute
WORKDIR /MiRPipe_Docker

# Set environment variable
# ENV http_proxy 'proxy.com'
ENV HOME /MiRPipe_Docker
ENV LC_ALL C
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

# Install dependencies like wget, unzip, java & useful packages like fastqc, samtools
RUN apt-get update -qq && apt-get install -y \            
      wget \
      curl \
      unzip \
      texlive-xetex \ 
      build-essential \
      default-jre \
      gcc \
      make \
      pigz \
      libbz2-dev \
      zlib1g-dev \
      libncurses5-dev \
      libncursesw5-dev \
      liblzma-dev \
      libcurl4-openssl-dev \
      libpq-dev \
      libssl-dev \
      python3 python3-pip python3-dev build-essential \
      bedtools \
      fastqc \
      pandoc \
      samtools

# Install additional python packages
RUN pip install --upgrade pip
RUN pip install --user -r requirements.txt

# Install multiqc
RUN pip3 install --user --upgrade multiqc && ln -sfn ~/.local/bin/multiqc /usr/bin/
RUN pip3 install --user --upgrade cutadapt && ln -sfn ~/.local/bin/cutadapt /usr/bin/
RUN pip3 install --user --upgrade pandoc


# Install R, DESeq2

RUN apt-get update && apt-get -y install --no-install-recommends --no-install-suggests \
       ca-certificates software-properties-common gnupg2 gnupg1 \
      && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
      && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
      && apt-get -y install r-base 

RUN apt-get install --yes libxml2-dev
RUN echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.rstudio.com"; options(repos=r)})' > ~/.Rprofile
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install(c("DESeq2","ShortRead","Biostrings"))'

# Download Bowtie Source Codes
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3/bowtie-1.2.3-linux-x86_64.zip -P $HOME/Tools
RUN unzip -o $HOME/Tools/bowtie-1.2.3-linux-x86_64.zip -d $HOME/Tools

# Install Mirdeep* (v38) and miRBase_v22
RUN wget --content-disposition  https://sourceforge.net/projects/mirdeepstar/files/MDS_command_line_v38.zip/download -O $HOME/Tools/MDS_command_line_v38.zip && unzip -o $HOME/Tools/MDS_command_line_v38.zip -d $HOME/Tools
RUN wget --content-disposition  https://sourceforge.net/projects/mirdeepstar/files/Index_files/hg38.zip/download -O $HOME/Tools/hg38.zip && unzip -o $HOME/Tools/hg38.zip -d $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome
RUN rm -r $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome/hg19
RUN rm $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase/*
RUN wget --content-disposition  ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz -P $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase
RUN pigz -p 5 -d $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase/hairpin.fa.gz
RUN wget --content-disposition  ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz -P $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase
RUN pigz -p 5 -d $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase/mature.fa.gz
RUN wget --content-disposition  ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3 -P $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase
RUN mv $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase/hsa.gff3 $HOME/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase/knownMiR.gff3

# Download Trim_Galore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.1.tar.gz -P $HOME/Tools && tar xvzf $HOME/Tools/0.6.1.tar.gz -C $HOME/Tools


# Install Fastx Tool
RUN wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 -P $HOME/Tools
RUN tar xjf $HOME/Tools/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 -C $HOME/Tools
RUN mv $HOME/Tools/bin $HOME/Tools/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64

# Install CD-HIT
RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz -P $HOME/Tools
RUN tar xvzf $HOME/Tools/cd-hit-v4.6.8-2017-1208-source.tar.gz -C $HOME/Tools
RUN make -C $HOME/Tools/cd-hit-v4.6.8-2017-1208

# Install RNAFold
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz -P $HOME/Tools
RUN tar xvzf $HOME/Tools/ViennaRNA-2.4.14.tar.gz -C $HOME/Tools  
WORKDIR $HOME/Tools/ViennaRNA-2.4.14
RUN ./configure
RUN make && make check && make install
WORKDIR $HOME


# Download BBMap Source Codes
RUN wget https://sourceforge.net/projects/bbmap/files/latest/download -O $HOME/Tools/bbmap.tar.gz
RUN tar xvzf $HOME/Tools/bbmap.tar.gz -C $HOME/Tools && rm $HOME/Tools/bbmap.tar.gz

# Install Firefox, geckodriver

# Reference: https://github.com/burningion/firefox-splinter-docker

RUN add-apt-repository -y ppa:mozillateam/firefox-next
RUN apt-get update && apt-get install -y firefox
RUN wget https://github.com/mozilla/geckodriver/releases/download/v0.26.0/geckodriver-v0.26.0-linux64.tar.gz && tar zxvf geckodriver-v0.26.0-linux64.tar.gz

# Clean-ups
RUN rm $HOME/Tools/*.zip && rm $HOME/Tools/*.gz
RUN rm $HOME/*.gz

# Expose port
EXPOSE 8888
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

ADD notebook.sh /

# Start notebook server
CMD ["/notebook.sh"]
