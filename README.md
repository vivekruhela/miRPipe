# MiRPipe: Open source RNA-Seq bioinformatics analysis docker for identification of miRNAs and piRNAs

![Graphical Abstract of miRPipe](miRPipe_Flowchart.png)

MiRPipe is the integrated, user-friendly jupyter notebook based RNA-Seq bioinformatics analysis docker. The aim of this docker is to help researchers perform miRNA identification and piRNA identification independently and with much ease.

## Running the RNA-SEQ Pipeline

We have created a Docker image packaging all the dependencies (command line tools, R and Python packages) for the pipeline, which are publically available. To run the Docker image on your local machine, you need to puoll the docker image. Please check the system requirements before installing the docker image.

### System Requirements

##### MAC Operating System:

System Requirements: Docker Desktop for Mac launches only if all of these requirements are met.

•Mac hardware must be a 2010 or newer model, with Intel’s hardware support for memory management unit (MMU) virtualization, including Extended Page Tables (EPT) and Unrestricted Mode. You can check to see if your machine has this support by running the following command in a terminal:
```
sysctl kern.hv_support
```
• macOS Sierra 10.12 and newer macOS releases are supported. We recommend upgrading to the latest version of macOS.

• At least 4GB of RAM

• VirtualBox prior to version 4.3.30 must NOT be installed (it is incompatible with Docker Desktop for Mac). If you have a newer version of VirtualBox installed, it’s fine.

##### WINDOWS Operating System:

System Requirements:

• Windows 10 64bit: Pro, Enterprise or Education (Build 15063 or later).

• Virtualization is enabled in BIOS. Typically, virtualization is enabled by default. This is different from having Hyper-V enabled. For more detail see Virtualization must be enabled in Troubleshooting.

• CPU SLAT-capable feature.

• At least 4GB of RAM.

##### LINUX Operating System:

System Requireents:

• 64bit, 8.00 GB RAM [OS version : Ubuntu 18.04 as MiRPipe has been developed in Ubuntu:18.04]

### Prerequisites

What things you need to install before installing the software-
##### For LINUX
a) docker
```
sudo apt-get update
sudo apt-get install docker-ce
```

##### For WINDOWS
Follow guidelines on this [webpage](https://docs.docker.com/v17.12/docker-for-windows/install/#install-docker-for-windows-desktop-app)  to install.


##### For Mac
Follown guidelines on this [webpage](https://docs.docker.com/v17.12/docker-for-mac/install/) to install.



Create an account on docker through the following [webpage](https://docs.docker.com/docker-id/), to run it


### Build Docker
After installing the docker, docker build is required to get the docker running.

To build, run the follwing command

```

docker pull docker.io/vivekruhela/mirpipe

```

# Execute docker

Before running docker, prepare the file `sample_list.csv` in the data directory in the following format:
```
sample           file          condition
Sample_A   Sample_A.fastq.gz    treated
Sample_B   Sample_B.fastq.gz    treated
Sample_C   Sample_C.fastq.gz    treated
   .              .                .
   .              .                .
   .              .                .
```
In `sample_list.csv`, all the sample name should be sorted in ascending sorting order. This is necessary because the final count file generated at the end of pipeline (and before differential expression analysis) have all the columns sorted in alpha-numerically ascending order. The order of columns in final count matrix generated in the end of pipeline and order in samples in `sample_list.csv` must be same in order to conduct correct differential expression analysis.

To run the docker with administration rights, run the following command at terminal.

```
docker run -p 8880:8888 \
           -e 'PASSWORD=password' \
           -e 'USE_HTTP=1' \
           -v /host/path/to/data:/MirPipe_Docker/data docker.io/vivekruhela/mirpipe

```

Explanation of execution step-

sudo docker run - Runs the docker with administrative rights.

-p 8880:8888 - Defines the port of the server and the port at which the docker will run.

-e 'PASSWORD=password' - Defines the password as "password", it can be changed as required by the user.

-e 'USE_HTTP=1' - To run in HTTP, we use USE_HTTP environment variable. Setting it to a non-zero value enables HTTP.

-v  /host/path/to/data:/MirPipe_Docker/data - Mount the host folder which contains raw fastq sequence file to data folder of MiRPipe docker

```
Example:

docker run -p 8880:8888 -e 'PASSWORD=password' -e 'USE_HTTP=1' -v /home/vivek/Small_fastq:/MiRPipe_Docker/data vivekruhela/mirpipe
```
After successfully running the docker execution command, open firefox browser, and type the follwing in the URL:

```
localhost:8880/mirpipe
```
Jupyter Notebook environment will be opened in Firefox browser.
Enter the PASSWORD mentioned in above command, by default it is set to- password
