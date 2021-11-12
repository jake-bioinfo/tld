# TLD -- Telomere Length Determination

## Table of Contents
* [General Information](#general-info)
* [Installation](#install)
* [Usage](#usage)
* [Citing TLD](#cite)
* [References](#ref)

## General Information
This application is used to determine telomere length distributions, given two sequencing samples and compare the length distributions of those samples. The application requires two fastq files and a reference fasta file. Given these three files telomere length distributions can be determined, denovo telomere motifs are found and chromosomal telomere read counts. 

There is a caveat in that the mean read length of the sequencing set must be longer than the predicted telomere length. For example, homosapien's telomeres are estimated to be ~20-30kb long; therefore, to accurately determine telomere length of homosapien's, the mean read length of the sequencing sets must be >30kb in length.

## Installation
To run this application, install docker and download the docker image

### Install docker [docker]: https://docs.docker.com/get-docker/

#### Install docker on ubuntu 20.04
```sh
$ sudo apt-get update
$ sudo apt-get install \
    ca-certificates \
    curl \
    gnupg \
    lsb-release

$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

$ echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

$ sudo apt-get update
$ sudo apt-get install docker-ce docker-ce-cli containerd.io

```

#### Download TLD docker image
```sh
$ docker pull jreed0pbsb/tld:latest
```

## Usage
In order to use this application you must first move your data files into the data directory of tld. The fastq files go into ```$ tld/data/fastq ``` and the reference files must go in ```$ tld/data/ref ```. 

```
        Usage: telo_pipe.sh -i <work_dir> -a <infq1> -f <infq2> -r <reference> -p <prefix> 
			    -d <result_dir> -s <sample_names> -m <medians> -n <platform> -t <threads>
			    -j <estimated_telomere_motif_length>

                -h|--help       	Prints out this dialogiue
		
		-w|--work_dir		Input working directory

                -a|--infq1      	Input fastq1 
		
		-f|--infq2		Input fastq2

                -p|--prefix     	Character prefix for file names, ex. "p_falciparum_telomeres"

                -r|--reference  	Reference assembly

                -d|--result_dir 	Full path to project results direcotry

		-n|--platform		Select sequencing platform, ont or pb, pb is default

		-s|--sample_names	samples names, ex. "wt,irr,KO,..."

		-m|--medians		medians for samples, ex. "5346,7895,9058..."

		-j|--motif		predicted motif length

		-l|--end_length		Length of end of chromosome for telomere analysis, DEFAULT = 200

		-e|--high_euk		Set for higher eukaryotes, DEFAULT = option unset (lower eukaryotes)

                -t|--threads    	integer for number of threads to use for operation, DEFAULT = max-2

                **              	-w,-p,-r,-d,-s,-m,-a,f,-j are !required!

```

### Comparing 2 yeast samples
#### Initialize docker container 
```sh
$ docker run -d --name <name_of_container> <name_of_img>
```
