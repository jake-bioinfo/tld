# TLD -- Telomere Length Determination

## Table of Contents
* [General Information](#info)
* [Installation](#install)
* [Usage](#usage)
* [Example -- Comparing Two Yeast Samples](#example)
* [Citing TLD](#cite)
* [References](#ref)

## <a name="info"></a>General Information
This application is used to determine telomere length distributions, given two sequencing samples and compare the length distributions of those samples. The application requires two fastq files and a reference fasta file. Given these three files telomere length distributions can be determined, denovo telomere motifs are found and chromosomal telomere read counts. 

There is a caveat in that the mean read length of the sequencing set must be longer than the predicted telomere length. For example, homosapien's telomeres are estimated to be ~20-30kb long; therefore, to accurately determine telomere length of homosapien's, the mean read length of the sequencing sets must be >30kb in length.

## <a name="install"></a>Installation
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

## <a name="usage"></a>Usage
In order to use this application you must first move your data files into the data directory of tld. The fastq files go into ```$ tld/data/fastq ``` and the reference files must go in ```$ tld/data/ref ```. 

```
        Usage: telo_pipe.sh -w <work_dir> -a <infq1> -f <infq2> -r <reference> -p <prefix> 
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

## <a name="example"></a>Example -- Comparing 2 Yeast Samples

### Install SRA-Toolkit conda

```sh
conda install -c bioconda sra-tools
```

or

### Install SRA-Toolkit, Download Yeast Samples and Run TLD docker
For example sake, TLD is installed in $HOME/tld

```sh
# Pull sra-tools docker image
docker pull ncbi/sra-tools

# Setup docker image and download yeast strains
docker run -id --name sra -v $HOME/tld/data:/dna ncbi/sra-tools:latest
docker exec -it --rm sra fastq-dump -v SRR13577847 -O /dna
mv $HOME/tld/data/SRR13577847.fastq $HOME/tld/data/s288c.fastq
docker exec -it --rm sra fastq-dump -v SRR13577846 -O /dna
mv $HOME/tld/data/SRR13577846.fastq $HOME/tld/data/cen-pk.fastq
docker container stop sra

# Get and unpack reference genome
wget -P $HOME/tld/data http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz
tar -xvf $HOME/tld/data/S288C_reference_genome_Current_Release.tgz -C $HOME/tld/data
gzip -d $HOME/tld/data/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa.gz
mv $HOME/tld/data/S288C_reference_genome_R64-3-1_20210421/S288C_reference_sequence_R64-3-1_20210421.fsa $HOME/tld/data/S288C_ref_genome.fasta
rm -rf $HOME/tld/data/S288C_reference_genome_R64-3-1_20210421
rm $HOME/tld/data/S288C_reference_genome_Current_Release.tgz

# Setup tld docker image and execute tld command
docker run -id --name tld -v $HOME/tld:/tld jreed0pbsb/tld:latest
docker exec -it --rm tld /tld/telo_pipe.sh -w /tld/data/w_dir -o /tld/data/o_dir \
	-a /tld/data/s288c.fastq \
	-f /tld/data/cen-pk.fastq \
	-r /tld/data/S288C_ref_genome.fasta \
	-p yeast \
	-s "s288c,cen-pk" \
	-m "1,1" -j 6 -l 100 -t 7

```

## <a name="cite"></a>Citing TLD
1. Reed J, Kirkman LA, Kafsack BF, Mason CE, Deitsch KW. Telomere length dynamics in response to DNA damage in malaria parasites. iScience. 2021 Jan 20;24(2):102082. doi: 10.1016/j.isci.2021.102082. PMID: 33644714; PMCID: PMC7887396.

## <a name="ref"></a>References
1. HOMER
2. SeqKit
3. SRA-ToolKit
4. Modes
5. R
6. ggplot2
7. BASH


