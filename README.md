# TLD -- Telomere Length Determination

## Table of Contents
* [General Information](#info)
* [Installation](#install)
* [Usage](#usage)
* [Example -- Comparing Two Yeast Samples](#example)
* [Citing TLD](#cite)
* [References](#ref)

## <a name="info"></a>General Information
This application is used to determine telomere length distributions, given two sequencing samples and compare the length distributions overall and by chromosome of those samples. The application requires two fastq files and a reference fasta file. Given these three files telomere length distributions can be determined, denovo telomere motifs are found, chromosomal telomere read counts and telomere length distributions. 

There is a caveat in that the mean read length of the sequencing set must be longer than the predicted telomere length. For example, *Homo sapien's* telomeres are estimated to be ~5-15kb long; therefore, to accurately determine telomere length of *Homo sapien* sample, the mean read length of the sequencing sets must be >15kb in length.

## <a name="install"></a>Installation
To run this application, install docker and download the docker image

### Install docker [docker]: https://docs.docker.com/get-docker/

#### Install docker on ubuntu 20.04
```sh
sudo apt-get update
sudo apt-get install \
    ca-certificates \
    curl \
    gnupg \
    lsb-release

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io

```

#### Download TLD docker image
```sh
sudo docker pull jreed0pbsb/tld:0.02
```

#### Clone TLD
```sh
mkdir <directory_to_install_tld>
git clone https://github.com/jake-bioinfo/tld.git <directory_to_install_tld>
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


#### Pull sra-tools docker image
```sh
sudo docker pull ncbi/sra-tools
```

#### Setup docker image and download yeast strains
```sh
# if executed within cloned directory
sudo docker run -id --name sra -v ./data:/dna ncbi/sra-tools:latest
# if cloned directory is in $HOME
# sudo docker run -id --name sra -v $HOME/tld/data:/dna ncbi/sra-tools:latest
sudo docker exec -it sra fastq-dump -v SRR13577847 -O /dna
# if executed within cloned directory
sudo mv ./data/SRR13577847.fastq ./data/s288c.fastq
# if cloned directory is in $HOME
# sudo mv $HOME/tld/data/SRR13577847.fastq $HOME/tld/data/s288c.fastq
sudo docker exec -it sra fastq-dump -v SRR13577846 -O /dna
# if executed within cloned directory
sudo mv ./data/SRR13577846.fastq ./data/cen-pk.fastq
# if cloned directory is in $HOME
# sudo mv $HOME/tld/data/SRR13577846.fastq $HOME/tld/data/cen-pk.fastq
sudo docker container stop sra
sudo docker container rm sra
```

#### Get and unpack reference genome
```sh
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz -O ./S288C_current_ref.tgz
tar -xvf S288C_current_ref.tgz
sudo mv S288C_ref* S288C_current_ref
gzip -d S288C_current_ref/*.fsa.gz
# if executed within cloned directory
sudo mv S288C_current_ref/*.fsa ./data/S288C_ref_genome.fasta
# if cloned directory is in $HOME
# sudo mv S288C_current_ref/*.fsa ./data/S288C_ref_genome.fasta
rm -rf S288C_current_ref
rm S288C_current_ref.tgz
```

#### Setup tld docker image and execute tld command
```sh
# if cloned directory is in $HOME
# sudo docker run -id --name tld -v $HOME/tld:/tld jreed0pbsb/tld:0.02
# if executed within cloned directory
sudo docker run -id --name tld -v ./:/tld jreed0pbsb/tld:0.02
# Print help
sudo docker exec -it tld /tld/telo_pipe.sh -h
# Set threads, default is nproc - 1
threads=$(echo "$(nproc) - 1" | bc)
sudo docker exec -it tld /tld/telo_pipe.sh -w /tld/data/w_dir -d /tld/data/o_dir \
	-a /tld/data/s288c.fastq \
	-f /tld/data/cen-pk.fastq \
	-r /tld/data/S288C_ref_genome.fasta \
	-p yeast \
	-s "s288c,cen-pk" \
	-m "1,1" -j 6 -l 100 -t ${threads}
```

## <a name="cite"></a>Citing TLD
1. Reed J, Kirkman LA, Kafsack BF, Mason CE, Deitsch KW. Telomere length dynamics in response to DNA damage in malaria parasites. iScience. 2021 Jan 20;24(2):102082. doi: 10.1016/j.isci.2021.102082. PMID: 33644714; PMCID: PMC7887396.

## <a name="ref"></a>References
1. Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432
2. W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE. doi:10.1371/journal.pone.0163962.
3. Sathish-deevi, modes-Package, (2016), Github repository, https://github.com/sathish-deevi/modes-Package
4. R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
5. Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
6. GNU, P. (2007). Free Software Foundation. Bash (3.2. 48)[Unix shell program].
7. Wickham H (2011). “The Split-Apply-Combine Strategy for Data Analysis.” Journal of Statistical Software, 40(1), 1–29. http://www.jstatsoft.org/v40/i01/.
8. Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2018). dplyr: A Grammar of Data Manipulation. R package version 0.7.6. https://CRAN.R-project.org/package=dplyr
9. Microsoft and Steve Weston (2020). foreach: Provides Foreach Looping Construct. R package version 1.5.1. https://CRAN.R-project.org/package=foreach
10. Microsoft Corporation and Steve Weston (2020). doParallel: Foreach Parallel Adaptor for the 'parallel' Package. R package version 1.0.16. https://CRAN.R-project.org/package=doParallel
11. H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2019). Biostrings: Efficient manipulation of biological strings. R package version 2.54.0.
12. Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. doi:10.1371/journal.pcbi.1003118GenomicAlignments 
13. Charif, D. and Lobry, J.R. (2007). seqinr. 
14. Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer
15. Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra
16. Claus O. Wilke (2020). cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'. R package version 1.1.1. https://CRAN.R-project.org/package=cowplot
