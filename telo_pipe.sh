#!/bin/bash
# Pipeline to run all scripts for telomere processing

# Initialize variables
in_dir=''
pre=''
ref=''
res_dir=''
sn=''
bf=''
med=''
t=''
infq_1=''
infq_2=''

cur_dir=`pwd`
bin_path=${cur_dir}/bin/

# Set help dialog
help_inf="

        Usage: telo_pipe.sh -i <work_dir> -a <infq1> -f <infq2> -r <reference> -p <prefix> 
			    -d <result_dir> -s <sample_names> -m <medians> -t <threads>

                -h|--help       	Prints out this dialogiue
		
		-w|--work_dir		Input working directory

                -a|--infq1      	Input fastq1 
		
		-f|--infq2		Input fastq2

                -p|--prefix     	Character prefix for file names, ex. \"p_falciparum_telomeres\"

                -r|--reference  	Reference assembly

                -d|--result_dir 	full path to project results direcotry

		-s|--sample_names	samples names, ex. \"wt,irr,KO,...\"

		-m|--medians		medians for samples, ex. \"5346,7895,9058...\"

                -t|--threads    	integer for number of threads to use for operation, DEFAULT = max-2

                **              	i,p,r,d,s,b,m are !required!

"

# Checkf for any args
if [[ -z $@ ]]; then
        echo -e "$help_inf"
        exit 1
fi

# Read options
ARGS=$( getopt -o h::w:a:f:p:r:d:s:m:t: -l "help::,work_dir:,infq1:infq2:,prefix:,reference:,result_dir:,sample_names:,medians:,threads" -n "telo_pipe.sh" -- "$@" );


eval set -- "$ARGS";

# extract options
while true; do
        case "$1" in
                -h|--help)
                        shift;
                        echo -e "${help_inf}";
                        exit 1;
                ;;
                -w|--work_dir)
                        shift;
                        if [[ -n $1 ]]; then
                                w_dir=$1;
                                shift;
                        fi
                ;;
		-a|--infq1)
                        shift;
                        if [[ -n $1 ]]; then
                                infq_1=$1;
                                shift;
                        fi
                ;;
		-f|--infq2)
                        shift;
                        if [[ -n $1 ]]; then
                                infq_2=$1;
                                shift;
                        fi
                ;;
                -p|--prefix)
                        shift;
                        if [[ -n $1 ]]; then
                                pre=$1;
                                shift;
                        fi
                ;;
                -r|--reference)
                        shift;
                        if [[ -n $1 ]]; then
                                ref=$1;
                                shift;
                        fi
                ;;
                -d|--result_dir)
                        shift;
                        if [[ -n $1 ]]; then
                                res_dir=$1;
                                shift;
                        fi
                ;;
		-s|--sample_names)
                        shift;
                        if [[ -n $1 ]]; then
                                sn=$1;
                                shift;
                        fi
                ;;
		-m|--medians)
                        shift;
                        if [[ -n $1 ]]; then
                                med=$1;
                                shift;
                        fi
                ;;
                -t|--threads)
                        shift;
                        if [[ -n $1 ]]; then
                                t=$1;
                                shift;
                        fi
                ;;
                --)
                shift;
                break;
                ;;
        esac
done

# Check required arguements are met
if [[ -z $w_dir || -z $pre || -z $ref || -z $res_dir \
|| -z $sn || -z $med || -z $infq_1} || -z $infq_2} ]]; then
        echo -e "\nRequired options (-w,-a,-f,-p,-r,-d,-s,-m) were not supplied, exiting\n"
        exit 1;
fi

# Check if threads set
if [[ -z $t ]]; then 
	echo -e "\nThreads not set, setting to max-2."
	t=$( echo `nproc --all`-2 | bc )
fi

# Echo initial variables and run telo count script
initial_vars="
	
	Working Directory: $w_dir

        Input fastq1: $infq_1

        Input fastq2: $infq_2

        Reference Assembly: $ref

	Threads: $t	
"

echo -e "\nThese are initial variables for telo csv script:
	${initial_vars}
	"

# Create telo length files
# Check if telo length file already exists
telo_csv=$( find ${w_dir} -name "200sw_telomere_ranges.perc.sorted.csv" )
if [[ -f ${telo_csv} ]]; then
        echo -e "\nTelo csv file already exists, continuing\n"

else
        echo -e "\nTelo csv file does not exists, starting telo script\n"
        ${bin_path}/telo_homer.sh -s ${infq_1} -a ${infq_2} -w ${w_dir} -r ${ref} -t ${t}

fi

# Assign files based on input directory 
in_fa_s=$( find ${w_dir}/output -name "*.sample.*" | grep -v "lenStats" | grep -v "slide" )
pref=$( echo ${in_fa_s} | xargs -L 1 basename | cut -d'.' -f1 )
in_fa=$( find ${w_dir}/output -name "*.fasta" | grep -v ${pref} )

# Set read out
readout_vars="

        Directory for telomere table results: $w_dir

	Telomere CSV: $telo_csv

	Input non-sample fasta: $in_fa

	Input sampled fasta: $in_fa_s        

	Prefix: $pre

        Reference Assembly: $ref

        Results directory: $res_dir

        Sample names: $sn

	Median values for normalization: $med

	Number of threads: $t

	Bin path: $bin_path

        "

# Echo vars for authentication
echo -e "\nThese are the variables:
        ${readout_vars}
        "
# Checking if results directory exists
if [[ -d ${res_dir} ]]; then
	echo -e "\nResult directory exists, continuing\n"

else
	echo -e "\nResult directory does not exist, creating it\n"
	mkdir ${res_dir}
fi 

# Checking for output directory
if [[ -d ${res_dir}/output ]]; then
	echo -e "\nOutput directory exists, continuing\n"

else
	echo -e "\nOutput directory does not exist, creating it\n"
	mkdir ${res_dir}/output

fi

# Telo length assignment
# Check if df file already created
df=$( find ${res_dir} -type f -name "${pre}.df.Rda" );
if [[ -f ${df} ]]; then 
	echo -e "\nDF file present, skipping length assessment."

else 
	echo -e "\nStarting telomere length assessment at `date`\n"
	${bin_path}/mod1_ln.R -i ${telo_csv} -o ${res_dir} -p ${pre} -m ${med} -r ${sn} -t ${t}

fi

# Telomere truncation
# Check if truncation file already created
trunc=$( find ${res_dir} -type f -name "${pre}.dna.trunc.fa.gz" );
if [[ -f ${trunc} ]]; then 
	echo -e "\nTruncation file already created, skipping."

else
	echo -e "\nStarting truncation of telomeres at `date`\n"
	${bin_path}/mod2_trunc.R -i ${in_fa} -s ${in_fa_s} -o ${res_dir} -p ${pre} -t ${t}

fi

# Minimap 2 alignment
# Check if alignment already created
align=$( find ${res_dir}/output -type f -name "${pre}.alignment.sorted.bam" );
if [[ -f ${align} ]]; then 
	echo -e "\nAlignment file already created, skipping."

else
	echo -e "\nStarting to align truncated reads to reference at `date`\n"
	${bin_path}/mm2_pbalign.sh -i ${res_dir}/${pre}.dna.trunc.fa.gz -p ${pre} -r ${ref} -d ${res_dir} -t ${t}

fi

# Assigning endedness based on minimap alignment
# Check if endednesss already assigned
end=$( find ${res_dir} -type f -name "${pre}.end.bam.csv" );
if [[ -f ${end} ]]; then 
	echo -e "\nEnd bam csv already created, skipping."

else
	echo -e "\nAssigning endedness at `date`\n"
	${bin_path}/mod3_bam.R -i ${res_dir}/output/${pre}.alignment.sorted.bam -o ${res_dir} -p ${pre} -t ${t}

fi

# Graph results
echo -e "\nGraphing results at `date`\n"
${bin_path}/mod4_graphing.R -r ${res_dir} -p ${pre} 

echo -e "JOBS DONE.... at `date`."
