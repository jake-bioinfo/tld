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

cur_dir=`pwd`
bin_path=${cur_dir}/bin/

# Set help dialog
help_inf="

        Usage: telo_pipe.sh -i <input> -r <reference> -p <prefix> -d <result_dir> -s <sample_names> 
			    -m <medians> -t <threads>

                -h|--help       	Prints out this dialogue

                -i|--input      	telomere results directory

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
ARGS=$( getopt -o h::i:p:r:d:s:m:t: -l "help::,input:,prefix:,reference:,result_dir:,sample_names:,medians:,threads" -n "telo_pipe.sh" -- "$@" );


eval set -- "$ARGS";

# extract options
while true; do
        case "$1" in
                -h|--help)
                        shift;
                        echo -e "${help_inf}";
                        exit 1;
                ;;
                -i|--input)
                        shift;
                        if [[ -n $1 ]]; then
                                in_dir=$1;
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
if [[ -z $in_dir || -z $pre || -z $ref || -z $res_dir || -z $sn || -z $med ]]; then
        echo -e "\nRequired options (-i,-p,-r,-d,-s,-m) were not supplied, exiting\n"
        exit 1;
fi

# Check if threads set
if [[ -z $t ]]; then 
	echo -e "\nThreads not set, setting to max-2."
	t=$( echo `nproc --all`-2 | bc )
fi

# Assign files based on input directory 
telo_csv=$( find ${in_dir} -name "200sw_telomere_ranges.perc.sorted.csv" )
in_fa_s=$( find ${in_dir}/output -name "*.sample.*" | grep -v "lenStats" )
pref=$( echo ${in_fa_s} | xargs -L 1 basename | cut -d'.' -f1 )
in_fa=$( find ${in_dir}/output -name "*.fasta" | grep -v ${pref} )


# Set read out
readout_vars="

        Directory for telomere table results: $in_dir

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
	echo -e "\nResult directory does not exists, creating\n"
	mkdir ${res_dir}
fi 

# Telo length assignment
echo -e "\nStarting telomere length assessment at `date`\n"
${bin_path}/mod1_ln.R -i ${telo_csv} -o ${res_dir} -p ${pre} -m ${med} -r ${sn} -t ${t}

# Telomere truncation
echo -e "\nStarting truncation of telomeres at `date`\n"
${bin_path}/mod2_trunc.R -i ${in_fa} -s ${in_fa_s} -o ${res_dir} -p ${pre} -t ${t}

# Minimap 2 alignment
echo -e "\nStarting to align truncated reads to reference at `date`\n"
${bin_path}/mm2_pbalign.sh -i ${res_dir}/${pre}.dna.trunc.fa.gz -p ${pre} -r ${ref} -d ${res_dir} -t ${t}

# Assigning endedness based on minimap alignment
echo -e "\nAssigning endedness at `date`\n"
${bin_path}/mod3_bam.R -i ${res_dir}/output/${pre}.alignment.sorted.bam -o ${res_dir} -p ${pre} -t ${t}

# Graph results
echo -e "\nGraphing results at `date`\n"
${bin_path}/mod4_graphing.R -r ${res_dir} -p ${pre} -o ${res_dir} 

echo -e "JOBS DONE.... at `date`."
