#!/bin/bash
# Align truncated pb reads to reference with minimap2

# Activate bioconda env
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ddocent_env

# Initial variables
t=''
in_t_fa=''
pre=''
res_dir=''
ref=''

# Set help info
help_inf="
	
	Usage: mm2_pbalign.sh -i <input> -r <reference> -p <prefix> -d <result_dir> -t <threads>
		
		-h|--help	Prints out this dialogue
		
		-i|--input	full path to truncated telomere reads
		
		-p|--prefix     Character prefix for file names, ex. \"p_falciparum_telomeres\"
		
		-r|--reference	Reference assembly
		
		-d|--result_dir	full path to project results direcotry
		
		-t|--threads	integer for number of threads to use for operation

		**		all options are !required!

"

# Checkf for any args
if [[ -z $@ ]]; then
	echo -e "$help_inf"
	exit 1
fi

# Read options
ARGS=$( getopt -o h::i:p:r:d:t: -l "help::,input:,prefix:,reference:,result_dir:,threads" -n "mm2_pbalign.sh" -- "$@" );


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
                                in_t_fa=$1;
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
if [[ -z $in_t_fa || -z $pre || -z $ref || -z $res_dir || -z $t ]]; then
        echo -e "\nRequired options (-i,-p,-r,-d,-t) were not supplied, exiting\n"
        exit 1;
fi

# Set read out
readout_vars="

        This input trucated telomere fasta: $in_t_fa

        Prefix: $pre

        Reference Assembly: $ref

        Results directory: $res_dir

        Number of threads: $t

        "
# Echo vars for authentication
echo -e "\nThese are the variables:
        ${readout_vars}
        "


# Sync data and setup defaults
threads=${t}
export work_dir=${res_dir}
mkdir ${work_dir}/output
export output=${work_dir}/output

# Run minimap2 and convert to bam
echo -e "\n\nStart minimap2 of trunc at `date`."
minimap2 -t ${threads} -ax map-pb ${ref} ${in_t_fa} > ${output}/${pre}.alignment.sam

echo -e "\n\nStart conversion of trunc to bam, sort and index at `date`."
samtools view -S -b ${output}/${pre}.alignment.sam > ${output}/${pre}.alignment.bam
samtools sort ${output}/${pre}.alignment.bam -o ${output}/${pre}.alignment.sorted.bam
samtools index ${output}/${pre}.alignment.sorted.bam

conda deactivate

# Print out bam stats
echo -e "\n\nCoverage and stats follow"
samtools coverage ${output}/${pre}.alignment.sorted.bam 

echo -e "\n\n"
samtools stats ${output}/${pre}.alignment.sorted.bam > ${output}/${pre}.alignment.stats

echo -e "\n\nJOBS DONE..... at `date`."
