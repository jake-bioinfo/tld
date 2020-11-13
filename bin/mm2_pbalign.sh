#!/bin/bash
# Align truncated pb reads to reference with minimap2

## Add pb or ont platform option and modify mm cmd

# Initial variables
t=''
in_t_fa=''
pre=''
res_dir=''
ref=''
plat=''

# Set help info
help_inf="
	
	Usage: mm2_pbalign.sh -i <input> -r <reference> -p <prefix> -d <result_dir> -t <threads> -n <platform>
		
		-h|--help	Prints out this dialogue
		
		-i|--input	full path to truncated telomere reads
		
		-p|--prefix     Character prefix for file names, ex. \"p_falciparum_telomeres\"
		
		-r|--reference	Reference assembly
		
		-d|--result_dir	Full path to project results direcotry

		-n|--platform	Select platform for experiment, pb or ont, default is pb
		
		-t|--threads	integer for number of threads to use for operation

		**		all options are except -n !required!

"

# Checkf for any args
if [[ -z $@ ]]; then
	echo -e "$help_inf"
	exit 1
fi

# Read options
ARGS=$( getopt -o h::i:p:r:d:t:n: -l "help::,input:,prefix:,reference:,result_dir:,threads:,platform:" -n "mm2_pbalign.sh" -- "$@" );


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
                -n|--platform)
                        shift;
                        if [[ -n $1 ]]; then
                                plat=$1;
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

# Check platform
if [[ $plat = "pb" ]]; then
        echo -e "\nPlatform selected is PacBio.\n"

elif [[ $plat = "ont" ]]; then
        echo -e "\nPlatform selected is Oxford Nanopore.\n"

elif [[ -z $plat ]]; then
        echo -e "\nPlatform not selected, defaulting to Pacbio.\n"
        plat='pb'

fi

# Set read out
readout_vars="

        This input trucated telomere fasta: $in_t_fa

        Prefix: $pre

        Reference Assembly: $ref

        Results directory: $res_dir

        Number of threads: $t

	Platform: $plat

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
echo -e "\n\nStart minimap2 of trunc data at `date`."

if [[ $plat == "ont" ]]; then
	minimap2 -t ${threads} -aY -x asm20 ${ref} ${in_t_fa} > ${output}/${pre}.alignment.sam
else
	minimap2 -t ${threads} -ax map-${plat} ${ref} ${in_t_fa} > ${output}/${pre}.alignment.sam

fi

echo -e "\n\nStart conversion of trunc to bam, sort and index at `date`."
samtools view -S -b ${output}/${pre}.alignment.sam > ${output}/${pre}.alignment.bam
samtools sort ${output}/${pre}.alignment.bam -o ${output}/${pre}.alignment.sorted.bam
samtools index ${output}/${pre}.alignment.sorted.bam

# Print out bam stats
echo -e "\n\nCoverage and stats follow"
samtools coverage ${output}/${pre}.alignment.sorted.bam 

echo -e "\n\n"
samtools stats ${output}/${pre}.alignment.sorted.bam > ${output}/${pre}.alignment.stats

echo -e "\n\nJOBS DONE..... at `date`."
