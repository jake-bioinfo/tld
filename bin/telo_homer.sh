#!/bin/bash
# Pull telomeres from irradiated and non-irradiated pfalci, sort, count length compare
# For new Sequel 2 data

## Create cases for ont vs pb
## ont read structure: >525bd77d-0100-42b4-9c36-03e6b7f75624 runid=46f564a18d50efd876dc52d7bb1d=FAK79249 protocol_group_id=R51-B11 sample_id=R51-B11
## PB read structure: >m54334U_200503_055320/1141/ccs



# Set initial variables
s1=''
s2=''
ref=''
w_dir=''

# Set help dialog
help_inf="

        Usage: telo_homer.sh -s <sample_1> -a <sample_2> -w <work_dir> -r <reference> -t <threads>

                -h|--help               Prints help for telo_homer

                -s|--sample_1		Input sample 1 fastq (accepts fastq or fastq.gz)             

		-a|--sample_2		Input sample 2 fastq (accepts fastq or fastq.gz)
		
		-w|--work_dir		Working directory

                -r|--reference          Reference genome

                -t|--threads            integer for number of threads to use for operation, DEFAULT = max-2

                **                      -s,-a,-w,-r,-n !required!

"

# Chec for any args
if [[ -z $@ ]]; then
        echo -e "$help_inf"
        exit 1
fi

# Read options
ARGS=$( getopt -o h::s:a:w:r:t: -l "help::,sample_1:,sample_2:,work_dir:,reference:,threads" -n "telo_pipe.sh" -- "$@" );


eval set -- "$ARGS";

# extract options
while true; do
        case "$1" in
                -h|--help)
                        shift;
                        echo -e "${help_inf}";
                        exit 1;
                ;;
                -s|--sample_1)
                        shift;
                        if [[ -n $1 ]]; then
                                s1=$1;
                                shift;
                        fi
                ;;
                -a|--sample_2)
                        shift;
                        if [[ -n $1 ]]; then
                                s2=$1;
                                shift;
                        fi
                ;;
                -w|--work_dir)
                        shift;
                        if [[ -n $1 ]]; then
                                w_dir=$1;
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


#Check required arguements are met
if [[ -z $s1 || -z $s2 || -z $ref || -z ${w_dir} ]]; then
        echo -e "\nRequired options (-s,-a,-w,-r) were not supplied, exiting\n"
        exit 1;
fi

# Check for working directory
if [[ -d ${w_dir} ]]; then
	echo -e "\nWorking directory already created, continuing...\n"

else
	echo -e "\nWorking directory not created, creating it\n"
	mkdir ${w_dir}; cd ${w_dir}

fi

# Check if threads set
if [[ -z $t ]]; then
        echo -e "\nThreads not set, setting to max-2."
        t=$( echo `nproc --all`-2 | bc )
fi

# Set read out
readout_vars="

        Sample 1: $s1

        Sample 2: $s2

	Reference genome: $ref

        Number of threads: $t

        "

# Echo vars for authentication
echo -e "\nThese are the variables:
        ${readout_vars}
        "
# Source anaconda environment
conda_base=$( find ~/ -maxdepth 3 -type d -name "anaconda3" | head -n 1 ) 
source ${conda_base}/etc/profile.d/conda.sh

# Export directories and file paths
s1floc=${s1}
s2floc=${s2}
reffloc=${ref}
export TMP=`pwd`

# Set number of threads
threads=${t}

# Check input directory and files exists
check_in_dir=$( ls ${TMP} | egrep "input" );

if [[ -n ${check_in_dir} ]]; then
	echo -e "\nInput directory already created.\n"
	export in=${TMP}/input
else 
	echo -e "\nInput directory not created, creating and exporting.\n"
	mkdir ${TMP}/input
	export in=${TMP}/input
fi

# Peform input file checks
check_in_f1=$( ls ${in} | egrep "fastq" | head -n 1 );
check_in_f2=$( ls ${in} | egrep "fastq" | head -n 2 | tail -n 1 );
check_in_ref=$( ls ${in} | egrep "fasta" );
s1f_bn=$( basename -- ${s1floc} );
s2f_bn=$( basename -- ${s2floc} );
ref_bn=$( basename -- ${ref} ); 

# Check if file 1 is synced already
if [[ -f ${check_in_f1} ]]; then
	echo -e "\nInput file 1 already synced.\n"
        export ni_r=${in}/${s1f_bn}
else
        echo -e "\nInput file 1 not synced, syncing now.\n"
        rsync -av --progress ${s1floc} ${in}
	export ni_r=${in}/${s1f_bn}
fi

# Check if file 2 is synced already
if [[ -f ${check_in_f2} ]]; then
	echo -e "\nInput file 2 already synced.\n"
        export i_r=${in}/${s2f_bn}
else
        echo -e "\nInput file 2 not synced, syncing now.\n"
        rsync -av --progress ${s2floc} ${in}
	export i_r=${in}/${s2f_bn}
fi

# Check for ref file
if [[ -f ${check_in_ref} ]]; then
        echo -e "\nInput reference file already synced.\n"
        export ref=${in}/${ref_bn}
else
        echo -e "\nInput reference file not synced, syncing now.\n"
        rsync -av --progress ${reffloc} ${in}
	export ref=${in}/${ref_bn}
fi

# Get sample files basename
nir_bn=$( basename -- ${ni_r} | sed "s:.fastq::g" );
ir_bn=$( basename -- ${i_r} | sed "s:.fastq::g" );


# Check if output directory is already created
check_out_dir=$( ls ${TMP} | egrep "output" );
if [[ -n ${check_out_dir} ]]; then
        echo -e "\nOutput directory already created.\n"
        export out=${TMP}/output
else
        mkdir ${TMP}/output
        export out=${TMP}/output
fi

# Create tmp directory
[ -d ${TMP}/tmp ] && \
{ echo -e "\ntmp directory present, proceeding."; export t_dir=${TMP}/tmp; } || \
{ echo -e "\ntmp directory not present, creating it."; mkdir ${TMP}/tmp; export t_dir=${TMP}/tmp; }

# Run homer to search for motifs
# Activate homer conda environment
conda activate homer

# Create ref index if not already there
check_ref_i=$( ls ${in} | egrep "Pfalciparum.genome.fasta.fai" );
[ -f ${check_ref_i} ] && \
{ echo -e "\nReference not present, creating it."; samtools faidx ${ref}; } || \
{ echo -e "\nReference present, proceeding."; }

# Create telo bed file
echo -e "\nProcessing index file to create telo.bed.\n"
touch ${in}/telo.bed

i=1
while IFS= read -r l; do

	# Create labels
	lab_p_5="seg_${i}_5p"
	echo -e "\n\tThis is label for 5' + strand: ${lab_p_5}"

	lab_n_5="seg_${i}_5n"
	echo -e "\n\tThis is label for 5' - strand: ${lab_n_5}"

        lab_p_3="seg_${i}_3p"
        echo -e "\n\tThis is label for 3' + strand: ${lab_p_3}"

        lab_n_3="seg_${i}_5n"
        echo -e "\n\tThis is label for 3' - strand: ${lab_n_3}"
	
	# Read chromosome header
	chr_head=$( echo -e  "${l}" | cut -f1 ); 
	echo -e "\n\tThis is chromosome header: ${chr_head}"
	
	# Print 5' + lines to telo.bed file
	echo -e "${lab_p_5}\t${chr_head}\t0\t1000\t0" >> ${in}/telo.bed
	# Print 5' - lines to telo.bed file
	echo -e "${lab_n_5}\t${chr_head}\t0\t1000\t1" >> ${in}/telo.bed

	# Create 3' chromosome locations
	end_loc=$( echo -e "${l}" | cut -f2 );
	st_loc=$( echo "scale=1; ${end_loc} - 1000" | bc )	

	# Print 3' + lines to telo.bed file
	echo -e "${lab_p_3}\t${chr_head}\t${st_loc}\t${end_loc}\t0" >> ${in}/telo.bed
	# Print 3' - lines to telo.bed file
	echo -e "${lab_n_3}\t${chr_head}\t${st_loc}\t${end_loc}\t1" >> ${in}/telo.bed
	
	# Increase iterator
	let i++

done < ${ref}.fai

# Check if homer has been run, if not run it. 
check_homer=$( ls ${out} | egrep "homerMotifs.all.motifs" );
[ -f ${check_homer} ] && \
{ echo -e "\nHomer result file not found, run homer.";
 
	# Run Homer on reference
	findMotifsGenome.pl ${in}/telo.bed ${ref} \
	${out} -len 7 -noweight -S 15 -maxN 2 -preparse -p ${threads}; } || \
{ echo -e "\nHomer result file found, proceeding."; }

# Parse homer output
# Create temp.motifs file with the first 5 motifs by multiplicity
# Check if motifs file already created
check_motif=$( ls ${out} | egrep "temp.motifs" );
[ -f ${check_motif} ] && \
{ echo -e "\nMotifs file not found, creating them."; 
	egrep "^>" ${out}/homerMotifs.all.motifs | egrep "CC|GG" | sed "s:,:\t:g" \
	| sed "s:>::g" | sort -k14 -r | cut -f1 | head -n 5 > ${out}/temp.motifs

	# Create reverse and reverse complements of all motifs
	touch ${out}/temp_rc.motifs
	while IFS= read -r l; do 
		rev_comp=$( echo "${l}" | rev | tr "[ATGCatgc]" "[TACGtacg]" );
		rev=$( echo "${l}" | rev );
		echo -e "${l}\t${rev_comp}\t${rev}" >> ${out}/temp_rc.motifs;
		
	done < ${out}/temp.motifs; } || \

{ echo -e "\nMotifs file found, proceeding."; }

# Create initial grep pull down statement
while IFS= read -r l; do
	orig=$( echo "${l}" | cut -f1 )
	rev_comp=$( echo "${l}" | cut -f2 )
	rev=$( echo "${l}" | cut -f3 )
	grep_exp+="${orig}${orig}|${rev}${rev}|${rev_comp}${rev_comp}|"
	ct_grep+="${orig}|${rev}|${rev_comp}|" 
done < ${out}/temp_rc.motifs

# check grep expression
grep_exp=$( echo "${grep_exp}" | head -c -2 );
#grep_exp=$( echo -e "${grep_exp}" | head -c -2 ); 
echo -e "\nThis is grep expression: ${grep_exp}\n";
ct_grep=$( echo "${ct_grep}" | sed "s:.$::g" );
echo -e "\nThis is count grep expression: ${ct_grep}\n";

# Check for telomere reads combined file
check_telo=$( ls ${out} | egrep "telo.fastq" | wc -l );
if [[ "${check_telo}" == 2 ]]; then
	echo -e "\nTelo files exist, not creating them.\n"
else
	echo -e "\nTelo files do not exist, creating them. \n"

	#parallel -j ${threads} --pipe --block 10M LC_ALL=C \
	cat ${ni_r} | egrep -i -B1 -A2 "${grep_exp}" \
 		| sed '/^--$/d' > ${out}/${nir_bn}.telo.fastq

        #parallel -j ${threads} --pipe --block 10M LC_ALL=C \
        cat ${i_r} | egrep -i -B1 -A2 "${grep_exp}" \
                | sed '/^--$/d' > ${out}/${ir_bn}.telo.fastq

	# Converting to fasta
	seqkit fq2fa -j ${threads} ${out}/${nir_bn}.telo.fastq -o ${out}/${nir_bn}.telo.fasta
	seqkit fq2fa -j ${threads} ${out}/${ir_bn}.telo.fastq -o ${out}/${ir_bn}.telo.fasta

fi

# Check number of reads in each combined telo file and random downsample larger file, should be irrf12 bc it had a coverage +30X
if [[ -f ${out}/${nir_bn}.telo.fasta && -f ${out}/${ir_bn}.telo.fasta ]]; then
	echo -e "\nBoth combined telo files exist. Determining number of reads in each. \n"
	
	check_if12=$( egrep "^>" ${out}/${ir_bn}.telo.fasta | wc -l )
	check_ni3d7=$( egrep "^>" ${out}/${nir_bn}.telo.fasta | wc -l )
	echo -e "\nNumber of reads in each telo follows.\n\tif12: ${check_if12}\n\tni3d7: ${check_ni3d7}\n"

	# Downsample based on which clone has more telo reads
	if [[ ${check_if12} -gt ${check_ni3d7} ]]; then
	
		echo -e "\nif12 has more telo reads than ni3d7; therefore, if12 will be downsampled to equal the number of reads in ni3d7.\n"
		seqkit sample -s 2 -n ${check_ni3d7} ${out}/${ir_bn}.telo.fasta -o ${out}/${ir_bn}.sample.telo.fasta
	else
	
		echo -e "\nni3d7 has more telo reads than if12; therefore, ni3d7 will be downsampled to equal the number of reads in if12.\n"
		seqkit sample -s 5 -n ${check_if12} ${out}/${nir_bn}.telo.fasta -o ${out}/${nir_bn}.sample.telo.fasta
	fi
else

	echo -e "\nNo telo files exist. Please create telo files. Exiting.\n"
	exit 1

fi

## Collect basic stats on fasta files
ck_tFa=$( ls ${out} | egrep ".telo.*fasta" | egrep -v "lenStats" | egrep -v "slide.fa" | wc -l )
if [[ "${ck_tFa}" == 3 ]]; then
	echo -e "\nAll combined telo files exist. \n"
	ss_fa=$( ls ${out} | egrep "sample" | egrep "fasta" | egrep -v "lenStat" | egrep -v "slide.fa" )
	ss_fa_pre=$( basename -s .sample.telo.fasta ${ss_fa} )
	echo -e "\n\tThe subsampled file: ${ss_fa}\n"
	ns_fa=$( ls ${out} | egrep -v "${ss_fa_pre}" | egrep "fasta" | egrep -v "lenStat" | egrep -v "slide.fa" )
	echo -e "\n\tThe nonsampled file: ${ns_fa}\n" 
	
	bioawk -c fastx '{ print $name, length($seq) }' < ${out}/${ss_fa} | cut -f2 | r_fasta_basicStats.r > ${out}/${ss_fa}.lenStats
	bioawk -c fastx '{ print $name, length($seq) }' < ${out}/${ns_fa} | cut -f2 | r_fasta_basicStats.r > ${out}/${ns_fa}.lenStats

	ss_medl=$( cat ${out}/${ss_fa}.lenStats | egrep "median" | cut -d' ' -f 3 | awk '{print int($1+0.5)}' )
	ns_medl=$( cat ${out}/${ns_fa}.lenStats | egrep "median" | cut -d' ' -f 3 | awk '{print int($1+0.5)}' )
	echo -e "\n\tThe median read length of the subsampled file: ${ss_fa} \n\t\t${ss_medl}\n"
	echo -e "\n\tThe median read length of the nonsubsampled file: ${ns_fa} \n\t\t${ns_medl}\n"

	if (( ${ss_medl} < ${ns_medl} )); then
		ss_ns_ratio=$( echo "scale=2; ${ss_medl}/${ns_medl}" | bc )
		echo -e "\n${ss_fa} median length less than ${ns_medl}, ratio: ${ss_ns_ratio}\n"
	else
		ns_ss_ratio=$( echo "scale=2; ${ns_medl}/${ss_medl}" | bc )
		echo -e "\n${ns_fa} median length less than ${ss_fa}, ratio: ${ns_ss_ratio}\n" 
	fi

else 
	echo -e "\nAll combined telo files do not exist. Exiting. \n"
	exit 1
fi

## Determine sliding window size and step
if [[ ! -z "${ss_ns_ratio}" || ! -z "${ns_ss_ratio}" ]]; then

	echo -e "\nRatio variable found, determining window size and step.\n"
else

	echo -e "\nRatio variable empty, exiting.\n"; exit 1
fi

#win_sz=$( expr `echo -e "${ss_medl}\n${ns_medl}" | sort -nr | head -n 1` / 10 )
#win_st=$( expr ${win_sz} / 5 )
win_sz=200
win_st=100
echo -e "\nWindow size is ${win_sz} and window step is ${win_st}.\n"

## Create sliding window files
if [[ ! -z "${win_sz}" && ! -z "${win_st}" ]]; then

	echo -e "\nWindow size and window step were found, starting sliding window.\n" 
	echo -e "${out}/${ss_fa}\n${out}/${ns_fa}" > ${TMP}/slide.fofn
	parallel -k -j 2 sliding_window.py -i {} -o {}.slide.fa -w ${win_sz} -s ${win_st} :::: ${TMP}/slide.fofn
else

	echo -e "\nWindow size and window step were not found, exiting.\n"; exit 1
fi

## Create files for each read
ck_tmp_f=$( ls ${t_dir} | egrep "*.fa" | wc -l )
if (( "${ck_tmp_f}" < 50000 )); then

	echo -e "\nCheck for split files was unsuccessful, starting to split files."
	if [[ -f ${out}/${ss_fa}.slide.fa && -f ${out}/${ns_fa}.slide.fa ]]; then

		echo -e "\nBoth slide files found, starting to parse files.\n"
		ls ${out} | egrep "slide" > ${TMP}/split.fofn
		parallel -k -j ${threads} split_fasta.sh -i ${out}/{} -o ${t_dir} :::: ${TMP}/split.fofn

	else

		echo -e "\nSlide files not found, exiting. \n"; exit 1
	fi

else

	echo -e "\nCheck for split files was successful, skipping split."
fi

# Deactivate bioconda
conda deactivate

## Check for .fa files and process them into large csv
file_num=$( ls ${t_dir} | egrep "*.fa" | wc -l )
f1ln=$( cat ${out}/`head -n 1 ${TMP}/split.fofn` | egrep "^>" | wc -l )
f2ln=$( cat ${out}/`tail -n 1 ${TMP}/split.fofn` | egrep "^>" | wc -l )
combine_ln=$( expr ${f1ln} + ${f2ln} )

if [[ "${file_num}" == "${combine_ln}" ]]; then
	echo -e "\nCorrect number of files found. Starting to count telomere percents."
	echo -e "\nCount grep pattern: ${ct_grep}."
	ls ${t_dir} | egrep "*.fa" > ${TMP}/fa.fofn

	# Run telo percent count 
	echo -e "\nThis is parallel cmd:\n\tparallel -k -j ${threads} ct_telo.sh -g \"${ct_grep}\" -i ${t_dir}/{} -o ${TMP}/telomere_ranges.perc.csv :::: ${TMP}/fa.fofn" 
	parallel -k -j ${threads} ct_telo.sh -g \"${ct_grep}\" -i ${t_dir}/{} -o ${TMP}/telomere_ranges.perc.csv :::: ${TMP}/fa.fofn
else 
	echo -e "\nCorrect number of files not found. Exiting."
	exit 1
fi

# Cleaning up
head -n 1 ${TMP}/telomere_ranges.perc.csv > ${TMP}/tmp.csv
tail -n +2 ${TMP}/telomere_ranges.perc.csv | sort >> tmp.csv
mv ${TMP}/tmp.csv ${TMP}/200sw_telomere_ranges.perc.sorted.csv
rm -rf ${TMP}/tmp

echo -e "\nJOBS DONE... at `date`."
