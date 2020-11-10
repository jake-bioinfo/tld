#!/bin/bash
# Pull telomeres from irradiated and non-irradiated pfalci, sort, count length compare
# For new Sequel 2 data

# Activate r_env 
source ~/software/anaconda3/etc/profile.d/conda.sh
conda activate r_env

# Export directories
ni_raw=/home/har2011/pfalci/sync/202005_seq2_pfalci/m54333U_200503_055320.Q20.fastq
i_raw=/home/har2011/pfalci/sync/202005_seq2_pfalci/m54333U_200508_174519.Q20.fastq
export TMP=`pwd`

# Create directories
check_out_dir=$( cd ${TMP} | ls | egrep "output" );
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

# Pull telomere reads from each file

# Check for telomere reads combined file
check_telo=$( ls ${out} | egrep "telo.fastq" | wc -l );
if [[ "${check_telo}" == 2 ]]; then
	echo -e "\nTelo files exist, not creating them.\n"
else
	echo -e "\nTelo files do not exist, creating them. \n"
	cat ${ni_raw} | egrep -i -B1 -A2 "TTCAGGGTTCAGGG|TTTAGGGTTTAGGG|CCCTAAACCCTAAA|CCCTGAACCCTGAA|TTCAGGGTTTAGGG|CCCTAAACCCTGAA" \
		| sed ":--:d" > ${out}/3d7_2s.telo.fastq
	seqtk seq -a ${out}/3d7_2s.telo.fastq > ${out}/3d7_2s.telo.fasta

	cat ${i_raw} | egrep -i -B1 -A2 "TTCAGGGTTCAGGG|TTTAGGGTTTAGGG|CCCTAAACCCTAAA|CCCTGAACCCTGAA|TTCAGGGTTTAGGG|CCCTAAACCCTGAA" \
		| sed ":--:d" > ${out}/irr_2s.telo.fastq
	seqtk seq -a ${out}/irr_2s.telo.fastq > ${out}/irr_2s.telo.fasta

fi

# Check number of reads in each combined telo file and random downsample larger file, should be irrf12 bc it had a coverage +30X

if [[ -f ${out}/3d7_2s.telo.fasta && -f ${out}/irr_2s.telo.fasta ]]; then
	echo -e "\nBoth combined telo files exist. Determining number of reads in each. \n"
	
	check_if12=$( egrep "^>" ${out}/irr_2s.telo.fasta | wc -l )
	check_ni3d7=$( egrep "^>" ${out}/3d7_2s.telo.fasta | wc -l )
	echo -e "\nNumber of reads in each telo follows.\n\tif12: ${check_if12}\n\tni3d7: ${check_ni3d7}\n"

	# Downsample based on which clone has more telo reads
	if [[ ${check_if12} -gt ${check_ni3d7} ]]; then
	
		echo -e "\nif12 has more telo reads than ni3d7; therefore, if12 will be downsampled to equal the number of reads in ni3d7.\n"
		seqtk sample -s 100 ${out}/irr_2s.telo.fasta \
		${check_ni3d7} ${check_ni3d7} > ${out}/irr_2s.sample.telo.fasta
	else
	
		echo -e "\nni3d7 has more telo reads than if12; therefore, ni3d7 will be downsampled to equal the number of reads in if12.\n"
		seqtk sample -s 100 ${out}/3d7_2s.telo.fasta \
		${check_if12} ${check_if12} > ${out}/3d7_2s.sample.telo.fasta
	fi
else

	echo -e "\nNo telo files exist. Please create telo files. Exiting.\n"
	exit 1

fi

# Sliding window of telomere reads

## Collect basic stats on fasta files
ck_tFa=$( ls ${out} | egrep ".telo.*fasta" | egrep -v "lenStats" | egrep -v "slide.fa" | wc -l )
if [[ "${ck_tFa}" == 3 ]]; then
	echo -e "\nAll combined telo files exist. \n"
	ss_fa=$( ls ${out} | egrep "sample" | egrep "fasta" | egrep -v "lenStat" | egrep -v "slide.fa" )
	ss_fa_prefix=$( echo ${ss_fa} | egrep -o ".......sample" | cut -c -3 )
	echo -e "\n\tThe subsampled file: ${ss_fa}\n"
	ns_fa=$( ls ${out} | egrep -v "${ss_fa_prefix}" | egrep "fasta" | egrep -v "lenStat" | egrep -v "slide.fa" )
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

# Deactivate r_env, activate bioconda
conda deactivate
conda activate ddocent_env

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
	parallel -k -j 4 sliding_window.py -i {} -o {}.slide.fa -w ${win_sz} -s ${win_st} :::: ${TMP}/slide.fofn
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
		parallel -k -j 8 split_fasta.sh -i ${out}/{} -o ${t_dir} :::: ${TMP}/split.fofn

	else

		echo -e "\nSlide files not found, exiting. \n"; exit 1
	fi

else

	echo -e "\nCheck for split files was successful, skipping split."
fi

## Check for .fa files and process them into large csv
file_num=$( ls ${t_dir} | egrep "*.fa" | wc -l )
f1ln=$( cat ${out}/`head -n 1 ${TMP}/split.fofn` | egrep "^>" | wc -l )
f2ln=$( cat ${out}/`tail -n 1 ${TMP}/split.fofn` | egrep "^>" | wc -l )
combine_ln=$( expr ${f1ln} + ${f2ln} )

if [[ "${file_num}" == "${combine_ln}" ]]; then
	echo -e "\nCorrect number of files found. Starting to count telomere percents."
	ls ${t_dir} | egrep "*.fa" > ${TMP}/fa.fofn
	parallel -k -j 64 ct_telo.sh -i ${t_dir}/{} -o ${TMP}/telomere_ranges.perc.csv :::: ${TMP}/fa.fofn
else 
	echo -e "\nCorrect number of files not found. Exiting."
	exit 1
fi

head -n 1 ${TMP}/telomere_ranges.perc.csv > ${TMP}/tmp.csv
tail -n +2 ${TMP}/telomere_ranges.perc.csv | sort >> tmp.csv
mv ${TMP}/tmp.csv ${TMP}/200sw_telomere_ranges.perc.sorted.csv

echo -e "\nJOBS DONE... at `date`."
conda deactivate
