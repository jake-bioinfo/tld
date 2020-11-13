#!/bin/bash
# Pull telomeres from irradiated and non-irradiated pfalci, sort, count length compare
# For new Sequel 2 data

# Set initial variables
s1=''
s2=''
ref=''
w_dir=''
plat=''
mot=''
bd=''
el=''
he=''

# Set help dialog
help_inf="

        Usage: telo_homer.sh -s <sample_1> -a <sample_2> -w <work_dir> -r <reference> -p <platform> -b <binpath> -t <threads>

                -h|--help               Prints help for telo_homer

                -s|--sample_1		Input sample 1 fastq, don't use underscores in name (accepts fastq or fastq.gz)             

		-a|--sample_2		Input sample 2 fastq, don't use underscores in name (accepts fastq or fastq.gz)
		
		-w|--work_dir		Working directory

                -r|--reference          Reference genome

		-b|--binpath		Path to necessary scripts
	
		-p|--platform		Takes arguements ont or pb, default is pb

		-m|--motif		Predicted motif length

		-l|--end_length		End length of chromosome to analyze for telomere motif, DEFAULT = 200 bp

		-e|--high_euk		If your organism is a higher eukaryote use this flag

                -t|--threads            integer for number of threads to use for operation, DEFAULT = max-2

                **                      -s,-a,-w,-r,-n !required!

"

# Chec for any args
if [[ -z $@ ]]; then
        echo -e "$help_inf"
        exit 1
fi

# Read options
ARGS=$( getopt -o h::s:a:w:r:b:p:m:l:e::t: -l "help::,sample_1:,sample_2:,work_dir:,reference:,\
binpath:,platform:,motif:,end_length:,high_euk::,threads:" -n "telo_homer.sh" -- "$@" );


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
                -b|--binpath)
                        shift;
                        if [[ -n $1 ]]; then
                                bd=$1;
                                shift;
                        fi
                ;;
		-p|--platform)
			shift;
			if [[ -n $1 ]]; then
				plat=$1;
				shift
			fi 
		;;
		-m|--motif)
			shift;
			if [[ -n $1 ]]; then
				mot=$1;
				shift
			fi
		;;
                -l|--end_length)
                        shift;
	                if [[ -n $1 ]]; then
				el=$1;
				shift
			fi
		;;
                -e|--high_euk)
			shift;
				he=1;
			shift
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
if [[ -z $s1 || -z $s2 || -z $ref || -z $w_dir || -z $bd || -z $mot ]]; then
        echo -e "\nRequired options (-s,-a,-w,-r,-b, -m) were not supplied, exiting\n"
        exit 1;
fi

# Check platform
if [[ $plat = "pb" ]]; then
	echo -e "\nPlatform set to PacBio, continuing...\n"

elif [[ $plat = "ont" ]]; then 
	echo -e "\nPlatform set to Oxford Nanopore, continuing...\n"

elif [[ -z $plat ]]; then
	echo -e "\nPlatform not specified, default is PacBio, continuing...\n"
        plat="pb"

fi

# Check for working directory
if [[ -d ${w_dir} ]]; then
	echo -e "\nWorking directory already created, continuing...\n"
	cd ${w_dir};

else
	echo -e "\nWorking directory not created, creating it\n"
	mkdir ${w_dir}; cd ${w_dir}

fi

# Check if threads set
if [[ -z $t ]]; then
        echo -e "\nThreads not set, setting to max-2.\n"
        t=$( echo `nproc --all`-2 | bc )
fi

# Set end length to default
if [[ -z $el ]]; then
	echo -e "\nEnd length was not set, setting to 200.\n"
	el=200
fi

# Set read out
readout_vars="

        Sample 1: $s1

        Sample 2: $s2

	Reference genome: $ref

	Bin path: $bd

	Platform: $plat

	Predicted motif length: $mot

	End length to analyze: $el

	Higher eukaryote set: $he

        Number of threads: $t

        "

# Echo vars for authentication
echo -e "\nThese are the variables:
        ${readout_vars}
        "

# Export directories and file paths
export TMP=`pwd`
export bin_path=${bd}

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
# Check if files are fastq and fasta for reference
ch_fq1=$( realpath ${s1} | egrep "*.fastq" )
ch_fq2=$( realpath ${s2} | egrep "*.fastq" )

# Check reference for fna or fasta
ch_fna=$( echo ${ref} | egrep "*.fna" );

if [[ -n ${ch_fna} ]]; then
	ch_ref=$( realpath ${ref} | egrep "*.fna" )
else
	ch_ref=$( realpath ${ref} | egrep "*.fasta" )
fi

if [[ -n ${ch_fq1} ]]; then 
	echo -e "\nFastq for sample 1 found, continuing."
	s1floc=$( realpath $( dirname ${s1} ))
	s1f_bn=$( basename -s .fastq ${s1} );
else 
	echo -e "\nFastq for sample 1 not found or is not fastq, exiting."
	exit 2
fi

if [[ -n ${ch_fq2} ]]; then
	echo -e "\nFastq for sample 2 found, continuing."
	s2floc=$( realpath $( dirname ${s2} ))
	s2f_bn=$( basename -s .fastq ${s2} );
else
        echo -e "\nFastq for sample 2 not found or is not fastq, exiting."
        exit 2
fi


if [[ -n ${ch_ref} ]]; then
        echo -e "\nFasta for reference found, continuing."
        reffloc=$( realpath $( dirname ${ref} ))
        
	if [[ -n ${ch_fna} ]]; then
		reff_bn=$( basename -s .fna ${ref} );
		reff_bn_s="fna"
	else
		reff_bn=$( basename -s .fasta ${ref} );
		reff_bn_s="fasta"
	fi
else
        echo -e "\nFasta for reference not found or is not fasta, exiting."
        exit 2
fi

# Check if higher eukaryotes set, if so reduce siz
# of reference for future analysis
if [[ -z ${he} ]]; then
	echo -e "\nHigher eukaryote option not set, continuing."
else
	echo -e "\nHigher eukaryote option set, reducing size of reference for future analysis"
	seqkit concat <(seqkit subseq -r 1:1000000 ${ref}) <(seqkit subseq -r -1000000:-1 ${ref}) \
		> ${in}/mod_${reff_bn}.${reff_bn_s}
	export ref=${in}/mod_${reff_bn}.${reff_bn_s}
fi

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

# Create ref index if not already there
check_ref_i=$( ls $( echo ${reffloc} ) | egrep "${reff_bn}.fasta.fai" );
[ -f ${check_ref_i} ] && \
{ echo -e "\nReference not present, creating it."; samtools faidx ${ref}; } || \
{ echo -e "\nReference present, proceeding."; }

# Create telo bed file
echo -e "\nProcessing index file to create telo.bed.\n"
touch ${out}/telo.bed

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
	echo -e "${lab_p_5}\t${chr_head}\t0\t${el}\t0" >> ${out}/telo.bed
	# Print 5' - lines to telo.bed file
	echo -e "${lab_n_5}\t${chr_head}\t0\t${el}\t1" >> ${out}/telo.bed

	# Create 3' chromosome locations
	end_loc=$( echo -e "${l}" | cut -f2 );
	st_loc=$( echo "scale=1; ${end_loc} - ${el}" | bc )	

	# Print 3' + lines to telo.bed file
	echo -e "${lab_p_3}\t${chr_head}\t${st_loc}\t${end_loc}\t0" >> ${out}/telo.bed
	# Print 3' - lines to telo.bed file
	echo -e "${lab_n_3}\t${chr_head}\t${st_loc}\t${end_loc}\t1" >> ${out}/telo.bed
	
	# Increase iterator
	let i++

done < ${ref}.fai

# Check if homer has been run, if not run it. 
check_homer=$( ls ${out} | egrep "homerMotifs.all.motifs" );
[ -f ${check_homer} ] && \
{ echo -e "\nHomer result file not found, run homer.";
	mot2=$( echo "2*${mot}" | bc ); 

	# Run Homer on reference
	findMotifsGenome.pl ${out}/telo.bed ${ref} \
	${out} -len ${mot} -noweight -S ${mot2} -maxN 2 -preparse -p ${threads}; } || \
{ echo -e "\nHomer result file found, proceeding."; }

# Parse homer output
# Create temp.motifs file with the first 5 motifs by multiplicity
# Check if motifs file already created
check_motif=$( ls ${out} | egrep "temp.motifs" );
[ -f ${check_motif} ] && \
{ echo -e "\nMotifs file not found, creating them."; 
	egrep "^>" ${out}/homerMotifs.all.motifs | egrep "CC|GG" | sed "s:,:\t:g" \
	| sed "s:>::g" | sort -k14 -r | cut -f1 | head -n 3 > ${out}/temp.motifs

	# Create reverse and reverse complements of all motifs
	touch ${out}/temp_rc.motifs
	while IFS= read -r l; do 
		rev_comp=$( echo "${l}" | rev | tr "[ATGCatgc]" "[TACGtacg]" );
		rev=$( echo "${l}" | rev );
		echo -e "${l}\t${rev_comp}\t${rev}" >> ${out}/temp_rc.motifs;
		
	done < ${out}/temp.motifs; } || \

{ echo -e "\nMotifs file found, proceeding."; }

# Create initial grep pull down statement
# Check if higher eukaryotes option was set, if yes add more grep expressions
[ -z ${he} ] && \
{ echo -e "\nHigher eukaryotes option not set, proceeding for lower eukaryotes\n"
	while IFS= read -r l; do
		orig=$( echo "${l}" | cut -f1 )
		rev_comp=$( echo "${l}" | cut -f2 )
		rev=$( echo "${l}" | cut -f3 )
		grep_exp+="${orig}${orig}|${rev}${rev}|${rev_comp}${rev_comp}|"
		ct_grep+="${orig}|${rev}|${rev_comp}|" 
	done < ${out}/temp_rc.motifs; } || \

{ echo -e "\nHigher eukarytoes option set, proceeding for higher eukaryotes\n"
	while IFS= read -r l; do
		orig=$( echo "${l}" | cut -f1 )
		rev_comp=$( echo "${l}" | cut -f2 )
		rev=$( echo "${l}" | cut -f3 )
		#grep_exp+="${orig}${orig}${orig}${orig}${orig}${orig}|${rev}${rev}${rev}${rev}${rev}${rev}|${rev_comp}${rev_comp}${rev_comp}${rev_comp}${rev_comp}${rev_comp}|"
		grep_exp+="${orig}${orig}${orig}${orig}|${rev}${rev}${rev}${rev}|${rev_comp}${rev_comp}${rev_comp}${rev_comp}|"
		grep_exp=$( echo ${grep_exp} | sed "s: ::g" );
		ct_grep+="${orig}|${rev}|${rev_comp}|" 
	done < ${out}/temp_rc.motifs; }

# check grep expression
grep_exp=$( echo "${grep_exp}" | head -c -2 );
echo -e "\nThis is grep expression: ${grep_exp}\n";
ct_grep=$( echo "${ct_grep}" | sed "s:.$::g" );
echo -e "\nThis is count grep expression: ${ct_grep}\n";

# Downsample fastq files based on GB
# Check if downsampled file already exists
if [[ -f ${s1}.s.fastq || -f ${s2}.s.fastq ]]; then
	echo -e "\nDownsampled file already exists, continuing..."

else

	# Calculate GB
	ch_s1gb=$( grep -v "^@" ${s1} | wc -c )
	ch_s2gb=$( grep -v "^@" ${s2} | wc -c )

	echo -e "\nBases in each file\n\t${s1}:${ch_s1gb}\n\t${s2}:${ch_s2gb}"

	if [[ ${ch_s1gb} -gt ${ch_s2gb} ]]; then
	
		# Downsampling file 1 to match file 2
		echo -e "\nDownsampling ${s1} to match GB in ${s2}"
		r=$( echo "${ch_s2gb}/${ch_s1gb}" | bc -l )
		seqkit sample -s 8 -j ${threads} -p ${r} ${s1} -o ${in}/${s1f_bn}.s.fastq
		s1=${in}/${s1f_bn}.s.fastq

	else

		# Downsampling file 2 to match file 1
		echo -e "\nDownsampling ${s2} to match GB in ${s1}"
        	r=$( echo "${ch_s1gb}/${ch_s2gb}" | bc -l )
        	seqkit sample -s 8 -j ${threads} -p ${r} ${s2} -o ${in}/${s2f_bn}.s.fastq
        	s2=${in}/${s2f_bn}.s.fastq

	fi

	# Recheck GB
	ch_s1gb=$( grep -v "^@" ${s1} | wc -c )
	ch_s2gb=$( grep -v "^@" ${s2} | wc -c )

	echo -e "\nBases in each file after downsample\n\t${s1}:${ch_s1gb}\n\t${s2}:${ch_s2gb}"

fi

# Check for telomere reads combined file
check_telo=$( ls ${out} | egrep "telo.fastq" | wc -l );
if [[ "${check_telo}" == 2 ]]; then
	echo -e "\nTelo files exist, not creating them.\n"
else
	echo -e "\nTelo files do not exist, creating them. \n"

	# Grep ends of reads for telomere repeats
	seqkit grep -R 1:1000 -j ${threads} -s -i -r -p "${grep_exp}" ${s1}> ${out}/${s1f_bn}.telo.fastq
	seqkit grep -R -1000:-1 -j ${threads} -s -i -r -p "${grep_exp}" ${s1}>> ${out}/${s1f_bn}.telo.fastq
	
	# Remove duplicates
	seqkit rmdup -j ${threads} ${out}/${s1f_bn}.telo.fastq > ${out}/${s1f_bn}.t.telo.fastq
	mv ${out}/${s1f_bn}.t.telo.fastq ${out}/${s1f_bn}.telo.fastq

	# Grep ends of reads for telomere repeats
	seqkit grep -R 1:1000 -j ${threads} -s -i -r -p "${grep_exp}" ${s2} > ${out}/${s2f_bn}.telo.fastq
	seqkit grep -R -1000:-1 -j ${threads} -s -i -r -p "${grep_exp}" ${s2}  >> ${out}/${s2f_bn}.telo.fastq
	
	# Remove duplicates
	seqkit rmdup -j ${threads} ${out}/${s2f_bn}.telo.fastq > ${out}/${s2f_bn}.t.telo.fastq
	mv ${out}/${s2f_bn}.t.telo.fastq ${out}/${s2f_bn}.telo.fastq

	# Converting to fasta
	seqkit fq2fa -j ${threads} ${out}/${s1f_bn}.telo.fastq -o ${out}/${s1f_bn}.telo.fasta
	seqkit fq2fa -j ${threads} ${out}/${s2f_bn}.telo.fastq -o ${out}/${s2f_bn}.telo.fasta

fi

# Check number of reads in each combined telo file and random downsample larger file
if [[ -f ${out}/${s1f_bn}.telo.fasta && -f ${out}/${s2f_bn}.telo.fasta ]]; then
	echo -e "\nBoth combined telo files exist. Determining number of reads in each. \n"
	
	check_s2=$( egrep "^>" ${out}/${s2f_bn}.telo.fasta | wc -l )
	check_s1=$( egrep "^>" ${out}/${s1f_bn}.telo.fasta | wc -l )
	echo -e "\nNumber of reads in each telo follows.\n\tif12: ${check_s2}\n\tni3d7: ${check_s1}\n"

	# Downsample based on which clone has more telo reads
	if [[ ${check_s2} -gt ${check_s1} ]]; then
	
		echo -e "\nS2 has more telo reads than s1; therefore, S2 will be downsampled to equal the number of reads in S1.\n"
		seqkit sample -s 8 -n ${check_s1} ${out}/${s2f_bn}.telo.fasta -o ${out}/${s2f_bn}.sample.telo.fasta
	else
	
		echo -e "\nS1 has more telo reads than S2; therefore, S1 will be downsampled to equal the number of reads in S2.\n"
		seqkit sample -s 8 -n ${check_s2} ${out}/${s1f_bn}.telo.fasta -o ${out}/${s1f_bn}.sample.telo.fasta
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
	
	bioawk -c fastx '{ print $name, length($seq) }' < ${out}/${ss_fa} | cut -f2 | ${bin_path}/r_fasta_basicStats.r > ${out}/${ss_fa}.lenStats
	bioawk -c fastx '{ print $name, length($seq) }' < ${out}/${ns_fa} | cut -f2 | ${bin_path}/r_fasta_basicStats.r > ${out}/${ns_fa}.lenStats

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

	# If ont, replace hyphens in read name
	if [[ $plat = "ont" ]]; then
		sed -i "s: .*::g" ${out}/${ss_fa}
		sed -i "s: .*::g" ${out}/${ns_fa}
		sed -i "s:-:_:g" ${out}/${ss_fa}
		sed -i "s:-:_:g" ${out}/${ns_fa}
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

win_sz=200
win_st=100
echo -e "\nWindow size is ${win_sz} and window step is ${win_st}.\n"

# Maintain input order of files
ss_bn=$( basename -s .sample.telo.fasta ${ss_fa} );
ns_bn=$( basename -s .telo.fasta ${ns_fa} );
s1_bn=$( basename -s .fastq ${s1f_bn} );
s2_bn=$( basename -s .fastq ${s2f_bn} );

if [[ "${ss_bn}" == "${s1_bn}" ]]; then
	export sample1=${ss_bn}.sample.telo.fasta
	export sample2=${ns_bn}.telo.fasta

else
	export sample1=${ns_bn}.telo.fasta
	export sample2=${ss_bn}.sample.telo.fasta

fi

## Create sliding window files
if [[ ! -z "${win_sz}" && ! -z "${win_st}" ]]; then

	echo -e "\nWindow size and window step were found, starting sliding window.\n" 
	#echo -e "${out}/${ss_fa}\n${out}/${ns_fa}" > ${TMP}/slide.fofn
	echo -e "${out}/${sample1}\n${out}/${sample2}" > ${TMP}/slide.fofn
	parallel -k -j 2 ${bin_path}/sliding_window.py -i {} -o {}.slide.fa -w ${win_sz} -s ${win_st} :::: ${TMP}/slide.fofn
else

	echo -e "\nWindow size and window step were not found, exiting.\n"; exit 1
fi

## Create files for each read
ck_tmp_f=$( ls ${t_dir} | egrep "*.fa" | wc -l )
if (( "${ck_tmp_f}" < 50000 )); then

	echo -e "\nCheck for split files was unsuccessful, starting to split files."
	if [[ -f ${out}/${ss_fa}.slide.fa && -f ${out}/${ns_fa}.slide.fa ]]; then

		echo -e "\nBoth slide files found, starting to parse files.\n"
		#ls ${out} | egrep "slide" > ${TMP}/split.fofn
		echo -e "${sample1}.slide.fa\n${sample2}.slide.fa" > ${TMP}/split.fofn
		parallel -k -j ${threads} ${bin_path}/split_fasta.sh -i ${out}/{} -o ${t_dir} :::: ${TMP}/split.fofn

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
	echo -e "\nCount grep pattern: ${ct_grep}."
	ls ${t_dir} | egrep "*.fa" > ${TMP}/fa.fofn

	# Run telo percent count 
	echo -e "\nThis is parallel cmd:\n\tparallel -k -j ${threads} ct_telo.sh -g \"${ct_grep}\" -i ${t_dir}/{} -o ${TMP}/telomere_ranges.perc.csv :::: ${TMP}/fa.fofn" 
	parallel -k -j ${threads} ${bin_path}/ct_telo.sh -g \"${ct_grep}\" -i ${t_dir}/{} -o ${TMP}/telomere_ranges.perc.csv :::: ${TMP}/fa.fofn
else 
	echo -e "\nCorrect number of files not found. Exiting."
	exit 1
fi

# Cleaning up
sm1=$( echo ${sample1} | cut -d'.' -f1 );
sm2=$( echo ${sample2} | cut -d'.' -f1 );

head -n 1 ${TMP}/telomere_ranges.perc.csv > ${TMP}/tmp.csv
tail -n +2 ${TMP}/telomere_ranges.perc.csv | egrep ${sm1} | sort >> tmp.csv
tail -n +2 ${TMP}/telomere_ranges.perc.csv | egrep ${sm2} | sort >> tmp.csv
mv ${TMP}/tmp.csv ${TMP}/200sw_telomere_ranges.perc.sorted.csv
rm -rf ${TMP}/tmp

echo -e "\nJOBS DONE... at `date`."
