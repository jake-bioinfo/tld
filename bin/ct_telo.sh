#!/bin/bash
# Reads fasta and creates csv file with read name, range of sequence in read analyzed, % of range 
# which contains telomere repeats.

# Get opts variables initialization
input_f=''
out_f=''

# Set help info
help_readout="

        Usage: ct_telo.sh -i <input_single_read.fa> -o <output_csv_file>

                -h      --help          Displays this message

                -i      --input         Input file, full path or local path to single fasta file

                -o      --output_file   Output file, if out file already exists append instead of write

		-g	--grep_exp	Grep expression for caculating telomere percentage

                *       Long form of options requires \"=\"
                        ex: ct_telo.sh -i what.fa --output_file=/path/to/out.csv

                **      -i, -o          Required arguments

        "

# Check for file
if [[ -z $@ ]]; then
        echo -e "$help_readout"
        exit 1
fi

# GetOpt
ARGS=$( getopt -o h::i:o:g: -l "help::,input_file:,output_file:,grep_exp:" -n "ct_telo.sh" -- "$@" );

eval set -- "$ARGS";

# extract options and their arguments into variables
while true; do
        case "$1" in
                -h|--help)
                shift;
                        echo -e "${help_readout}";
                exit 1;
                ;;
                -i|--input_file)
                shift;
                if [[ -n $1 ]]; then
                        input_f=$1;
                shift;
                fi
                ;;
                -o|--output_file)
                shift;
                if [[ -n $1 ]]; then
                        out_f=$1;
                shift;
                fi
                ;;
		-g|--grep_exp)
		shift;
		if [[ -n $1 ]]; then
			grep_x=$1;
		shift;
		fi
		;;
                --)
                shift;
                break;
                ;;
        esac
done

# Check req args
if [[ -z $input_f ]]; then
        echo -e "

                -i is required and was not supplied, exiting

                "
        exit 1;
fi

# Check if out file does not exist, create it if it doesn't
if [[ ! -f ${out_f} ]]; then

	echo -e "${out_f} does not exist, creating it.\n";
	echo -e "Sample.name,Read.name,Start.range,End.range,Seq.ln,Telomere.count,Percent.telomere.repeat" > ${out_f};

fi

# Perform telomere count
while read line
do
    if [[ ${line:0:1} == '>' ]]; then
	# Get read information and write
	sample_nm=$( echo -e "${line}" | awk -F'__' '{print $1}' )
        read_nm=$( echo -e "${line#>}" | sed -e "s:/::g" | sed -e "s: :_:g" \
		| sed -e "s:|:_:g" | sed -e "s/:/-/g" | awk -F'__' '{print $2}' )
	rng_st=$( echo "${read_nm}" | cut -d'-' -f 2 )
	rng_end=$( echo ${read_nm} | cut -d'-' -f 3 | sed "s:.fa::" ) 
	out_ln1=$( echo "${sample_nm},${read_nm},${rng_st},${rng_end}" )

    elif [[ ${line:0:4} != 'Read' ]]; then
	# Count sequence length, telomere repeats in sequence, telomere percent
        seq=$( echo ${line} )
	seq_ln=$( echo ${seq} | egrep -o "." | tr -d '\n' | wc -c )
	telo_ct=$( echo ${seq} | egrep -o "${grep_x}" \
	| tr -d '\n' | wc -c )
	telo_perc=$( echo "scale=2; (${telo_ct}*100)/${seq_ln}" | bc )
	out_ln2=$( echo "${seq_ln},${telo_ct},${telo_perc}" )
    fi

done < ${input_f}

# Append to end of the file
echo -e "${out_ln1},${out_ln2}" >> ${out_f}
