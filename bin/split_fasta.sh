#!/bin/bash
# split_fasta.sh
# splits long fasta file with multiple sequences into files of each sequence

# Get opts variables initialization
input_f=''
out_d=''

# Set help info
help_readout="
	
	Usage: split_fasta.sh -i <input_multi_fasta> -o <output_directory>
	
		-h	--help		Displays this message

		-i	--input		Input file, full path or local path

		-o	--output_dir	Output directory, full path or local path
		
		*	Long form of options requires \"=\" 
			ex: split_fasta.sh -i what.fa --output_dir=/path/to/out
		
		**	-i, -o		Required arguments

	"

# Check for file
if [[ -z $@ ]]; then
	echo -e "$help_readout"
	exit 1
fi

# GetOpt
ARGS=$( getopt -o h::i:o: -l "help::,input_file:,output_dir:" -n "split_fasta.sh" -- "$@" );

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
		-o|--output_dir)
                shift;
                if [[ -n $1 ]]; then
                        out_d=$1;
                shift;
                fi
                ;;
		--)
                shift;
                break;
                ;;
        esac
done
 
echo -e "\nPast optget"

# Check req args
if [[ -z $input_f || -z $out_d ]]; then
	echo -e "

		-i, -o are required arguments, exiting

		"
	exit 1;
fi

echo -e "\nInput file is: $input_f

Output dir is: $out_d

"

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
	out_f=$( echo $outfile | sed -e "s:/::g" | sed -e "s: :_:g" | sed -e "s:|:_:g" | sed -e "s/:/-/g" )
#	echo -e "\nProcessing file: $out_f"
        echo $line > ${out_d}/"$out_f"
    else
        echo $line >> ${out_d}/"$out_f"
    fi
done < ${input_f}
