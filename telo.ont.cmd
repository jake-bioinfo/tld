./telo_pipe.sh -w ~/pfalci/201908_telomere_lengths/data/202011_ont_processing \
	-a ~/pfalci/shared_data/2020_fastq/20200208_combined_minion_3d7_wt.fastq \
	-f ~/pfalci/shared_data/2020_fastq/rad51irrB11.merge.fastq \
	-r ~/pfalci/shared_data/refGenomes/ncbi-genomes-2020-01-14/GCA_000002765.3_GCA_000002765_genomic_chr_lab.fna \
	-d ~/pfalci/201908_telomere_lengths/results/202011_ont \
	-p "202011" -s "ont_wt,ont_r51irr" -m "2572,1136" -n ont -t 6 
