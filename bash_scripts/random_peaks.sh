#!/bin/bash

####### Yolanda Guillen ######
######### April 2020 #########

##### Estimate the significance of matching TF peak regions ####

REFGENOME="/Users/yguillen/Desktop/temp/beta_catenin_project/db_files_ciber"

# Input arguments:
# first argument = external bed to be compared
# second argument = our set of peaks
# third argument = number of iterations. Recommended: lower than 1000
# fourth argument = Reference genome (hg19 or hg38 for human, mm9 or mm10 for mouse)
# fifth argument = Label for the output results



if [[ $# -ne 5 ]]; then
	echo $'\n\nYOU DID NOT RUN IT CORRECTLY\n\nRemember to input two bed files, number of iterations, genome, output_name and run. Example:\n\n bash random_peaks.sh set_to_compare.bed your_set.bed num_iterations hg38/mm10/hg19/mm9 comparison_name\n\nCheck the genome version, the coordinates of both beds need to be based on the same version\n'
	
	echo $'Try again\n\nYolanda\n\n'
	exit 1
else

	echo $'Remember to check the genome version of the bed files\n'
targetbed=$1
bcatpeaks=$2
it=$3
genome=$4
outname=$5

echo $'\nCheck your reference genome, human or mouse, and the version\n'

# Create as many random peaks as peaks in the set to compare, of wid length (length distribution external peaks), 100 iteration

	echo $'Intersecting both sets of peaks...\n'
	# order coordinates input beds and choose unique peaks
	sort -u -k4,4  ${targetbed} | sort -k1,1 -k2,2n | awk '{print$1"\t"$2"\t"$3"\t"$4}' > sorted_input.bed
	sort -u -k4,4 ${bcatpeaks} | sort -k1,1 -k2,2n | awk '{print$1"\t"$2"\t"$3"\t"$4}' > sorted_bcat.bed

	# Range of peak length in external bed
	widmin=$(awk '{print$3-$2}' sorted_input.bed | sort -n | head -n 1)
	widmax=$(awk '{print$3-$2}' sorted_input.bed | sort -rn | head -n 1)

	widmin2=$(awk '{print$3-$2}' sorted_bcat.bed | sort -n | head -n 1)
	widmax2=$(awk '{print$3-$2}' sorted_bcat.bed | sort -rn | head -n 1)
	
	echo $'Range of peak length in external bed (first bed)--> min: '$widmin' max: '$widmax''
	echo $'Range of peak length in targe bed (second bed)--> min: '$widmin2' max: '$widmax2''

	#num=$(wc -l sorted_input.bed)
	#alnum=$(wc -l sorted_bcat.bed)
	
	num2=$(wc -l sorted_input.bed | awk '{print$1}')
	alnum2=$(wc -l sorted_bcat.bed | awk '{print$1}')
	echo $'Number of unique peaks from external bed -->\t'${num2}''
	echo $'Number of unique peaks from your bed -->\t'${alnum2}''
	echo $'\n'

	# intersect to check how many peaks are in common
	bedtools closest -a sorted_input.bed -b sorted_bcat.bed -d | awk '($NF<=100 && $NF>=0) {print}' > ${outname}_common_sets_d100.bed
	bedtools closest -a sorted_input.bed -b sorted_bcat.bed -d | awk '($NF<=1 && $NF>=0) {print}' > ${outname}_common_sets_d1.bed

	dhun=$(wc -l ${outname}_common_sets_d100.bed)
	don=$(wc -l ${outname}_common_sets_d1.bed)
	echo $'Number of common peaks at a distance 100 bp -->\t'${dhun}''
	echo $'Number of common peaks overlapping minimum one base -->\t'${don}''

#set empty freq files
if [ -f "freq_common_d*txt" ]; then
	echo "remove existing previous freqs..."
	rm freq_common_d*txt
	echo $'\nComputing '${num2}' random peaks '${it}' times...'
else
	echo $'\nComputing '${num2}' random peaks '${it}' times...'
fi



for i in $(seq 1 ${it})
do	
	# Create random (num) random peaks, of (wid) length, i times and add number of iteration in peak name
	# hg38_chrom_filter.size or mm10_chrom_filter.size
	
	##bedtools random -l ${widran} -n ${num2} -g ${REFGENOME}/${genome}_chrom_filter.size | sort -k1,1 -k2,2n > random_peaks_${i}.bed
	bedtools shuffle -i sorted_input.bed -g ${REFGENOME}/${genome}_chrom.size -incl ${REFGENOME}/merge_oregano_${genome}.bed -noOverlapping | sort -k1,1 -k2,2n > random_peaks_${i}.bed
	awk -v ran="$i" '{print$1"\t"$2"\t"$3"\t""random_"ran"_"$4}' random_peaks_${i}.bed > random_peaks_${i}_name.bed
	rm random_peaks_${i}.bed

	# intersect and extract common
	#Commons at 100 bp
	bedtools closest -a sorted_bcat.bed -b random_peaks_${i}_name.bed -d | awk '($NF<=100 && $NF>=0) {print}' > common_random_${i}_d100.bed
	wc -l common_random_${i}_d100.bed | awk '{print$1"\t"$2}' >> freq_common_d100.txt
	rm common_random_${i}_d100.bed

	# commons at 1 bp
	bedtools closest -a sorted_bcat.bed -b random_peaks_${i}_name.bed -d | awk '($NF<=1 && $NF>=0) {print}' > common_random_${i}_d1.bed
	wc -l common_random_${i}_d1.bed | awk '{print$1"\t"$2}' >> freq_common_d1.txt
	rm common_random_${i}_d1.bed

	rm random_peaks_${i}_name.bed

done

sed -i '1i intersect_100bp	iter_100' freq_common_d100.txt
sed -i '1i intersect_1bp	iter_1' freq_common_d1.txt

paste freq_common* > ${outname}_dist_iter_peaks.txt
rm freq_common*


### alternative measure bedtools reldist

#generate random peak list to compare with our dataset


##bedtools random -l ${widran} -n ${num2} -g ${REFGENOME}/${genome}_chrom_filter.size | sort -k1,1 -k2,2n > random_test_peaks.bed
bedtools shuffle -i sorted_input.bed -g ${REFGENOME}/${genome}_chrom.size -incl ${REFGENOME}/merge_oregano_${genome}.bed -noOverlapping | sort -k1,1 -k2,2n > random_test_peaks.bed

bedtools reldist -a sorted_bcat.bed -b sorted_input.bed | awk '{print"observed""\t"$1"\t"$2"\t"$3"\t"$4}' > frequency_relative_distance_obser.txt
bedtools reldist -a sorted_bcat.bed -b random_test_peaks.bed | awk '{print"expected""\t"$1"\t"$2"\t"$3"\t"$4}' | awk NR\>1 > frequency_relative_distance_expect.txt

cat frequency_relative*obser* frequency_relative_*expect* | awk NR\>1 > ${outname}_plot_distance_dist.txt
sed -i '1i data\treldist\tcount\ttotal\tfraction' ${outname}_plot_distance_dist.txt

rm frequency_relative_distance*
rm random_test_peaks.bed


rm sorted_input.bed
rm sorted_bcat.bed

fi

echo $'DONE!\n\nYolanda\n\n'
