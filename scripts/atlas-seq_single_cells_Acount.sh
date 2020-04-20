#!/bin/bash

####### This is a stable version. Do not modify without changing the name
# should be run from atlas-clustering.sh result folder

####### Dependencies
# bedtools		2.25.0			https://github.com/arq5x/bedtools2
# gnu grep/awk

####### Release history
# 1.0 release notes:
#  	- compatible with atlas-clustering.sh script v2.2
#	- use ...insertions.full.bed files
#	- takes into account the density of reads to identify the best cluster within a supercluster, but also the % of bad reads (spurious+irregular) within a cluster

####### load preferences and define global settings
script_name="atlas-single-cell"
script_version='1.0'

start_time=`date +%s`
day=$(date +"[%d-%m-%Y] [%T]")

# Load location of picard tools, reference genome, and result folder from external ".atlas.conf" parameter file
configuration_file="~/.atlas.conf"
if [ -f "$configuration_file" ];
	then
		echo -e "\nMissing '.atlas.conf' file in home directory.\n";
		exit 1
fi
while read line
do
    eval $( awk '$1!~/^#/ {print $0}' )
done < ~/.atlas.conf

# Store usage explanations
###### to REWRITE
USAGE="\
${script_name} v${script_version}:\tL1 peak calling and annotation for single-cell experiments. \n\
usage:\t$( basename $0 ) [options] -a annotation_list_file.txt input_file_to_annotate.bed \n\
options:\n\
\t-h Print this help menu. \n\
\t-v What version of $script_name are you using? \n\
\t-f <5' ATLAS.cluster.true.bed file> to cross with a 5' ATLAS-seq experiment and identify potential full length L1. \n\
remark: the annotation list and the annotation .bed files should be located in the same folder.\n\
"

# define experiment-specific variables (can be changed by user)
sample_list="E6C1 E6C3 E6T"
# sample_list="E6T"
exp_type="3atlas 5atlas"
# exp_type="3atlas"
atlas_ver="2.2forktest_nosoft"
window=10000 # window to merge call L1 insertion
annotation_dir="${bioinfo_root}/annotations/atlas-seq_annot"
annot5="170530_annotations_sc_5atlas.imac.txt"
annot3="170530_annotations_sc_3atlas.imac.txt"

# define global variables
file_list=""
sample_list_comma=""
window_name=$(( $window / 1000 ))
header_cluster=$( echo -e "#CHR\tSTART\tEND\tCLUSTER_ID\tNB_NRR\tCLUSTER_STRAND" )
header_insertion=$( echo -e "#CHR\tINS_START\tINS_END\tCLUSTER_ID\tNB_NRR\tINS_STRAND" )

# create result directory
current_dir=$( pwd )
result_dir="L1call_$( date +"%y%m%d_%H%M%S" )"
mkdir ${result_dir}
cd ${result_dir}

# load samples to process in array
IFS=' ' read -a sample_name <<< $sample_list
nb_samples=${#sample_name[@]}

####### L1 peak calling

for experiment in ${exp_type[@]}
do
	for sample in ${sample_name[@]}
	do
		echo -e "*** Calling ${experiment} peaks for ${sample}"

			# calculate read density and % of (irregular+spurious) reads for each cluster

			sort -k1,1 -k2,2n "../${sample}.${experiment}.v${atlas_ver}.clusters.full.bed" \
			| awk '
				($1~/^#/){
					OFS="\t" ;
					$NF=$NF ;
					print $1,$2,$3,$4,$5,$6,"GOOD_CLUSTER_DENSITY" ;
				}
				($1!~/^#/){
					OFS="\t" ;
					$NF=$NF ;
#					density=$5*1000/($3-$2) ;
					density=($7+$8+$9)*1000/($3-$2) ;
					printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t%.0f\n",density ;
				}' \
			> "${sample}.${experiment}.v${atlas_ver}.clusters.density.bed"


			# make non-stranded superclusters from short-range clusters ($window distance)

			echo -e "#CHR\tSUPCLUS_START\tSUPCLUS_END\tSUPCLUS_ID\tCLUSTER_ID\tCLUSTER_NB_NRR\tCLUSTER_STRAND\tCLUSTER_DENSITY" \
			> "${sample}.${experiment}.superclusters.${window_name}kb.non-stranded.bed"

			sort -k1,1 -k2,2n "${sample}.${experiment}.v${atlas_ver}.clusters.density.bed" \
			| bedtools merge -d $window -c 4,5,6,7 -o collapse,collapse,collapse,collapse -i - \
			| sort -k1,1 -k2,2n \
			| awk -v sample=${sample} -v expe=${experiment} -v win=${window_name} '{
				OFS="\t" ;
				printf $1 "\t" $2 "\t" $3 "\t" sample "_" win "kb_" expe "_ID_%.6d\t",k ;
				printf $4 "\t" $5 "\t" $6 "\t" $7 "\n" ;
				k++ ;
			}' \
			>> "${sample}.${experiment}.superclusters.${window_name}kb.non-stranded.bed"

			# identify for each supercluster the best short cluster based on density of reads with an expected structure

			awk -v sample=${sample} -v expe=${experiment} -v win=${window_name} '($1!~/^#/){
				OFS="\t" ;
				split($5,name,",") ;
				split($8,density,",") ;
				max=0 ;
				best="" ;
				for (i in name) {
					if (density[i]>max) {
						max=density[i] ;
						best=name[i] ;
					} ;
				} ;
				if (best !="") {print best} ;
			}' "${sample}.${experiment}.superclusters.${window_name}kb.non-stranded.bed" \
			| sort -k1,1 \
			> "${sample}.${experiment}.v${atlas_ver}.best.clusters.${window_name}kb.list"

			# extract best insertions using cluster names
			awk '$1~/^#/' "../${sample}.${experiment}.v${atlas_ver}.insertions.full.bed" \
			> "${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.full.bed"

			# http://unix.stackexchange.com/questions/110645/select-lines-from-text-file-which-have-ids-listed-in-another-file
			grep -Fwf "${sample}.${experiment}.v${atlas_ver}.best.clusters.${window_name}kb.list" "../${sample}.${experiment}.v${atlas_ver}.insertions.full.bed" \
			>> "${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.full.bed"

			cut -f1-6 "${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.full.bed" \
			> "${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.true.bed"
	done
done

# generate pooled 5' and 3' atlas-seq datasets to annotate insertions (eg, when a 5' atlas cluster is only found in some samples but not all)
for experiment in ${exp_type[@]}
do
	echo "" > tmp.txt
	# pool all samples for a given experimental type (5' or 3' atlas-seq)
	for sample in ${sample_name[@]}
	do
		cat "${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.true.bed" \
		>> tmp.txt
	done

	sort -k1,1 -k2,2n tmp.txt \
	| uniq \
	| awk '($0!="") {OFS="\t"; print $1,$2,$3,$4,$5,$6}' \
	> "pooled.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.true.bed"

	cp "pooled.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.true.bed" ../
done

# cleanup
cp *.v${atlas_ver}.best.insertions.${window_name}kb.true.bed ../
cp *.v${atlas_ver}.best.insertions.${window_name}kb.full.bed ../
cp pooled*.bed ../
rm tmp*

########## Part 2: peak annotation

cd ${current_dir}
# for experiment in ${exp_type[@]}
for experiment in "3atlas"
do
	for sample in ${sample_name[@]}
	do
		echo -e "*** Annotating ${experiment} peaks for ${sample}"

		# annotate best insertion
		if [ ${experiment} == "5atlas" ];
			then
				atlas-annotate-multi.sh \
					-f "${current_dir}/pooled.3atlas.v${atlas_ver}.best.insertions.${window_name}kb.true.bed" \
					-a "${annotation_dir}/${annot5}" \
					"${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.full.bed"
			else
				atlas-annotate-multi.sh \
					-f "${current_dir}/pooled.5atlas.v${atlas_ver}.best.insertions.${window_name}kb.true.bed" \
					-a "${annotation_dir}/${annot3}" \
					"${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.full.bed"
				# identify potential stretches of A above putative L1 insertions

				echo -e "" > tmp.tab
				bedtools slop -s -l 25  -r 0 -i "${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.true.bed" -g "${ref_genome_dir}/hg19.genome" \
				| bedtools getfasta -s -name -tab -fi "${ref_genome_dir}/${ref_genome}" -bed - -fo tmp.tab

				awk '{OFS="\t"; print $1, gsub(/[Aa]{5}/,"",$2)}' tmp.tab \
				| sort -k1,1 \
				> tmp.Acount.tab

				sort -k4,4 "${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.full.annotated.bed" \
				> tmp.L1.bed

				echo -e $( awk '($1~/^#/){OFS="\t"; $NF=$NF; print $0,"A_stretch"}' tmp.L1.bed ) \
				| awk '{OFS="\t"; $NF=$NF; print $0}' \
				> "${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.full.annotated.bed"

				join -t $'\t' -1 4 -2 1 tmp.L1.bed tmp.Acount.tab \
				| awk '($1!~/^#/){OFS="\t"; printf $2 "\t" $3 "\t" $4 "\t" $1; for (i=5; i<=NF; i++) {printf "\t" $i}; printf "\n"}' \
				| sort -k1,1 -k2,2n \
				>> "${sample}.${experiment}.v${atlas_ver}.best.insertions.${window_name}kb.full.annotated.bed"
		fi
	done
done
rm tmp*
