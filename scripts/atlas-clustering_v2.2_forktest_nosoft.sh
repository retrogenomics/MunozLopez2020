#!/bin/bash

####### This is a stable version. Do not modify without changing the name

####### Dependencies:
# seqtk			1.0-r75-dirty	https://github.com/lh3/seqtk
# cutadapt		1.9.2.dev0		https://github.com/marcelm/cutadapt
# bwa			0.7.10-r789		https://github.com/lh3/bwa
# picard tools	1.136			http://broadinstitute.github.io/picard/
# bedtools		2.25.0			https://github.com/arq5x/bedtools2
# gnu grep/awk

####### External .atlas.conf file should be placed in home folder
# Example:

	#	# default folder and dependency files
	#	# this file should be located in the home directory folder of the user
	#
	#	# where picard tools executable files have been installed
	#	path_to_picard="/usr/local/picard"
	#
	# 	# where MIRA assembler binaries have been installed
	# 	path_to_mira="/Applications/Biotools/mira_4.0.2/bin"
	#
	#	# where useful scripts are saved
	#	path_to_useful_scripts="$HOME/Lab/bioinfo/pipelines/useful_scripts"
	#
	#	# where the reference genome files and indexes are stored
	#	ref_genome_dir="$HOME/Lab/bioinfo/references/human"
	#	ref_genome="hg19.fa"
	#
	#	# where the results are stored (a new folder is created for each analysis in the following folder)
	#	# !!!!! It should be created in the USER home folder (write rights)
	#	results_root="$HOME/Lab/bioinfo/results"

####### Release history
# 2.2 release notes (29/04/2016): stable version.
#	- add an additional step after read cluster generation to define L1 insertion point (L1 call) for 5' ATLAS (now both 5' and 3' ATLAS generate insertion point coordinates and strand)
# 	- multiple corrections for the new de novo assembly option to refine breakpoint:
#		- verify is cluster might corresponds to reference L1 (using -r argument file). If yes, do not perform de novo assembly
# 		- verify that contig could be made. If not, use cluster instead
# 		- verify that contig maps as cluster. If not, use cluster instead
# 		- verify that contig is softclipped. If not, use cluster instead
#	- output files are more homogenous with and without the -a option
#
# 2.1 release notes (12/04/2016)
#	- add an additional step after read cluster generation to define L1 insertion point (L1 call) for 3' ATLAS only
#	- for 3' atlas, two possibilities are offered (through the -a option):
#			* if -a is not used (default, quick): L1 insertion point is defined as the 3' end of cluster (insertion point is not defined at nucleotide resolution)
#			* if -a is used (much slower, reference L1 bed file required): L1 insertion point is refined for each cluster through:
#				- ATLAS-seq read de novo assembly;
#				- consensus sequence generation; and
#				- remapping with soft-clipping (insertion point is defined at nucleotide resolution).
#	- new argument (-r) introduced to specify a reference L1 .bed file (required if -a used, optional if -a not used, useless if -n)
#
# 2.0 release notes (17/03/2016)
#	- major reorganization of the script (each sample is processed entirely before passing to the next sample)
#	- major rewriting of 3' ATLAS-seq read structure analysis (read names are tagged by cutadapt upon removal or not of L1 or pA, and no longer mapped independently)
#	- demultiplexing of barcoded samples is now achieved only by cutadapt (and no longer by first fastx sorting, then cutadapt trimming)
#	- replace fastx by seqtk to reverse complement sequences
#	- add an option for subsampling (-s)
#	- cleanup the script comments and output
#	- output a script.log file with the script run to obtain the results
#	- add compatibility for bedtools v2.25 (the merge bedtools with -s option now reports strand as 4th column; and coverage tool has a modified behaviour)
#	- display pipeline running time as HH:MM:SS instead of seconds
#
# 1.8 release notes (25/05/2015)
#	- load default folders from parameter files atlas.param
#	- calculate read statistics directly by counting reads in output files instead of parsing cutadapt log
#	- added option to control multi-threading (mostly for bwa)
#
# 1.7 release notes (03/05/2015)
#	- adapt to run on ProfessorX (portability improvements)
#
# 1.6 release notes (13/03/2015)
#	- display spurious priming statistics (100 read structure)
#	- change read structure display (NB_NRR instead of %)
# 	- change cutadapt parameters for polyA trimming allowing 1 mismatch for trailing (T)8
#	- change cutadapt quality trimming q17 to q12 (to avoid loosing pA tail due to loss of quality in homopolymers)
#	- remove the decimals for TPM RPM display
#
# 1.5 release notes
# 	- added the -n for L1-neo experiments
#
# 1.4 release notes
#	- changed the way read structure is recorded and displayed in each cluster
# 	- calculate a normalized value of cluster size (RPM=reads assigned per millions of mapped reads (redundant), TPM=unique tags assigned per millions of mapped non-redundant)
#
# 1.3 release notes
#	- name of the script changed to atlas-clustering
#	- name of the full output bed file changed from ...full.bed to ...clusters.full.bed
#	- name of the standard output bed file changed from ...clusters.display.bed to ...clusters.true.bed
#	- added a header with track type in the outputed .bedgraph file
#	- corrected a bug leading to a missing column in the clusters.full.bed file
#
# 1.2_rc1 release notes
#	- input data file and barcode file as argument
#
# 1.1_rc4 release notes
# 	- clean up comments
#	- split the IRREGULAR reads into ARTEFACTUAL (100) and IRREGULAR (010,110,101)
#
# 1.1_rc3 release notes
#	- change the primer,L1,polyA counts to aberrant (100,010,110,101), neutral (000,001), high-confidence (011,111) reads in "${barcode_name[$i]}.12.aligned.noduplicate.triminfo.pooled.bed"
#	- update bedtools syntax to take into account changes from versions >2.2 and deprecated options.
#
# 1.1_rc2 release notes
#	- update to ensure compatibility with bedtools versions > 2.20 (bedtools merge options changed)
#	- added a step which excludes clusters containing reads with spurious priming (reads containing a trimmed primer but not L1 or pA sequence) and keep only the aligned reads from these clusters
#
# 1.0 release notes
# 	- barcode trimming with cutadapt has now no error allowed in the barcode since barcode splitting is done with 0 mismatch allowed
#	- primer, L1 and polyA trimming from 3' end now requires at least 8 overlapping nucleotides (-O 8 option in cutadapt)
#	- BWA verbose level changed from 3=message (default) to 1=error (-v 1 option)
#	- for all Picard tools lines, verbosity level option was set to VERBOSITY=ERROR
#	- intermediate files are now generated in a tmp folder, which is deleted at the end of the analysis to save space.
#	- corrected a bug that add XX:Z: flags (trimming information) to header of sam files
#	- save read and cluster statistics throughout the pipeline and generate a log file
#	- corrected a bug that duplicates the header of the final bed file for display
#	- compute bedgraph coverage files from non-redundant and mapped reads

####### known issues
#	- if a 5' ATLAS-seq fastq file is used but the -f option has been forgotten, the script is still running but errors appear at clustering due to empty files
#	- -a option not (yet) implemented for 5' atlas (with -f option).

####### load preferences and define global settings
script_name="atlas-clustering"
script_version='2.2forktest_nosoft'

start_time=`date +%s`
day=$(date +"[%d-%m-%Y] [%T]")

LC_NUMERIC_OLD=$LC_NUMERIC
export LC_NUMERIC="en_US.UTF-8"

starline=$( printf "*%.0s" {1..105} ) # separator for output
step=1	# store progress through the pipeline

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

###### script argument parsing and processing

# set defaults argument values
atlas5p=0
neo=0
distance=100 # window (bp) for merging reads into clusters
barcode_file=""
threads=4
sampling=1
assembly=0
reference_L1=""
use_L1_ref=0

# store usage explanations
USAGE="\
$script_name v$script_version:\tanalysis of ATLAS-seq runs from barcode demultiplexing to clustering. \n\n\
usage:\t$( basename $0 ) [options] -b barcode_file.txt input_file.fastq \n\
options:\n\
\t-h Print this help menu. \n\
\t-v What version of $script_name are you using? \n\
\t-f To indicate a 5' ATLAS-seq experiment  \n\
\t-n To indicate an L1-neo 3' ATLAS-seq experiment \n\
\t-d Maximum distance between reads to be merged into cluster [default=$distance] \n\
\t-s Subsampling of input fastq file (no=1; or indicate fraction of reads to consider, e.g. 0.01) [default=$sampling] \n\
\t-a Perform de novo local assembly of untrimmed reads to refine insertion point \n\
     (slow, requires -r argument, only used for 3' ATLAS-seq experiments) [default=not used] \n\
\t-r Reference L1 .bed file (required with -a option) \n\
\t-t Number of threads used for mapping [default=$threads]\n\
"

# parse script arguments
while getopts 'hvb:d:t:s:r:fnq' opt ; do
	case $opt in
		f) atlas5p=1 ;;
		n) neo=1 ;;
		a) assembly=1 ;;
		d) distance=$OPTARG ;;
		b) barcode_file=$OPTARG ;;
		t) threads=$OPTARG ;;
		s) sampling=$OPTARG ;;
		r) reference_L1=$OPTARG ;;
		h) echo -e "\n$USAGE"; exit 1 ;;
		v) echo -e "${script_name} v${script_version}" ; exit 1 ;;
		\?) echo -e "\nInvalid option: -$OPTARG\n" >&2; echo -e $USAGE; exit 1 ;;
	esac
done

# skip over the processed options
shift $((OPTIND-1))
input_file="$1"

# check for mandatory positional parameters

if [[ -z "${input_file}" || ! -f "${input_file}" ]];
then
	echo -e "\nInput file not specified or not existing.\n";
	echo -e $USAGE ; exit 1
fi

if [[ -z "${barcode_file}" || ! -f "${barcode_file}" ]];
then
	echo -e "\nBarcode file not specified or not existing.\n";
	echo -e $USAGE ; exit 1
fi

if [ $neo -eq 1 ];
then
	if [[ -n "${reference_L1}" ]];
	then
		if [[ -f "${reference_L1}" ]];
		then
			echo -e "\nYou specified a reference L1 file, but it won't be used with the -n option.\n";
		else
			echo -e "\nThe reference L1 file (${reference_L1}) could not be found, but it is useless anyway with the -n option.\nError ignored.\n";
		fi
	fi
else
	if [ $assembly -eq 0 ];
	then
		if [[ -n "${reference_L1}" ]];
		then
			if [[ -f "${reference_L1}" ]];
			then
				use_ref_L1=1;
			else
				echo -e "\nThe reference L1 file (${reference_L1}) could not be found.\nL1 calling will be performed without the -r option.\n";
			fi
		fi
	else
		if [[ -n "${reference_L1}" ]];
		then
			if [[ -f "${reference_L1}" ]];
			then
				use_ref_L1=1;
			else
				echo -e "\nThe reference L1 file (${reference_L1}) could not be found.\n";
				echo -e $USAGE ; exit 1
			fi
		else
			echo -e "\nReference L1 file (-r argument) not specified.\n";
			echo -e $USAGE ; exit 1
		fi
	fi
fi

# define ATLAS-seq primer and expected sequences, as well as minimal amplicon size
# ATLAS specific sequences
linker='GTGGCGGCCAGTATTCGTAGGAGGGCGCGTAGCATAGAACGT' # RBMSL2 +T from A-tailing
polyA='TTTTTTTTTTTT' # poly(A) sequence (reverse complement)
RB5PA2='TGGAAATGCAGAAATCACCG'
L1_5end_rc='TCTTCTGCGTCGCTCACGCTGGGAGCTGTAGACCGGAGCTGTTCCTATTCGGCCATCTTGGCTCCTCCCCC'
linker_rc='ACGTTCTATGCTACGCGCCCTCCTACGAATACTGGCCGCCAC'

if [ $atlas5p -eq 1 ];
then
	exp_name="5"
	min_amplicon=$RB5PA2$L1_5end_rc$linker_rc
	min_amplicon_size=${#min_amplicon}
else
	if [ $neo -eq 1 ];
	then
		exp_name="neo.3"
		RB3PA1='CGATACCGTAAGCCGAATTG'		# between SV40 promoter (Neo cassette) and end of L1 3'UTR (LOUXXX)
		RB3PA1_rc='CAATTCGGCTTACGGTATCG'
		L1_3end_rc='CGAACCCTGACGTCTTTATTATACTTTAAGTTTTAGGGTACATGTGCACATTGCC' # 3' end of L1 sequence from JM101/L1.3 Delta2 (reverse complement)
	else
		exp_name="3"
		RB3PA1='ATACCTAATGCTAGATGACACA'		# RB3PA1 primer from Badge et al. AJHG 2003, ends at the L1-Ta diagnotisc polymorphism (ACA)
		RB3PA1_rc='TGTGTCATCTAGCATTAGGTAT'
		L1_3end_rc='ATTATACTCTAAGTTTTAGGGTACATGTGCACATTGTGCAGGTTAGTTACATATGTATACATGTGCCATGCTGGTGCGCTGCACCCACTAA' # 3' end of L1 sequence without RB3PA1 (reverse complement)
	fi
	min_amplicon=$linker
	min_amplicon_size=$(( ${#min_amplicon} + 25 ))
fi

# find and print the directory name for data, reference and barcode files

data_dir=$( cd "$( dirname "${input_file}" )"; pwd )
data_name=$( basename "${input_file}" )

barcode_dir=$( cd "$( dirname "${barcode_file}" )"; pwd )
barcode=$( basename "${barcode_file}" )
nb_barcodes=$( grep -v -e '^#' "${barcode_dir}/$barcode" | wc -l )

refL1_dir=$( cd "$( dirname "${reference_L1}" )"; pwd )
refL1_name=$( basename "${reference_L1}" )

# obtain the name and sequence of barcodes and store them in bash tables

i=1
grep -v -e '^#' "${barcode_dir}/$barcode" \
> barcode_cleaned.txt

while read aLine ;
do
	barcode_name[$i]=$( echo $aLine | awk '{print $1}' )
	barcode_seq[$i]=$( echo $aLine | awk '{print $2}' )
	sample_name[$i]=$( echo $aLine | awk '{print $3}' )
	i=$(($i+1))
done < barcode_cleaned.txt
rm barcode_cleaned.txt

# create result directory with date and name of samples
results_dir="${results_root}/$( date +"%y%m%d_%H%M%S_" )${exp_name}atlas_${script_version}"
for i in `seq 1 $nb_barcodes`
do
	results_dir="${results_dir}_${sample_name[$i]}"
done

mkdir "${results_dir}"
mkdir "${results_dir}/tmp"

# print header
echo -ne "\n\
$starline\n\
$day ATLAS-seq ANALYSIS PIPELINE v$script_version \n\
$starline\n\
${exp_name}' ATLAS-seq experiment \n\
Sequencing data:\t${data_dir}/${data_name} \n\
Reference genome file:\t${ref_genome_dir}/${ref_genome} \n\
Barcodes:\t\t${barcode_dir}/${barcode} \n\
Samples: \n"

for i in $( seq 1 ${nb_barcodes} )
do
	echo -e "\t- ${barcode_name[$i]}: ${sample_name[$i]}"
done
echo -e "$starline"

# start generating statistic output file
output_stat="\
$starline\n\
$day ATLAS-seq ANALYSIS PIPELINE v${script_version} \n\
$starline\n\
${exp_name}' ATLAS-seq experiment \n\
Sequencing data:\t${data_dir}/${data_name} \n\
Reference genome file:\t${ref_genome_dir}/${ref_genome} \n\
Barcodes:\t\t${barcode_dir}/${barcode} \n\
"

####### test if subsampling required and if yes, generate subsampled input file
printf "[Step $step - Data loading]\n"
if [ $( echo "$sampling<1" | bc -l ) -eq 1 ];
then
	input="${data_dir}/${data_name}"
	data_dir="${results_dir}/tmp"
	data_name="sub.${data_name}"
	seqtk sample "$input" "$sampling" > "${data_dir}/${data_name}"
fi
(( step++ ))

####### barcode demultiplexing using cutadapt
cd "${results_dir}/tmp"
printf "[Step $step - Demultiplexing]\n"

# create cutadapt adapter line
adapt=""
for i in `seq 1 $nb_barcodes`
do
	adapt+=" -g ${barcode_name[$i]}=^${barcode_seq[$i]} "
done

# split fastq file according to barcode sequence (and remove barcode) with cutadapt
cutadapt -e 0.10 -q 10 $adapt \
	--untrimmed-output=discarded_missing_bc.fastq \
	-o {name}.01.trimmed_bc.fastq \
	"${data_dir}/${data_name}" \
&> /dev/null
touch discarded_missing_bc.fastq

# output statistics
total_reads=$(( $( wc -l "${data_dir}/${data_name}" | awk '{print $1}' ) / 4 ))
discarded_missing_bc=$(( $( wc -l discarded_missing_bc.fastq | awk '{print $1}' ) / 4 ))

output_stat+="\
$starline\n\
Total processed reads:\t${total_reads} \n\
  - no barcode found:\t${discarded_missing_bc}/${total_reads} \n\
"

for i in `seq 1 $nb_barcodes`
do
	barcode_out[$i]=$(( $( wc -l "${barcode_name[$i]}.01.trimmed_bc.fastq" | awk '{print $1}') / 4 ))
	output_stat+="  - ${sample_name[$i]} (${barcode_name[$i]}):\t${barcode_out[$i]}/${total_reads}\n"
done
output_stat+="$starline\n"

(( step++ ))

####### pipeline branching depending on type of ATLAS experiment (5' or 3')
if [ $atlas5p -eq 1 ];
then
	####### 5' ATLAS-seq
	####### process each sample one after the other
	for i in $( seq 1 $nb_barcodes );
	do
		cd "${results_dir}/tmp"
		printf "[Step $step - Processing sample ${sample_name[$i]}]\n"

		# remove L1-specific primer sequence from the reads with cutadapt and too short
		printf "  - Trim target-specific primer..."
		cutadapt -e 0.12 -q 10 -m ${min_amplicon_size} -g "^$RB5PA2" \
			--info-file="${barcode_name[$i]}.02b.primer.trimming.tab" \
			--too-short-output="${barcode_name[$i]}.02c.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.02d.discarded.missing.primer.fastq" \
			-o "${barcode_name[$i]}.02a.primer.trimming.fastq" \
			"${barcode_name[$i]}.01.trimmed_bc.fastq" \
		&> /dev/null

		# store statistics
		primer_in[$i]=${barcode_out[$i]}
		primer_out[$i]=$(( $( wc -l ${barcode_name[$i]}.02a.primer.trimming.fastq | awk '{print $1}') / 4 ))
		primer_too_short[$i]=$(( $( wc -l ${barcode_name[$i]}.02c.discarded_tooshort.fastq | awk '{print $1}') / 4 ))
		size_selected_read_count[$i]=$(( ${primer_in[$i]} - ${primer_too_short[$i]} ))

		printf "Done \n"

		# remove 5' L1 sequence from the reads with cutadapt (only for reads in which an L1-specific primer was found).
		printf "  - Trim L1 5' end..."
		cutadapt -e 0.12 -q 10 -m 25 -g $L1_5end_rc \
			--info-file="${barcode_name[$i]}.03b.5prime.trimming.tab" \
			--too-short-output="${barcode_name[$i]}.03c.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.03d.discarded.missing.L1.5end.fastq" \
			-o "${barcode_name[$i]}.03a.5prime.trimming.fastq" \
			"${barcode_name[$i]}.02a.primer.trimming.fastq"  \
		&> /dev/null

		# store statistics
		L1_end_in[$i]=${primer_out[$i]}
		L1_end_out[$i]=$(( $( wc -l ${barcode_name[$i]}.03a.5prime.trimming.fastq | awk '{print $1}') / 4 ))

		printf "Done \n"

		# trim ATLAS-linker at 3' end of reads with cutadapt
			# note that reads not reaching the linker sequence are discarded since they cannot be used to eliminate PCR duplicates
		printf "  - Trim ATLAS linker..."
		cutadapt -e 0.12 -q 10 -m 25 -a $linker_rc \
			--info-file="${barcode_name[$i]}.04b.3prime.trimming.tab" \
			--too-short-output="${barcode_name[$i]}.04c.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.04d.discarded.missing.adapter.fastq" \
			-o "${barcode_name[$i]}.04a.3prime.trimming.fastq" \
			"${barcode_name[$i]}.03a.5prime.trimming.fastq" \
		&> /dev/null

		# store statistics
		linker_in[$i]=${L1_end_out[$i]}
		linker_out[$i]=$(( $( wc -l ${barcode_name[$i]}.04a.3prime.trimming.fastq | awk '{print $1}') / 4 ))

		printf "Done \n"

		# Reverse complement using seqtk
			# This step is necessary to remove PCR duplicates using Picard.
			# 5' ATLAS reads reaching the linker are reversed complemented to start from the linker junction.
			# Picard will next eliminate those starting from the same position (after mapping), assuming they correspond to the same DNA fragmentation/linker ligation event.
			# Note: Reads not reaching the linker are discarded, thus lowering the coverage of the 5' junctions.
		printf "  - Reverse complement reads..."
		seqtk seq -r "${barcode_name[$i]}.04a.3prime.trimming.fastq" > "${barcode_name[$i]}.04f.3prime.trimming.rc.fastq"

		printf "Done \n"

		# map reads with bwa-mem (-t:nb of core used; -M:for compatibility with Picard MarkDuplicates; -v = verbosity mode)
		printf "  - Map reads on ${ref_genome}..."
		bwa mem -t $threads -M -v 1 "${ref_genome_dir}/${ref_genome}" "${barcode_name[$i]}.04f.3prime.trimming.rc.fastq" \
		2>/dev/null \
		| samtools view -q 20 -F 260 -bu - \
		| java -Xmx8g -jar "${path_to_picard}/picard.jar" SortSam \
			INPUT=/dev/stdin \
			OUTPUT="${barcode_name[$i]}.05.aligned.bam" \
			SORT_ORDER=coordinate \
			QUIET=true \
			VERBOSITY=ERROR \
			VALIDATION_STRINGENCY=LENIENT \
			CREATE_INDEX=true \

		# test if the alignment .sam file is empty (only header, no reads) for samtools compatibility
		if [[ $( samtools view ${barcode_name[$i]}.05.aligned.bam | tail -1 ) =~ ^@.* ]];
			then
				mapped_count[$i]=0
			else
				# calculate and store the count of uniquely mapped reads after trimming of 5' linker and 3' primer-L1-pA in an array
				mapped_count[$i]=$( samtools view -c "${barcode_name[$i]}.05.aligned.bam" )
		fi
		echo

		printf "  - Map reads on ${ref_genome}...Done \n"

		# mark and remove duplicate reads (identical start position)
		printf "  - Remove PCR duplicates..."
		java -Xmx8g -jar "${path_to_picard}/picard.jar" MarkDuplicates \
			INPUT="${barcode_name[$i]}.05.aligned.bam"\
			OUTPUT="${barcode_name[$i]}.06.aligned.noduplicate.bam" \
			METRICS_FILE="/dev/null/" \
			REMOVE_DUPLICATES=true \
			DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES\
			ASSUME_SORTED=true \
			VALIDATION_STRINGENCY=LENIENT \
			QUIET=true \
			VERBOSITY=ERROR \
			CREATE_INDEX=true

		# test if the alignment .sam file after duplicate removal is not empty (only header, no reads) for samtools compatibility. The "!" after the "if" is to inverse the condition
		if ! [[ $( samtools view "${barcode_name[$i]}.06.aligned.noduplicate.bam" | tail -1 ) =~ ^@.* ]];
		then
			unique_read_count[$i]=$( samtools view -c "${barcode_name[$i]}.06.aligned.noduplicate.bam" )
		fi

		printf "Done \n"

		# generate clusters (using bamtobed then merge, with number and name of non-redundant reads)
		printf "  - Generate clusters (.bed)..."

		# test if the .bam file for a given barcode exists and create bed files
		if [[ -f "${barcode_name[$i]}.06.aligned.noduplicate.bam" ]];
			then
				# create a bed file with an entry per read and sort it
				# options:
				#	-k1,1 sort alphabetically on the first field [using -k1 only would use from field 1 to end of line] and resolve ties with second field (-k2,2n) sorted numerically.
				# This .bed file contains only uniquely mapped, non-redundant reads, for which L1-specific primer, L1 5' end and linker were trimmed, and longer than 25nt.
				bedtools bamtobed -i "${barcode_name[$i]}.06.aligned.noduplicate.bam" \
				| sort -k1,1 -k2,2n \
				> "${barcode_name[$i]}.07.aligned.noduplicate.sorted.bed"

				# create a bed file with an entry per read and sort it (including redundant reads for rpm calculation)
				samtools view -q 20 -bu -F 260 "${barcode_name[$i]}.05.aligned.bam" \
				| bedtools bamtobed -i - \
				| sort -k1,1 -k2,2n \
				> "${barcode_name[$i]}.07.aligned.sorted.bed"
			else
				# create an empty .bed file if no .bam is found
				touch "${barcode_name[$i]}.07.aligned.noduplicate.sorted.bed"
				touch "${barcode_name[$i]}.07.aligned.sorted.bed"
		fi

		# merge overlapping reads or reads distant from less than $distance into clusters and sort them
			# warning: bedtools merge -s output has changed in bedtools v2.25
			# options used:
			# 	-s = force strandness (only merge reads in the same orientation)
			# 	-c = columns to operate
			# 	-o = operations to process on columns
			#	-d = max distance between 2 entries to allow merging
		bedtools merge -i "${barcode_name[$i]}.07.aligned.noduplicate.sorted.bed" -s -d $distance -c 4,4,6 -o distinct,count,distinct \
		| awk '$1!~/^#/ {OFS="\t"; print $1,$2,$3,$5,$6,$7}' \
		| sort -k1,1 -k2,2n \
		> "${barcode_name[$i]}.08a.clusters.sorted.bed"

		# calculate nb of reads (redundant) per cluster
			# warning: bedtools coverage input file convention has changed (-a and -b are like other tools since bedtools v2.24 )
		bedtools coverage -sorted -s -counts -a "${barcode_name[$i]}.08a.clusters.sorted.bed" -b "${barcode_name[$i]}.07.aligned.sorted.bed"\
		> "${barcode_name[$i]}.08b.clusters.sorted.with_total_reads.bed"

		# calculate and store the count of clusters
		cluster_count[$i]=$( grep -v -e "^#" ${barcode_name[$i]}.08a.clusters.sorted.bed | wc -l )

		# create a header for the final bed file
			# 	NB_READS = number of redundant reads
			#	RPM = reads assigned per millions of mapped reads (redundant)
			# 	TPM = unique tags assigned per millions of mapped non-redundant reads
		echo -e "#CHR\tSTART\tEND\tCLUSTER_ID\tNB_NRR\tCLUSTER_STRAND\tNB_READS\tRPM\tTPM" \
		> "${barcode_name[$i]}.09.numbered_clusters.bed"

		# replace the read names by a unique cluster ID in the name field ($4) and calculate RPM and TPM
		awk -v k=1 -v n=${mapped_count[$i]} -v u=${unique_read_count[$i]} '{OFS="\t"; print $1,$2,$3,k,$5,$6,$7,$7*1000000/n,$5*1000000/u; k++}' "${barcode_name[$i]}.08b.clusters.sorted.with_total_reads.bed"\
		| awk -v sample=${sample_name[$i]} '{OFS="\t"; printf $1 "\t" $2 "\t" $3 "\t"; printf sample"_5ATLAS_"; printf "%.4d\t",$4; printf $5 "\t" $6 "\t" $7 "\t"; printf "%.0f\t",$8; printf "%.0f\n",$9}'\
		>> "${barcode_name[$i]}.09.numbered_clusters.bed"

		# call insertion point (for 5' ATLAS-seq, clusters are in the same orientation as L1)
		awk 'BEGIN {printf "#CHR\tINS_START\tINS_END\tCLUSTER_ID\tNB_NRR\tINS_STRAND\tNB_READS\tRPM\tTPM\n"}
			($1!~/^#/){
				OFS="\t";
				if ($6=="+") {
					printf $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6;
					for (i=7; i<=NF; i++) {printf "\t"$i}; printf "\n"}
				else {
					printf $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6;
					for (i=7; i<=NF; i++) {printf "\t"$i}; printf "\n"}
				}' "${barcode_name[$i]}.09.numbered_clusters.bed" \
		> "${barcode_name[$i]}.11.insertions.bed"

		# generate minimal bed files for visualization purpose
		cut -f1-6 "${barcode_name[$i]}.09.numbered_clusters.bed" \
		> "${barcode_name[$i]}.10.cluster.display.bed"

		cut -f1-6 "${barcode_name[$i]}.11.insertions.bed" \
		> "${barcode_name[$i]}.12.insertions.true.bed"

		printf "Done \n"

		# generate a bedgraph file for visualisation purpose
		printf "  - Compute coverage (.bedgraph)..."

		# plus strand bedgraph
		bedgraph_head_p="track type=bedGraph name=${sample_name[$i]}.${exp_name}atlas(+) description="" visibility=full color=150,0,0 priority=20 autoScale=off alwaysZero=on maxHeightPixels=32 graphType=bar viewLimits=0:300 yLineMark=0 yLineOnOff=on windowingFunction=mean smoothingWindow=3"
		echo -e "$bedgraph_head_p" \
		> "${barcode_name[$i]}.06.aligned.noduplicate.plus.bedgraph"
		bedtools genomecov -ibam "${barcode_name[$i]}.06.aligned.noduplicate.bam" -bg -strand + \
		>> "${barcode_name[$i]}.06.aligned.noduplicate.plus.bedgraph"

		# minus strand bedgraph
		bedgraph_head_m="track type=bedGraph name=${sample_name[$i]}.${exp_name}atlas(-) description="" visibility=full color=150,0,0 priority=20 autoScale=off alwaysZero=on maxHeightPixels=32 graphType=bar viewLimits=0:300 yLineMark=0 yLineOnOff=on windowingFunction=mean smoothingWindow=3"
		echo -e "$bedgraph_head_m" \
		> "${barcode_name[$i]}.06.aligned.noduplicate.minus.bedgraph"
		bedtools genomecov -ibam "${barcode_name[$i]}.06.aligned.noduplicate.bam" -bg -strand - \
		>> "${barcode_name[$i]}.06.aligned.noduplicate.minus.bedgraph"

		printf "Done \n"

		# generates the final result and log files and clean up the temporary files
		printf "  - Calculate statistics (.log) and cleanup..."
		output_stat[$i]="\
Sample:\t\t\t${sample_name[$i]} \n\
Barcode:\t\t\t${barcode_name[$i]} \n\
Reads in the run:\t\t\t${total_reads} \n\
Barcode sorting:\t${total_reads}\t->\t${barcode_out[$i]} \n\
Size selection:\t${barcode_out[$i]}\t->\t${size_selected_read_count[$i]} \n\
L1-specific primer trimming:\t${size_selected_read_count[$i]}\t->\t${primer_out[$i]} \n\
L1 5' end trimming:\t${L1_end_in[$i]}\t->\t${L1_end_out[$i]} \n\
Linker trimming:\t${linker_in[$i]}\t->\t${linker_out[$i]} \n\
Unambigously mapped:\t${linker_out[$i]}\t->\t${mapped_count[$i]} \n\
Non-redundant reads:\t${mapped_count[$i]}\t->\t${unique_read_count[$i]} \n\
Clusters:\t${unique_read_count[$i]}\t->\t${cluster_count[$i]} \n\
"
		echo -e "${output_stat[$i]}" > "${barcode_name[$i]}.log"

		# move and rename files to keep
		cd "${results_dir}"
		echo -e "$output_stat" > "global.${exp_name}atlas.v${script_version}.log"
		mv "tmp/${barcode_name[$i]}.log" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.log"
		mv "tmp/${barcode_name[$i]}.06.aligned.noduplicate.bam" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.bam"
		mv "tmp/${barcode_name[$i]}.06.aligned.noduplicate.bai" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.bai"
		mv "tmp/${barcode_name[$i]}.09.numbered_clusters.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.clusters.full.bed"
		mv "tmp/${barcode_name[$i]}.10.cluster.display.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.clusters.true.bed"
		mv "tmp/${barcode_name[$i]}.11.insertions.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.insertions.full.bed"
		mv "tmp/${barcode_name[$i]}.12.insertions.true.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.insertions.true.bed"
		mv "tmp/${barcode_name[$i]}.06.aligned.noduplicate.plus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.plus.bedgraph"
		mv "tmp/${barcode_name[$i]}.06.aligned.noduplicate.minus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.minus.bedgraph"

		printf "Done \n"
		(( step++ ))
	done
else
	####### 3' ATLAS-seq
	####### process each sample one after the other
	for i in $( seq 1 $nb_barcodes );
	do
		cd "${results_dir}/tmp"
		printf "[Step $step - Processing sample ${sample_name[$i]}]\n"

		# trim ATLAS-linker at 5' end of reads with cutadapt
		printf "  - Trim ATLAS linker..."
		cutadapt -e 0.12 -q 10 -m ${min_amplicon_size} -g "^$linker" \
			--info-file="${barcode_name[$i]}.02b.trimmed_linker.tab" \
			--too-short-output="${barcode_name[$i]}.02c.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.02d.discarded_missing_adapter.fastq" \
			-o "${barcode_name[$i]}.02a.trimmed_linker.fastq" \
			"${barcode_name[$i]}.01.trimmed_bc.fastq" \
		&> /dev/null

		# store statistics
		linker_in[$i]=${barcode_out[$i]}
		linker_out[$i]=$(( $( wc -l "${barcode_name[$i]}.02a.trimmed_linker.fastq" | awk '{print $1}') / 4 ))
		linker_too_short[$i]=$(( $( wc -l "${barcode_name[$i]}.02c.discarded_tooshort.fastq" | awk '{print $1}' ) / 4 ))
		size_selected_read_count[$i]=$(( ${linker_in[$i]} - ${linker_too_short[$i]} ))

		printf "Done \n"

		# trim L1-specific primer, L1 sequence and polyA at 3' end of reads using cutadapt
		printf "  - Trim target-specific primer, L1 3' end, and poly(A)..."

		# look for and remove the L1-specific primer sequence, then the L1 sequence, and finally the polyA at the 3' end of the reads with cutadapt
			# Each possibility is encoded in the file name by a 3-character code:
			# 	- first character corresponds to L1-specific primer sequence
			# 	- second character corresponds to L1 sequence
			#	- third character corresponds to polyA
			# For each of the 3 positions: x=unknown, X=interrogated, 0=absent, 1=present.
			# Ex: 	in the following command, ...Xxx... = interrogation of presence/absence of L1-specific primer
			# 		1xx = reads with L1-specific primer found and trimmed
			#		0xx = reads with no L1-specific primer found
			#		etc...

		cutadapt -e 0.12 -m 25 -q 10 -O 6 -a "${RB3PA1_rc}" \
			--info-file="${barcode_name[$i]}.04.Xxx.3prime_trimming.tab" \
			--too-short-output="${barcode_name[$i]}.04.Xxx.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.04.0xx.fastq" \
			-o "${barcode_name[$i]}.04.1xx.fastq" \
			"${barcode_name[$i]}.02a.trimmed_linker.fastq" \
		&> /dev/null

		cutadapt -e 0.12 -m 25 -q 10 -O ${#L1_3end_rc} -a "${L1_3end_rc}" \
			--info-file="${barcode_name[$i]}.04.1Xx.3prime_trimming.tab" \
			--too-short-output="${barcode_name[$i]}.04.1Xx.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.04.10x.fastq" \
			-o "${barcode_name[$i]}.04.11x.fastq" \
			"${barcode_name[$i]}.04.1xx.fastq" \
		&> /dev/null

		cutadapt -e 0.12 -m 25 -q 10 -O 8 -a "T{100}" \
			--info-file="${barcode_name[$i]}.04.11X.3prime_trimming.tab" \
			--too-short-output="${barcode_name[$i]}.04.11X.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.04.110.fastq" \
			-o "${barcode_name[$i]}.04.111.fastq" \
			"${barcode_name[$i]}.04.11x.fastq" \
		&> /dev/null

		cutadapt -e 0.12 -m 25 -q 10 -O 8 -a "T{100}" \
			--info-file="${barcode_name[$i]}.04.10X.3prime_trimming.tab" \
			--too-short-output="${barcode_name[$i]}.04.10X.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.04.100.fastq" \
			-o "${barcode_name[$i]}.04.101.fastq" \
			"${barcode_name[$i]}.04.10x.fastq" \
		&> /dev/null

		cutadapt -e 0.12 -m 25 -q 10 -O 6 -a "${L1_3end_rc}" \
			--info-file="${barcode_name[$i]}.04.0Xx.3prime_trimming.tab" \
			--too-short-output="${barcode_name[$i]}.04.0Xx.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.04.00x.fastq" \
			-o "${barcode_name[$i]}.04.01x.fastq" \
			"${barcode_name[$i]}.04.0xx.fastq" \
		&> /dev/null

		cutadapt -e 0.12 -m 25 -q 10 -O 8 -a "T{100}" \
			--info-file="${barcode_name[$i]}.04.01X.3prime_trimming.tab" \
			--too-short-output="${barcode_name[$i]}.04.01X.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.04.010.fastq" \
			-o "${barcode_name[$i]}.04.011.fastq" \
			"${barcode_name[$i]}.04.01x.fastq" \
		&> /dev/null

		cutadapt -e 0.12 -m 25 -q 10 -O 6 -a "T{100}" \
			--info-file="${barcode_name[$i]}.04.00X.3prime_trimming.tab" \
			--too-short-output="${barcode_name[$i]}.04.00X.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.04.000.fastq" \
			-o "${barcode_name[$i]}.04.001.fastq" \
		"${barcode_name[$i]}.04.00x.fastq" \
		&> /dev/null

		# store statistics
			# store number of reads too short after trimming (below the 25nt threshold, excluded from the output file due to the -m 25 option)
			# loop to extract and store the number of reads too short after trimming for each of the possible read structures
			# note: not all these statistics are actually displayed or written in the .log files (kept for historical reasons and in case needed in the future)
		L1pA_too_short[$i]=0
		for j in "Xxx" "1Xx" "11X" "10X" "0Xx" "01X" "00X" ;
		do
			L1pA_too_short[$i]=$(( ${L1pA_too_short[$i]} + $( wc -l "${barcode_name[$i]}.04.$j.discarded_tooshort.fastq" | awk '{print $1}') / 4 ))
		done

		# calculate and store the number of reads spanning the primer-L1-pA junction with a remaining flanking sequence ≥ 25 bp
		full_amplicon_count[$i]=$(( $( wc -l "${barcode_name[$i]}.04.111.fastq" | awk '{print $1}' ) / 4 ))

		# calculate and store the number of reads spanning the L1-pA junction with a remaining flanking sequence ≥ 25 bp, but not reaching the primer
		L1pA_junction_count[$i]=$(( $( wc -l "${barcode_name[$i]}.04.011.fastq" | awk '{print $1}' ) / 4 ))

		# calculate and store the number of reads reaching the pA junction with a remaining flanking sequence ≥ 25 bp, but not reaching the primer and the L1
		pA_junction_count[$i]=$(( $( wc -l "${barcode_name[$i]}.04.001.fastq" | awk '{print $1}' ) / 4 ))

		# calculate and store the number of reads reaching a canonical L1-pA or pA junction with a remaining flanking sequence ≥ 25 bp
		canonical_junction_count[$i]=$(( ${full_amplicon_count[$i]} + ${L1pA_junction_count[$i]} + ${pA_junction_count[$i]} ))

		# calculate and store the number of reads not reaching the pA junction with a flanking sequence ≥ 25 bp
		unreached_junction_count[$i]=$(( $( wc -l "${barcode_name[$i]}.04.000.fastq" | awk '{print $1}' ) / 4 ))

		# calculate and store the number of reads with weird structure (with primer but no L1 or pA, with L1 but without polyA, etc...) with a remaining flanking sequence ≥ 25 bp
		unexpected_trimming_count[$i]=$(( ( $( wc -l "${barcode_name[$i]}.04.100.fastq" | awk '{print $1}' ) + $( wc -l "${barcode_name[$i]}.04.010.fastq" | awk '{print $1}' ) + $( wc -l "${barcode_name[$i]}.04.110.fastq" | awk '{print $1}' ) + $( wc -l "${barcode_name[$i]}.04.101.fastq" | awk '{print $1}' ) ) / 4 ))

		# reads remaining after junction trimming
		junction_out[$i]=$(( ${canonical_junction_count[$i]} + ${unreached_junction_count[$i]} + ${unexpected_trimming_count[$i]} ))

		# output basic statistics in the shell
		# echo -e "Reads remaining after linker trimming:\t${linker_out[$i]}/$total_reads"
		# echo -e "Canonical L1-pA or pA junction were found and trimmed in:\t${canonical_junction_count[$i]}/${linker_out[$i]} reads."
		# echo -e "Junction was not reached in:\t${unreached_junction_count[$i]}/${linker_out[$i]} reads."
		# echo -e "Unexpected read structure was found in:\t${unexpected_trimming_count[$i]}/${linker_out[$i]} reads."
		# echo -e "Reads too short after junction trimming and discarded:\t${L1pA_too_short[$i]}/${linker_out[$i]}."
		# echo -e "Reads remaining after junction trimming:\t${junction_out[$i]}/${linker_out[$i]}."
		# echo

		printf "Done \n"

		# append tag (read structure) to read names and concatenate fastq files
		printf "  - Analyze read structure..."
		echo -ne "" > ${barcode_name[$i]}.05.triminfo.fastq

		for j in "000" "001" "011" "111" "010" "110" "101" "100";
		do
			# add an optional tag "XX" reflecting the elements (RB3PA1, L1end and PolyA) that have been trimmed in 3', encoded in a string of 3 booleans.
			# Ex: "001" means: RB3PA1 not trimmed, L1end not trimmed, polyA trimmed.
			awk -v j=$j '{if (NR%4==1) {print $1" XX:Z:"j} else {print $0}}' ${barcode_name[$i]}.04.$j.fastq \
			>> ${barcode_name[$i]}.05.triminfo.fastq
		done

		printf "Done \n"

		# map reads with bwa-mem
			# options: 	-t = nb of core used
			#			-M = for compatibility with Picard MarkDuplicates
			#			-v = verbosity mode
			#			-C = append FASTA/FASTQ comment to output
		printf "  - Map reads on ${ref_genome}..."
		bwa mem -t $threads -C -M -v 1 "${ref_genome_dir}/${ref_genome}" "${barcode_name[$i]}.05.triminfo.fastq" \
		2>/dev/null \
		| awk '($1~/^@/) || ($6!~/[0-9]+S/)' \
		| samtools view -q 20 -F 260 -bu - \
		| java -Xmx8g -jar "${path_to_picard}/picard.jar" SortSam \
			INPUT=/dev/stdin \
			OUTPUT="${barcode_name[$i]}.06.triminfo.aligned.bam" \
			SORT_ORDER=coordinate \
			QUIET=true \
			VERBOSITY=ERROR \
			VALIDATION_STRINGENCY=LENIENT \
			CREATE_INDEX=true \

		# test if the alignment .sam file is empty (only header, no reads) for samtools compatibility
		if [[ $( samtools view "${barcode_name[$i]}.06.triminfo.aligned.bam" | tail -1 ) =~ ^@.* ]];
			then
				mapped_count[$i]=0
			else
				# calculate and store the count of uniquely mapped reads after trimming of 5' linker and 3' primer-L1-pA in an array
				mapped_count[$i]=$( samtools view -c "${barcode_name[$i]}.06.triminfo.aligned.bam" )
		fi

		printf "Done \n"

		# remove PCR duplicates
		printf "  - Remove PCR duplicates..."

		# mark and remove duplicate reads (identical start position), option java -Xmx8g related to memory usage (4g= 4Gb)
		java -Xmx8g -jar "${path_to_picard}/picard.jar" MarkDuplicates \
			INPUT="${barcode_name[$i]}.06.triminfo.aligned.bam"\
			OUTPUT="${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" \
			METRICS_FILE="/dev/null/" \
			REMOVE_DUPLICATES=true \
			DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
			ASSUME_SORTED=true \
			VALIDATION_STRINGENCY=LENIENT \
			QUIET=true \
			VERBOSITY=ERROR \
			CREATE_INDEX=true

		# test if the alignment .sam file after duplicate removal is not empty (only header, no reads) for samtools compatibility. The "!" after the "if" is to inverse the condition
		if ! [[ $( samtools view "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" | tail -1 ) =~ ^@.* ]];
		then
			unique_read_count[$i]=$( samtools view -c "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" )
		fi

		printf "Done \n"

		# generate clusters (using bamtobed then merge, with number and name of non-redundant reads)
		printf "  - Generate clusters (.bed)..."

		# test if the .bam file for a given barcode exists and create bed files
		if [[ -f "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" ]];
		then
			# create a bed file with an entry per non-redundant read and sort it
			echo -ne "" > "${barcode_name[$i]}.08.triminfo.aligned.noduplicate.bed"
			for j in "000" "001" "011" "111" "010" "110" "101" "100";
			do
				samtools view -h "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" \
				| grep -E "(^@|XX\:Z\:$j)" \
				| samtools view -h -b - \
				| bedtools bamtobed -i - \
				| awk -v j=$j '{OFS="\t"; print $0, j}' \
				>> "${barcode_name[$i]}.08.triminfo.aligned.noduplicate.bed"
			done
			sort -k1,1 -k2,2n "${barcode_name[$i]}.08.triminfo.aligned.noduplicate.bed" \
			> "${barcode_name[$i]}.09a.triminfo.aligned.noduplicate.sorted.bed"

			# create a bed file with an entry per read and sort it (including redundant reads for rpm calculation)
			samtools view -q 20 -bu -F 260 "${barcode_name[$i]}.06.triminfo.aligned.bam" \
			| bedtools bamtobed -i - \
			| sort -k1,1 -k2,2n \
			> "${barcode_name[$i]}.09b.triminfo.aligned.sorted.bed"
		else
			# create an empty .bed file if no .bam is found
			touch "${barcode_name[$i]}.09a.triminfo.aligned.noduplicate.sorted.bed"
			touch "${barcode_name[$i]}.09b.triminfo.aligned.sorted.bed"
		fi

		# merge overlapping reads or reads distant from less than $distance into clusters and sort them
			# warning: bedtools merge -s output has changed in bedtools v2.25
			# options used:
			# 	-s = force strandness (only merge reads in the same orientation)
			# 	-c = columns to operate
			# 	-o = operations to process on columns
			#	-d = max distance between 2 entries to allow merging
		bedtools merge -i "${barcode_name[$i]}.09a.triminfo.aligned.noduplicate.sorted.bed" -s -d $distance -c 4,4,6 -o distinct,count,distinct \
		| awk '$1!~/^#/ {OFS="\t"; print $1,$2,$3,$5,$6,$7}' \
		| sort -k1,1 -k2,2n \
		> "${barcode_name[$i]}.10.clusters.sorted.bed"

		# calculate nb of reads (redundant) per cluster
			# warning: bedtools coverage input file convention has changed (-a and -b are like other tools since bedtools v2.24 )
		bedtools coverage -sorted -s -counts -a "${barcode_name[$i]}.10.clusters.sorted.bed" -b "${barcode_name[$i]}.09b.triminfo.aligned.sorted.bed"\
		> "${barcode_name[$i]}.10a.clusters.sorted.bed"

		# sequentially add to each cluster the count of non-redundant reads with a given read structure (010, 110, etc...)
		k="a"
		for j in "010" "110" "101" "100" "000" "001" "011" "111";
		do
			# Take the pseudo bed file with all clusters and annotate sequentially with the count of non-redundant reads with each type of read structure
				# intersect reports overlaps between two feature files
				# options used:
				# -wa	Write the original entry in A for each overlap
				# -c	For each entry in A, report the number of overlaps with B.
				# -s	Force strandness
			samtools view -h "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" \
			| grep -E "(^@|XX\:Z\:$j)" \
			| samtools view -h -b - \
			| bedtools intersect -wa -c -s -a "${barcode_name[$i]}.10$k.clusters.sorted.bed" -b - \
			> "${barcode_name[$i]}.10`echo $k | tr 'a-y' 'b-z'`.clusters.sorted.bed"

			# increment the letter for file name
			k=$( echo $k | tr 'a-y' 'b-z' )
		done

		# create a reduced cluster file with the number of reads falling in the different trimming categories (grouped by type)
			# [L1-pA-flank] = 111 ($15), 011 ($14)
			# [pA-flank] = 001 ($13)
			# [flank] = 000 ($12)
			# [irregular] = 010 ($8), 110 ($9), 101 ($10)
			# [spurious] = 100 ($11)

		awk 'BEGIN {OFS="\t"; print "#CHR","START","END","READS_ID","NB_NRR","CLUSTER_STRAND","L1-pA-flank","pA-flank","flank","irregular","spurious","NB_READS"} \
			$0!~/^#/ {print $1,$2,$3,$4,$5,$6,$14+$15,$13,$12,$8+$9+$10,$11,$7}' "${barcode_name[$i]}.10i.clusters.sorted.bed" \
			| sort -k1,1 -k2,2n \
			> "${barcode_name[$i]}.11.clusters.bed"

		# create a header for the final bed file
			# NB_READS = number of redundant reads
			# RPM = reads assigned per millions of mapped reads (redundant)
			# TPM=unique tags assigned per millions of mapped non-redundant reads
		echo -e "#CHR\tSTART\tEND\tCLUSTER_ID\tNB_NRR\tCLUSTER_STRAND\tL1-pA-flank\tpA-flank\tflank\tirregular\tspurious\tNB_READS\tRPM\tTPM" > "${barcode_name[$i]}.12.numbered_clusters.bed"

		# replace the read names by a unique cluster ID in the name field ($4)
		awk -v k=1 -v n=${mapped_count[$i]} -v u=${unique_read_count[$i]} '$0!~/^#/ {OFS="\t"; print $1,$2,$3,k,$5,$6,$7,$8,$9,$10,$11,$12,$12*1000000/n,$5*1000000/u; k++}' "${barcode_name[$i]}.11.clusters.bed" \
		| awk -v sample=${sample_name[$i]} '{OFS="\t"; printf $1 "\t" $2 "\t" $3 "\t"; printf sample"_3ATLAS_";printf "%.4d\t",$4; printf $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t"; printf "%.0f\t",$13; printf "%.0f\n",$14}' \
		>> "${barcode_name[$i]}.12.numbered_clusters.bed"

		# calculate and store the count of clusters
		cluster_count[$i]=$( grep -c -v -e "^#" ${barcode_name[$i]}.12.numbered_clusters.bed )

		# generate a minimal bed file for visualization purpose
		echo -e "#CHR\tSTART\tEND\tCLUSTER_ID\tNB_NRR\tCLUSTER_STRAND" \
		> "${barcode_name[$i]}.13.cluster.display.bed"

		awk '$0!~/#/ {OFS="\t"; print $1,$2,$3,$4,$5,$6}' "${barcode_name[$i]}.12.numbered_clusters.bed" \
		>> "${barcode_name[$i]}.13.cluster.display.bed"

		# generate a bam file and index for reads in clusters
		samtools view -bu -h -L "${barcode_name[$i]}.13.cluster.display.bed" "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" \
		| java -Xmx8g -jar "${path_to_picard}/picard.jar" SortSam \
			INPUT=/dev/stdin \
			OUTPUT="${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.bam" \
			SORT_ORDER=coordinate \
			QUIET=true \
			VERBOSITY=ERROR \
			VALIDATION_STRINGENCY=LENIENT \
			CREATE_INDEX=true

		printf "Done \n"

		# generate a bedgraph file for visualisation purpose
		printf "  - Compute coverage (.bedgraph)..."

		# plus strand bedgraph
		bedgraph_head_p="track type=bedGraph name=${sample_name[$i]}.${exp_name}atlas(+) description="" visibility=full color=0,150,0 priority=20 autoScale=off alwaysZero=on maxHeightPixels=32 graphType=bar viewLimits=0:300 yLineMark=0 yLineOnOff=on windowingFunction=mean smoothingWindow=3"
		echo -e "$bedgraph_head_p" \
		> "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.plus.bedgraph"

		bedtools genomecov -ibam "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.bam" -bg -strand + \
		>> "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.plus.bedgraph"

		# minus strand bedgraph
		bedgraph_head_m="track type=bedGraph name=${sample_name[$i]}.${exp_name}atlas(-) description="" visibility=full color=0,150,0 priority=20 autoScale=off alwaysZero=on maxHeightPixels=32 graphType=bar viewLimits=0:300 yLineMark=0 yLineOnOff=on windowingFunction=mean smoothingWindow=3"
		echo -e "$bedgraph_head_m" \
		> "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.minus.bedgraph"

		bedtools genomecov -ibam "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.bam" -bg -strand - \
		>> "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.minus.bedgraph"

		printf "Done \n"

		# L1 calling
		if [ $assembly -eq 0 ];
		then
			####### Default quick option selected (extremity of cluster)
			# add insertion point coordinates instead of cluster coordinates and sorting
			# (for 3' ATLAS-seq, clusters are in the opposite orientation as compared to L1)
			printf "  - Call L1 insertions from ATLAS-seq clusters..."

			awk '{
				OFS="\t";
				if (NF < 6) {$NF=$NF; print $0}
				else if (NF == 6) {
					if ($1 ~ /^#/) {printf $1"\tINS_START\tINS_END\t"$4"\t"$5"\tINS_STRAND\t" ; for (j=7;j<=NF-1;++j) printf $j"\t"; printf $NF"\n"}
					else {
						if ($6=="+") {print $1,$3,$3,$4,$5,"-"}
						else {print $1,$2,$2,$4,$5,"+"}
						}
				}
				else {
					if ($1 ~ /^#/) {printf $1"\tINS_START\tINS_END\t"$4"\t"$5"\tINS_STRAND\t" ; for (j=7;j<=NF-1;++j) printf $j"\t"; printf $NF"\n"}
					else {
						if ($6=="+") {printf $1"\t"$3"\t"$3"\t"$4"\t"$5"\t-\t"; for (j=7;j<=NF-1;++j) printf $j"\t"; printf $NF"\n"}
						else {printf $1"\t"$2"\t"$2"\t"$4"\t"$5"\t+\t"; for (j=7;j<=NF-1;++j) printf $j"\t"; printf $NF"\n"}
						}
				}
			}' "${barcode_name[$i]}.12.numbered_clusters.bed" \
			| sort -k1,1 -k2,2n \
			> "${barcode_name[$i]}.16a.L1_insertions.breakpoint.bed"

			cut -f1-6 "${barcode_name[$i]}.16a.L1_insertions.breakpoint.bed" \
			| sort -k1,1 -k2,2n \
			> "${barcode_name[$i]}.16a.L1_insertions.breakpoint.true.bed"

			printf "Done \n"

		else
			####### Slow option (-p option, insertion point refined precisely through de novo assembly of junctions)
			# generate the list of clusters that will pass through reassembly

			bedtools closest -D a -iu -S -a "${barcode_name[$i]}.13.cluster.display.bed" -b "${refL1_dir}/${refL1_name}" \
			| awk '$7!="." && $13<=200 {print $4}' \
			>"${barcode_name[$i]}.putative.ref.clusters.list"

			bedtools closest -D a -iu -S -a "${barcode_name[$i]}.13.cluster.display.bed" -b "${refL1_dir}/${refL1_name}" \
			| awk '$7=="." || $13>200 {print $4}' \
			>"${barcode_name[$i]}.putative.nonref.clusters.list"

			j=1

			while read aLine ;
			do
				cluster_name[$j]=$aLine
				j=$(( $j + 1 ))
			done < "${barcode_name[$i]}.putative.nonref.clusters.list"

			nb_clusters=$( wc -l "${barcode_name[$i]}.putative.nonref.clusters.list" | awk '{print $1}' )
			cluster_coord=""

			# initialize progress bar
			bar=$( printf '=%.0s' {1..25} )
			total_bar_length=${#bar}
			printf "  - De novo assembly of junctions...\n"

			# initialize contig fasta file
			echo "" > "${barcode_name[$i]}.contig.fasta"

			# loop for each cluster
			for j in $( seq 1 ${nb_clusters} ); do

				# extract cluster coordinates
				cluster_coord=$( grep -e ${cluster_name[$j]} "${barcode_name[$i]}.13.cluster.display.bed" | awk '{print $1":"$2"-"$3}' )

				# extract read names
				samtools view -S "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.bam" "${cluster_coord}" \
				| awk '{print $1}' \
				> reads_list.txt

				# extract sequences from linker_trimmed fastq files based on read names
				seqtk subseq "${barcode_name[$i]}.02a.trimmed_linker.fastq" reads_list.txt | seqtk sample - 100 > reads_cluster.fastq

				# generate contigs & rename contigs to include cluster name (single read clusters are considered as contigs)
				reads_nb=$(( $( wc -l reads_cluster.fastq | awk '{print $1}' ) / 4 ))

				# increment progress bar
				current_bar_length=$(( $j * ${total_bar_length} / ${nb_clusters} ))
				printf "\r    [%${#nb_clusters}u/${nb_clusters}] [%-${total_bar_length}s] [${cluster_name[$j]}: %3u read(s) considered for assembly]" "$j" "${bar:0:${current_bar_length}}" "${reads_nb}"

				# perform de novo assembly of untrimmed reads in the same cluster with MIRA (unless a single read is present in the cluster)
				echo -ne "" > "${barcode_name[$i]}.no.contig.assembled.list"	# store list of clusters for which de novo assembly did not produced any contig
				echo -ne "" > "${barcode_name[$i]}.contig.assembled.list"		# store list of clusters for which de novo assembly produced one or more contigs (includes single-read contigs)
				if [ $reads_nb -gt 1 ] ;
					then
						# create manifest file for MIRA
						manifest_content="\
						project = atlas \n\
						job = denovo,genome,accurate \n\
						parameters = -NAG_AND_WARN:check_maxreadnamelength=no -NAG_AND_WARN:check_average_coverage=warn -NW:cnfs=warn IONTOR_SETTINGS -AS:mrpc=2 -AS:mrl=25 -CL:lccf=no -CL:lccb=no -CL:cpat=no -CL:qc=yes -CL:bsqc=yes -CL:bsqcmq=15 -CL:bsqcwl=10 -CL:qcwl=10 -CL:qcmq=15 -CL:pec=no \n\
						\n\
						readgroup = ${cluster_name[$j]} \n\
						data = reads_cluster.fastq \n\
						technology = iontor"
						echo -e "${manifest_content}" > manifest.txt

						# run mira
						"${path_to_mira}"/mira manifest.txt 2>error.log 1>output.log
						path_to_mira_output="atlas_assembly/atlas_d_results"

						# test assembly error (no contig produced)
						if [ ! -s "${path_to_mira_output}/atlas_out.maf" ];
							then
								echo -e " warning: no MIRA output for ${cluster_name[$j]}"
								echo -e "${cluster_name[$j]}" >> "${barcode_name[$i]}.no.contig.assembled.list"
								touch "${path_to_mira_output}/atlas_out.maf"
						fi

						# extract contigs and keep the one made with the highest number of reads
						awk -v n=${cluster_name[$j]} 'BEGIN{OFS="\t"}
							($0~/^CO/) { printf(n"|%s\t", $2) }
							($0~/^(NR|CS)/) { printf("%s\t", $2) }
							($0~/^CQ/) { printf("%s", $2) }
							($0~/^\/\//) { printf("\n") }' "${path_to_mira_output}/atlas_out.maf" \
						| sort -k2,2rn \
						| head -1 \
						| awk '{OFS="\n"; split($1,a,"|"); gsub(/\*/,"",$3); print ">"a[1]"|"$2,$3}'\
						>> "${barcode_name[$i]}.contig.fasta"
						echo -e "${cluster_name[$j]}" >>  "${barcode_name[$i]}.contig.assembled.list"
					else
						# if cluster contains a single read, consider it as a contig and add it to contig multifasta file
						awk -v n=${cluster_name[$j]} 'BEGIN{OFS="\t"; printf(n"|atlas_c1\t1\t") }
							NR % 4 == 2 {printf("%s\t", $1)}
							NR % 4 == 0 {printf("%s\n", $1)}' reads_cluster.fastq \
						| awk '{OFS="\n"; split($1,a,"|"); gsub(/\*/,"",$3); print ">"a[1]"|1",$3}'\
						>> "${barcode_name[$i]}.contig.fasta"
						echo -e "${cluster_name[$j]}" >>  "${barcode_name[$i]}.contig.assembled.list"
				fi
			done

			printf "\n    ...Done\n"

			# mapping contig on hg19 with BWA
			printf "  - Identify junction breakpoints..."

			bwa mem -t $threads -M -v 1 "${ref_genome_dir}/${ref_genome}" "${barcode_name[$i]}.contig.fasta" \
			2>/dev/null \
			| samtools view -bu - \
			| java -Xmx8g -jar "${path_to_picard}/picard.jar" SortSam \
			INPUT=/dev/stdin \
			OUTPUT="${barcode_name[$i]}.14.contig.aligned.bam" \
			SORT_ORDER=coordinate \
			QUIET=true \
			VERBOSITY=ERROR \
			VALIDATION_STRINGENCY=LENIENT \
			CREATE_INDEX=true

			# extract list of clusters for which mapping of contig and cluster differs
			grep -f "${barcode_name[$i]}.contig.assembled.list" "${barcode_name[$i]}.13.cluster.display.bed" \
			| bedtools intersect -s -wa -loj -a - -b "${barcode_name[$i]}.14.contig.aligned.bam" \
			| awk '$7=="." {print $4}' \
			>  "${barcode_name[$i]}.contig.mismapped.list"

			# extract list of clusters for which contig was not soft-clipped (reference L1 or junction not reached)
			samtools view -F260 -q20 "${barcode_name[$i]}.14.contig.aligned.bam" \
			| grep -v -f "${barcode_name[$i]}.contig.mismapped.list" \
			| awk '{OFS="\t"; if (($2==16 && $6!~/^[0-9]+S/) || ($2==0 && $6!~/[0-9]+S$/)) {split($1,contig_name,"|"); print contig_name[1]}}' \
			>  "${barcode_name[$i]}.contig.not_softclipped.list"

			# extract list of clusters for which contig was soft-clipped (non-reference L1)
			samtools view -F260 -q20 "${barcode_name[$i]}.14.contig.aligned.bam" \
			| grep -v -f "${barcode_name[$i]}.contig.mismapped.list" \
			| awk '{OFS="\t"; if (($2==16 && $6~/^[0-9]+S/) || ($2==0 && $6~/[0-9]+S$/)) {split($1,contig_name,"|"); print contig_name[1]}}' \
			>  "${barcode_name[$i]}.contig.softclipped.list"

			# generate bed files for the different categories of clusters/contigs

			# -if assembly was not performed (cluster <= 200 bp distant from reference L1HS)
			grep -f "${barcode_name[$i]}.putative.ref.clusters.list" "${barcode_name[$i]}.13.cluster.display.bed" \
			| awk '{OFS="\t"; if ($6=="+") {print $1,$3,$3,$4,$5,"-","assembly_not_done|putative_ref_L1HS"} else {print $1,$2,$2,$4,$5,"+","assembly_not_done|putative_ref_L1HS"}}' \
			> "${barcode_name[$i]}.15a.noassembly.bed"

			# - if no contig made, use cluster coordinates to define insertion point
			grep -f "${barcode_name[$i]}.no.contig.assembled.list" "${barcode_name[$i]}.13.cluster.display.bed" \
			| awk '{OFS="\t"; if ($6=="+") {print $1,$3,$3,$4,$5,"-","no_contig_assembled"} else {print $1,$2,$2,$4,$5,"+","no_contig_assembled"}}' \
			> "${barcode_name[$i]}.15b.nocontig.bed"

			# - if contig made but not properly mapped, use cluster coordinates to define insertion point
			grep -f "${barcode_name[$i]}.contig.mismapped.list" "${barcode_name[$i]}.13.cluster.display.bed" \
			| awk '{OFS="\t"; if ($6=="+") {print $1,$3,$3,$4,$5,"-","contig_assembled|mismapped"} else {print $1,$2,$2,$4,$5,"+","contig_assembled|mismapped"}}' \
			> "${barcode_name[$i]}.15c.mismapped_contig.bed"

			# - if contig made but mapped contig was not softclipped at the putative junction (i.e. reference L1 or junction not reached), use cluster coordinates to define insertion point
			grep -f "${barcode_name[$i]}.contig.not_softclipped.list" "${barcode_name[$i]}.13.cluster.display.bed" \
			| awk '{OFS="\t"; if ($6=="+") {print $1,$3,$3,$4,$5,"-","contig_assembled|not_soft-clipped"} else {print $1,$2,$2,$4,$5,"+","contig_assembled|not_soft-clipped"}}' \
			> "${barcode_name[$i]}.15d.not_softclipped_contig.bed"

			# - if contig made and softclipped at the putative junction (i.e. putative non-reference L1), use contig breakpoint to define insertion point
			# generate a bed file with insertion points (and L1 orientation instead of cluster orientation)
			samtools view -b -F260 -q20 "${barcode_name[$i]}.14.contig.aligned.bam" \
			| bedtools bamtobed \
			| grep -f "${barcode_name[$i]}.contig.softclipped.list" \
			| awk '{OFS="\t"; split($4,a,"|"); if ($6=="+") {print $1,$3,$3,a[1],a[2],"-","contig_assembled|softclipped"} else {print $1,$2,$2,a[1],a[2],"+","contig_assembled|softclipped"}}' \
			| sort -k4,4 \
			> "${barcode_name[$i]}.15e.softclipped_contig.bed"

			# for soft-clipped contig, correct score in bed file (nb of reads used for de novo assembly replaced by NB_NRR)
			join -t$'\t' -o 1.1 1.2 1.3 1.4 2.5 1.6 1.7 -1 4 -2 4 "${barcode_name[$i]}.15e.softclipped_contig.bed" "${barcode_name[$i]}.13.cluster.display.bed" \
			> "${barcode_name[$i]}.15e.softclipped_contig.corrected_score.bed"

			# concatenate all L1 insertion coordinates
			( echo -e "#CHR\tINS_START\tINS_STOP\tCLUSTER_ID\tNB_NRR\tINS_STRAND\tBREAKPOINT" ;\
			cat "${barcode_name[$i]}.15a.noassembly.bed" \
				"${barcode_name[$i]}.15b.nocontig.bed" \
				"${barcode_name[$i]}.15c.mismapped_contig.bed" \
				"${barcode_name[$i]}.15d.not_softclipped_contig.bed" \
				"${barcode_name[$i]}.15e.softclipped_contig.corrected_score.bed" )\
			| sort -k4,4 \
			| join -t$'\t' -o 1.1 1.2 1.3 1.4 1.5 1.6 2.7 2.8 2.9 2.10 2.11 2.12 2.13 2.14 1.7 -1 4 -2 4 - "${barcode_name[$i]}.12.numbered_clusters.bed" \
			| sort -k1,1 -k2,2n \
			> "${barcode_name[$i]}.16a.L1_insertions.breakpoint.bed"

			cut -f1-6 "${barcode_name[$i]}.16a.L1_insertions.breakpoint.bed" \
			| sort -k1,1 -k2,2n \
			> "${barcode_name[$i]}.16a.L1_insertions.breakpoint.true.bed"

			# generate a file with additional information (CIGAR, sequence, etc) associated to insertion ID
			# sense L1 (antisense cluster, sam flag 16): contig starts with softclipped sequence, output contig sequence and softclipped sequence
			# antisense L1 (sense cluster, sam flag 0): contig ends with softclipped sequence, output reverse complement of contig sequence and rev. comp. of softclipped seq.
			samtools view -F260 -q20 "${barcode_name[$i]}.14.contig.aligned.bam" \
			| awk 'BEGIN{
				tr["A"]="T";
				tr["T"]="A";
				tr["G"]="C";
				tr["C"]="G";
				tr["N"]="N";
				tr["Y"]="R";
				tr["R"]="Y";}
				{
				OFS="\t";
				if ($2==16) {
					if ($6~/^[0-9]+S/) {
						a=match($6,"S");
						l=substr($6,1,a-1);
						softclipped_seq=substr($10,1,l);
					}
					else {
						softclipped_seq=".";
					}
					split($1,contig_name,"|");
					print contig_name[1],$6,$10,softclipped_seq;
				}
				else {
					if ($2==0) {
						rev_comp_contig="";
						length_contig=split(toupper($10),seq,"");
						for (i=length_contig; i>=1; i--) {
							if (tr[seq[i]]!="") {
								rev_comp_contig=rev_comp_contig tr[seq[i]]
							}
							else {
								rev_comp_contig=rev_comp_contig "?"
							}
						}
						if ($6~/[0-9]+S$/) {
							a=match($6,/[0-9]+S$/);
							b=match($6,/S$/);
							length_softclipped=substr($6,a,b-a);
							rev_comp_softclipped=substr(rev_comp_contig, 1, length_softclipped);
						}
						else {
							rev_comp_softclipped=".";
						}
						split($1,contig_name,"|");
						print contig_name[1],$6,rev_comp_contig,rev_comp_softclipped;
					}
				}
			}'\
			| sort -k1,1 \
			> "${barcode_name[$i]}.16b.L1_insertions.junction.sequences.tab"

			# join L1 insertion points with sequence information
			join -t$'\t' -o 1.1 1.2 1.3 1.4 1.5 1.6 2.2 2.3 -1 4 -2 1 "${barcode_name[$i]}.16a.L1_insertions.breakpoint.true.bed" "${barcode_name[$i]}.16b.L1_insertions.junction.sequences.tab" \
			> "${barcode_name[$i]}.16c.L1_insertions.junction.sequences.forward.bed"

			printf "Done \n"

			# obtain sequence of +/-25bp around insertion site
# 			bedtools slop -b 25 -i refined_insertion_point.bed -g $ref_genome_dir/hg19.genome \
# 			| bedtools getfasta -s -name -fi $ref_genome_dir/$ref_genome -bed - -fo refined_insertion_point.fasta
		fi

		# generates the final result files
		printf "  - Calculate statistics (.log) and cleanup..."

		output_stat[$i]="\
Sample:\t\t\t${sample_name[$i]} \n\
Barcode:\t\t\t${barcode_name[$i]} \n\
Reads in the run:\t\t\t${total_reads} \n\
Barcode sorting:\t${total_reads}\t->\t${barcode_out[$i]} \n\
Size selection:\t${barcode_out[$i]}\t->\t${size_selected_read_count[$i]} \n\
Linker trimming:\t${size_selected_read_count[$i]}\t->\t${linker_out[$i]} \n\
L1 and pA trimming:\t${linker_out[$i]}\t->\t${junction_out[$i]}\t(including ${canonical_junction_count[$i]} possible junctions)\n\
Unambigously mapped:\t${junction_out[$i]}\t->\t${mapped_count[$i]} \n\
Non-redundant reads:\t${mapped_count[$i]}\t->\t${unique_read_count[$i]} \n\
Clusters:\t${unique_read_count[$i]}\t->\t${cluster_count[$i]} \n\
"
		echo -e "${output_stat[$i]}" > "${barcode_name[$i]}.log"

		# move and rename files to keep
		cd "${results_dir}"
		mv "tmp/${barcode_name[$i]}.log" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.log"
		mv "tmp/${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.bam" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.bam"
		mv "tmp/${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.bai" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.bai"
		mv "tmp/${barcode_name[$i]}.12.numbered_clusters.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.clusters.full.bed"
		mv "tmp/${barcode_name[$i]}.13.cluster.display.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.clusters.true.bed"
		mv "tmp/${barcode_name[$i]}.16a.L1_insertions.breakpoint.true.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.insertions.true.bed"
		mv "tmp/${barcode_name[$i]}.16a.L1_insertions.breakpoint.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.insertions.full.bed"
		mv "tmp/${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.plus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.plus.bedgraph"
		mv "tmp/${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.minus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.minus.bedgraph"

		printf "Done \n"
		(( step++ ))
	done
fi

####### restore system variables
export LC_NUMERIC=$LC_NUMERIC_OLD

####### display log files
cd "${results_dir}"
for i in $( seq 1 $nb_barcodes );
do
	output_stat+="\
$( cat "${sample_name[$i]}.${exp_name}atlas.v${script_version}.log" ) \n\
$starline \n\
"
done

# calculate runtime for the whole pipeline
end_time=`date +%s`
runtime=$( date -u -d @$(( end_time - start_time )) +"%T" )
day=$(date +"[%d-%m-%Y] [%T]")

output_stat+="\
$day \tRunning time: $runtime (hh:mm:ss) \n\
$starline \n\
"

echo -e "$output_stat" | tee "global.${exp_name}atlas.v${script_version}.log"
cat "$0" > "$( basename "$0" ).script.log"

####### delete intermediate files
rm -r tmp
