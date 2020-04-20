# MunozLopez2020
Scripts to call L1 insertions from single-cell ATLAS-seq experiments.

## Before starting

### Dependencies
- [cutadapt](https://github.com/marcelm/cutadapt) (tested version: 1.14)
- [bwa](https://github.com/lh3/bwa) (tested version: 0.7.16a)
- [Picard tools](http://broadinstitute.github.io/picard/) (tested version: 1.136)
- [bedtools](https://github.com/arq5x/bedtools2) (tested version: 2.25.0)
- [seqtk](https://github.com/lh3/seqtk) (tested version: 1.0; note that seqtk is only required if subsampling of sequencing data is used - to reduce time of analysis in tests)
- [GNU parallel](https://www.gnu.org/software/parallel/) (tested version: 20150522)
- GNU grep/awk

### Prepare the general configuration file and organize your project folders
Place and edit the `.atlas.conf` file in your home folder according to your configuration.
Prepare directories as follow:
```
project/
|-- annotations
|-- data
|   `-- atlas-seq
|-- results
`-- scripts
```
The `data/atlas-seq` folder should contain the .fastq files.

### Other requirements
- A human reference genome sequence (ex:`hg19.fa`) and its bwa index (ex: `hg19.fa.amb, .ann, .bwt, .fai, .pac, .sa`). Their location is indicated in the `.atlas.conf` file.
- adjust the paths of annotations files in the `annotations/annotations_sc_5atlas.txt` and `annotations/annotations_sc_5atlas.txt` files.

## Procedure

### Step 1: For each sequencing run, run the atlas-clustering script:
Starting from a sequencing run, this script:
- demultiplex the sequencing reads based on samples
- trim the reads and map them on the reference genome provided
- cluster sequencing reads and identify potential break points

1. example for a 3'-atlas-seq run
The barcode file `.bc` is a tabular text file with 3 columns (index name, index sequence, sample name). An example is provided in the annotation folder.
```
cd project/results
../scripts/atlas-clustering_v2.2_forktest_nosoft.sh \
	-d 0 \
	-b ../annotations/3pp.bc \
	../data/atlas-seq/single_cell_embryos/R05_INS-203.ATLAS-seq.E6T_E6C1_3prime.fastq
```

2. example for a 5'-atlas-seq run
The barcode file `5pp.bc` has 3 columns: index name, index sequence, sample name.
```
cd project/results
../scripts/atlas-clustering_v2.2_forktest_nosoft.sh \
	-d 0 \
	-f \
	-b ../annotations/5pp.bc \
	../data/atlas-seq/single_cell_embryos/R07_INS-208.ATLAS-seq_E6C1_5prime_E6C3_5prime.fastq
```

### Step 2: move all result files into same folder (../150319_pooled_encode_cell_lines_1.6)
```
mkdir -p project/results/pooled_single_cells
cp project/results/*_*atlas*/* project/results/pooled_single_cells/
```

### Step 3: run the peak calling and annotation script from the pooled folder
This script calls L1 peaks within large clusters of reads due to local amplification during the whole genome amplification process, and annotate them.
```
cd results/pooled_single_cells
project/scripts/atlas-seq_single_cells_Acount.sh
```

### Step 4: additional filters
This filtering step was added to remove artefactual amplifications in the WGA procedure.

```
cd results/pooled_single_cells
for file in *.3atlas.v2.2forktest_nosoft.best.insertions.10kb.full.annotated.bed ;
do
	name=$( echo -e $file | awk -F "." '{printf $1}' ) ;

	awk '($1~/^#/) || (($5>=3 && $8+$7>=1) && (($19~/L1HS\|Ta/) || (($19=="." && $21!~/L1PA/) && ($17=="." || $11==0 ) && ($NF<=1))))' ${name}.3atlas.v2.2forktest_nosoft.best.insertions.10kb.full.annotated.bed \
	> ${name}.3atlas.v2.2forktest_nosoft.best.insertions.10kb.full.annotated.filtered.bed;

	cut -f1-6 ${name}.3atlas.v2.2forktest_nosoft.best.insertions.10kb.full.annotated.filtered.bed \
	> ${name}.3atlas.v2.2forktest_nosoft.best.insertions.10kb.filtered.true.bed;
done
```
