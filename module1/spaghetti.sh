#############
# ARGUMENTS #
#############

# 1: input directory (required)

########################
# Read pre-processment #
########################

echo "Starting analysis"
echo "Pre-processing the reads"
echo "------------------------"
echo ""

###################
# Adapter removal #
###################

echo ""
echo "1. Adapter removal"
echo "------------------"

time (for i in $1*fastq
do
	porechop -t 12 -i $i -o ${i::-6}-porechop.fastq
done)

#################################
# Length filtering 1200-1800 pb #
#################################

echo ""
echo "2. Length filtering"
echo "-------------------"

time (for i in $1*porechop*
do
	cat $i | NanoFilt -l 1200 --maxlength 1800 > ${i::-6}-nanofilt.fastq
done)

######
# QC #
######

echo ""
echo "3. Quality check"
echo "----------------"

time(for i in $1*porechop-nanofilt*
do 
	NanoStat -t 12 --fastq $i > ${i::-6}-NanoStat.txt
done)

###################
# Chimera removal #
###################

echo ""
echo "4. Chimera removal"
echo "------------------"

time(for i in $1*nanofilt.fastq
do 
	minimap2 -x ava-ont -g 500 -t 12 $i $i > ${i%%.*}.paf
	yacrd -i ${i%%.*}.paf -o ${i%%.*}.yacrd -c 4 -n 0.4 scrubb -i $i -o ${i%%.*}.scrubb.fastq

	rm ${i%%.*}.paf
done)

###########################
# Read mapping (minimap2) #
###########################

# The script is based on minimap2's page (https://github.com/lh3/minimap2)
# Minimap2 has been used by other authors for 16S rRNA sequencing with ONT
# See: https://f1000research.com/articles/7-1755
# Also: https://elifesciences.org/articles/61504
# Also: https://www.sciencedirect.com/science/article/pii/S2001037019303745

# Silva 138 has been already indexed
# Silva 138 was download in Qiime2 format directly from: https://docs.qiime2.org/2020.8/data-resources/#marker-gene-reference-databases

#################
# Map sequences #
#################

echo ""
echo "5. Mapping"
echo "----------"

# Silva 138
time (for i in $1*scrubb.fastq
do echo "Mapping" $i; minimap2 -x map-ont -t 12 --secondary=no -K 10M ~/Descargas/Databases/Silva_138_qiime2/dna-sequences.mmi $i > $i.paf
done)

# For understanding PAF format: https://github.com/lh3/miniasm/blob/master/PAF.md

#############
# Filtering #
#############

echo ""
echo "6. Filtering the PAF files"
echo "--------------------------"

# Although secondary alignments are turned off, some query reads have been mapped to multiple database seqs.
# Let's filter the output file to just keep one alignment per read.
# We will keep the largest alignment.
# We will remove alignments <500 pb

time (for i in $1*.paf
do filterPAF.py -i $i > ${i::-4}-f.paf
done)

rm -r $1filteredPAFs
mkdir $1filteredPAFs
mv $1*-f.paf $1filteredPAFs

########################################
# Summarize & Merge filtered PAF files #
########################################

echo ""
echo "7. Summarizing the PAF files"
echo "----------------------------"

# It creates a table with the number of sequences assigned to each Database's ID for each sample
# It's something like a QIIME (v1) OTU table.

merfePAF.py -i $1filteredPAFs/ > $1otu_table.csv

####################################
# Create a phyloseq taxonomy table #
####################################

echo ""
echo "8. Creating the taxonomy table"
echo "------------------------------"

# This script creates the taxonomy table needed for loading the data into the phyloseq package.

taxonomyTable.py -i $1otu_table.csv -t ~/Descargas/Databases/Silva_138_qiime2/taxonomy.tsv > $1phyloseq_taxonomy.csv

