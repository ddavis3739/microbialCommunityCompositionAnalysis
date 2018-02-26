#!/bin/bash

## STEPS, following (https://link.springer.com/protocol/10.1007/8623_2016_228) for quality control
# 1) unzip
# 2) trim reads with sickle
# 3) error correct with SPAdes
# 4) join reads with PandaSeq
# 5) split based on primer sequence (should be done first because of trimming approaches)
# 6) create OTU tables

#create associated directories
mkdir unzipped
mkdir qualTrimmed
mkdir qualTrimmed/errorCorrected
mkdir allErrorCorrected
mkdir pairedSamples

cd pairedSamples
mkdir dereplicated
mkdir sortedSequences
mkdir centroids
mkdir chimeras
mkdir nonchimeras
mkdir OTU
mkdir rdp_assigned_taxonomy
mkdir log
cd ..

# 1) Copy and unzip files, must be in .fastq.gz format
cp */*/*.fastq.gz ./unzipped
gunzip unzipped/*.gz

# 2) Quality trim using both read directions (paired)
shopt -s globstar

for file in unzipped/*_R1.fastq;
do
	basefile=$(echo ${file#unzipped/}) 
	sickle pe -f $file \
	-r ${file%_R1.fastq}_R2.fastq \
	-t sanger -o qualTrimmed/${basefile%.fastq}_trimmed.fastq \
	-p qualTrimmed/${basefile%_R1.fastq}_R2_trimmed.fastq \
	-s qualTrimmed/${basefile%_R1.fastq}_singles.fastq
done

cd qualTrimmed

# 3) perform error correction with Spades BayesHammer 
# (outputs a folder (sample_n_errorCorrected) for each sample
for file in ./*_R1_trimmed.fastq;
do
	basefile=$(basename $file)
  	spades.py -o ./errorCorrected/${basefile%_R1_trimmed.fastq}_errorCorrected \
		--only-error-correction -1 $file -2 ${file%_R1_trimmed.fastq}_R2_trimmed.fastq \
		-t 30 -m 64 --disable-gzip-output;
done

# move files to more assesible location

cp */*/corrected/*R[1–2]*.fastq ../allErrorCorrected

# edit sequence headers for assembly into single reads

cd ../allErrorCorrected

for f in *R1*.fastq;
do 
	sed -i "/$(grep -m 1 "^@" $f | egrep -o "^[^:]+")/s/$/ 1:N:0:/" $f;
done

for f in *R2*.fastq;
do 
	sed -i "/$(grep -m 1 "^@" $f | egrep -o "^[^:]+")/s/$/ 2:N:0:/" $f;
done

# 4) paired end joining using pandaseq

for f in *R1*.fastq;
do 
	basefile=$(basename $file)
	pandaseq -f $f -r ${f%R1_trimmed.00.0_0.cor.fastq}R2_trimmed.00.0_0.cor.fastq \
		-A pear -B -p CCTACGGGNGGCWGCAG -q GACTACHVGGGTATCTAATCC \
		-w ../pairedSamples/${f%_R1_trimmed.00.0_0.cor.fastq}_paired.fasta \
		-g ../pairedSamples/${f%_R1_trimmed.00.0_0.cor.fastq}_log.txt;
done

# 5) Edit fasta header to include sample

cd ../pairedSamples

for f in *.fasta;
do 
	sample=">$(echo $f | awk -F "." '{print $1}')_"
	sed -i "s/>/$sample/" $f
done

# combine into one

mkdir combined

cat *.fasta > combined/seqs.fna

# trim out blank seqs
cd combined

awk '/^$/{for(x=NR-1;x<=NR;x++)d[x];}   
    {a[NR]=$0} 
    END{for(i=1;i<=NR;i++)if(!(i in d))print a[i]}'\
    seqs.fna > test.fna 

mv test.fna seqs.fna

cd ..

# 6) OTU de-novo on paired data using vsearch
for file in ./combined/*.fna
do
	basefile=$(basename $file)

	# type all commands in on one line
	# first dereplicate your QC'ed sequences. The minuniquesize parameter indicates
	# the smallest OTUs to accept, so a value of 2 would discard singleton sequences
	
	printf "\nDereplicating\n"
	vsearch --derep_fulllength $file --output dereplicated/$basefile --minuniquesize 2 --relabel OTU --fasta_width 0 --sizeout

	# the command also sorts sequences by their abundance, from most to least abundant

	printf "\nSorting\n"
	vsearch --sortbysize ./dereplicated/$basefile --output ./sortedSequences/$basefile

	# this command picks your OTU centroid sequences, use the sequences that are
	# output to assign taxonomy

	printf "\nClustering\n"
	vsearch --cluster_size ./sortedSequences/$basefile --id 0.97 --centroids ./centroids/$basefile

	# de novo chimera check

	printf "\nChimera filtering\n"
    	vsearch --uchime_denovo ./centroids/$basefile --chimeras ./chimeras/$basefile --nonchimeras ./nonchimeras/$basefile

	# this command takes your original QC'ed sequences and assigns them to the OTU
	# centroids that were created in the previous step.

	printf "\nOTU Table creation\n"
	vsearch --usearch_global $file --db ./nonchimeras/$basefile --id 0.97 --log ./log/${basefile%.fna}.log -biomout ./OTU/${basefile%.fna}.biom
	
	printf "\nAssign taxonomy\n" 
	java -Xmx1g -jar /usr/local/rdp_classifier/dist/classifier.jar classify \
		-o ./rdp_assigned_taxonomy/${basefile%.fna}_tax.txt -f fixrank \
		./nonchimeras/$basefile

	# Convert to .txt

	printf "\nConvert .biom to .txt\n"
	biom convert -i ./OTU/${basefile%.fna}.biom -o ./OTU/${basefile%.fna}_OTUtable.txt --to-tsv #--header-key taxonomy 
	
	# read into R and combine OTU table and taxonomy
	
	otuTable=$(echo OTU/${basefile%.fna}_OTUtable.txt)
	taxonomy=$(echo rdp_assigned_taxonomy/${basefile%.fna}_tax.txt)

	Rscript ../otu_table_merge.R $otuTable $taxonomy ${basefile%.fna}

done

# 7) deblur workflow, picks at higher threshold than OTU methods
# requires that deblur is loaded into cluster

cd pairedSamples