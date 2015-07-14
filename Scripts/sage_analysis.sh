#!/usr/bin/env bash

#Compute Transcripts with extended 3' regions if possible up to EXTENDED_LENGTH
SUPERCONTIGS_FASTA=../phytophthora_data/original/phytophthora_infestans_t30-4_1_supercontigs.fasta
TRANSCRIPTS_FASTA=../phytophthora_data/original/phytophthora_infestans_t30-4_1_transcripts.fasta
ANNOTATION_FILE=../phytophthora_data/original/phytophthora_infestans_t30-4_1_transcripts.gtf
EXTENDED_LENGTH=1000

perl Extract_targetRegions.pl -g $SUPERCONTIGS_FASTA -t $TRANSCRIPTS_FASTA -f $ANNOTATION_FILE -d $EXTENDED_LENGTH

##Phytophthora transcripts and the 1000bps downstream

#Move generated extended transcripts file to generated directory
 EXTENDED_TRANSCRIPTS=$(ls -tr | tail -1)
 mv $EXTENDED_TRANSCRIPTS ../phytophthora_data/generated

: <<'END'


#Perform SAGE in silico

#Generate targets of 30bps
IN_FILE=../phytophthora_data/generated/$EXTENDED_TRANSCRIPTS
OUT_FILE=../phytophthora_data/generated/phytophthora_infestans_restricted_regions.fasta
ENZYME_NAME=NlaIII
RECOGNITION_SEQUENCE=CATG
CUT_SITE=4
FRAGMENT_SIZE=30

perl SAGE_insilico.pl -i $IN_FILE -o $OUT_FILE -e $ENZYME_NAME -c $CUT_SITE -r $RECOGNITION_SEQUENCE -f $FRAGMENT_SIZE 

#Create Blast database

#Blast part is missing :(
BLAST_DB=


#Filtering
## Just take blast results with qcoverage equal to lenght of query 
QUERY_COVERAGE = 100
perl filter_blastResults_qcoverage.pl QUERY_COVERAGE blastresultfile

## We could change it! :P Do not make the filtering before creating the sql database
#Targets have lenght of 30bps



END
