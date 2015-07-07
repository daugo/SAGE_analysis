#!/usr/bin/env bash

#Compute Transcripts with extended 3' regions if possible up to EXTENDED_LENGTH
SUPERCONTIGS_FASTA=phytophthora_data/original/phytophthora_infestans_t30-4_1_supercontigs.fasta
TRANSCRIPTS_FASTA=phytophthora_data/original/phytophthora_infestans_t30-4_1_transcripts.fasta
ANNOTATION_FILE=phytophthora_data/original/phytophthora_infestans_t30-4_1_transcripts.gtf
EXTENDED_LENGTH=1000

perl Extract_targetRegions.pl -g $SUPERCONTIGS_FASTA -t $TRANSCRIPTS_FASTA -f $ANNOTATION_FILE -d $EXTENDED_LENGTH
#move generated extended transcripts file to generated directory
 EXTENDED_TRANSCRIPTS=$(ls -tr | tail -1)
 mv $EXTENDED_TRANSCRIPTS phytophthora_data/generated

#Perform SAGE in silico
IN_FILE=phytophthora_data/generated/$EXTENDED_TRANSCRIPTS
OUT_FILE=phytophthora_data/generated/AA
ENZYME_NAME=NlaIII
RECOGNITION_SEQUENCE=CATG
CUT_SITE=4
FRAGMENT_SIZE=30

 perl SAGE_insilico.pl -i $IN_FILE -o $OUT_FILE -e $ENZYME_NAME -c $CUT_SITE -r $RECOGNITION_SEQUENCE -f $FRAGMENT_SIZE 
