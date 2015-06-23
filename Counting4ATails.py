# Lala
import os, re
from Bio import SeqIO
os.chdir("/Users/eliana/Documents/SAGE_analysis")
handle = open("./phytophthora_data/phytophthora_infestans_t30-4_1_supercontigsSAGE_Targets_1000.fasta", "rU")

frequences = open("./phytophthora_data/count_of_cutting_sequence_in_extended_part.csv", "w")
indexes = open("./phytophthora_data/indexes_of_cutting_sequence_in_extended_part.csv", "w")

for record in SeqIO.parse(handle, "fasta") :
    # Transcript info
    #print record.id
    #print len(record.seq)
    #print record.seq
    
    # Sequence for cutting
    seq = "CATG"
    
    # Taking just the extended part (1000 nucleotides)
    description = str(record.description)
    searchObj = re.search( r'transcript_len:(.*)\|extend_len', description, re.M|re.I)
    extended = int(searchObj.group(1))
    #print extended   
    extended = record.seq[extended:]
    #print extended
    #print len(extended)

    #print extended.find(seq)
    #print record.seq.find(seq)
    index = [m.start() for m in re.finditer(seq, str(extended))]
    #print indexes
    #print len(indexes)
    
    count =  str(record.id) +" ,"+ str(len(extended)) + ", "+ str(len(index))
    print count
    frequences.write(count + "\n")
    print str(record.id) + ", " + str(len(extended)) + ", " + str(index)[1:-1]
    

    #break
     
handle.close()
frequences.close()
indexes.close()