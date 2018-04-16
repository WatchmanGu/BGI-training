# -*- coding: utf-8 -*-
import textwrap
import time


#record start time
start = time.clock()

#function to reverse sequence
def sequence_reverse(string):
    """
    A function to reverse a sequence.
    e.g. ATCG --> GCTA.
    Remind: '+', positive, 5'->3';'-', negative, 3'->5'.
    """
    return string[::-1]

#function to get the complementary sequence
def complementary_sequence(string):
    """
    a function to get the complementary sequence.
    e.g. ATCG --> TAGC
    """
    return string.replace("A","t").replace("T","a").replace("C","g").replace("G","c").upper()
    
#Chromosome sequence --> a string like "ATCCGAT..."
with open("chr13.fa","r") as chromosome_file:
    chromosome_seq = ""
    for line in chromosome_file:
        if line[0] != ">":
            chromosome_seq = chromosome_seq + line.strip()

"""
Read chromosome_gene.gff file and map the CDS to the chromosome
Than output as a file, including the gene_id and the sequence
"""
with open("chromosome_gene.gff","r") as chromosome_anotation_file:
    cds_gene_ids = []
    dict_cds = {}
    for line in chromosome_anotation_file:
        gff_fields = line.strip("\n").split("\t")
        if gff_fields[2] == "CDS":
            cds_id = gff_fields[8].strip(";").split("=")[1]#delete "parent"
            if cds_id not in cds_gene_ids:
                cds_gene_ids.append(cds_id)
                cds_seq = ""
                cds_end = 0
            cds_beg = int(gff_fields[3]) - 1#- 1 for list index 0-x
            if cds_beg >= cds_end:
                cds_end = int(gff_fields[4])
                cds_seq = cds_seq + chromosome_seq[cds_beg:cds_end]
            else:
                cds_end = int(gff_fields[4])
                cds_seq = chromosome_seq[cds_beg:cds_end] + cds_seq
            if gff_fields[6] == "+":
                dict_cds[cds_id] = cds_seq
            if gff_fields[6] == "-":
                dict_cds[cds_id] = complementary_sequence(sequence_reverse(cds_seq))
#For checking if the len(CDS)%3 not equal to 0
#for cds_id in cds_gene_ids:
    #if len(dict_cds[cds_id])%3 != 0:
        #print("Something wrong with these CDS:",cds_id)
#Output to cds.fa     
with open("cds.fa","w") as cds_seq_file:
    for cds_id, cds_seq in dict_cds.items():
        cds_seq_file.write(">" + cds_id + "\n")
        cds_seq_lines = textwrap.wrap(cds_seq, width = 50)
        for cds_seq_line in cds_seq_lines:
            cds_seq_file.write(cds_seq_line + "\n")



#record end time, and show running time
end = time.clock()
print("Running time:%s seconds"%(end - start))
print("Results --> cds.fa")

