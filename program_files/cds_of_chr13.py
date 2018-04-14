# -*- coding: utf-8 -*-
import textwrap
def sequence_reverse(string):
    """
    A function to reverse a sequence.
    e.g. ATCG --> GCTA.
    Remind: '+', positive, 5'->3';'-', negative, 3'->5'.
    """
    return string[::-1]
#Chromosome sequence --> a string like "ATCCGAT..."
with open("chr13.fa","r") as chromosome_file:
    chromosome_seq = ""
    for line in chromosome_file:
        if line[0] != "<":
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
            cds_beg = int(gff_fields[3]) - 1#- 1 for list index 0-x
            cds_end = int(gff_fields[4]) - 1
            cds_seq = cds_seq + chromosome_seq[cds_beg:cds_end + 1]
#how to with frame 1,2,3?!!!!!!
            if gff_fields[6] == "+":
                dict_cds[cds_id] = cds_seq
            if gff_fields[6] == "-":
                dict_cds[cds_id] = sequence_reverse(cds_seq)
#For checking if the len(CDS)%3 not equal to 0
for cds_id in cds_gene_ids:
    if len(dict_cds[cds_id])%3 != 0:
        print(cds_id)
#Output to cds.fa     
with open("cds.fa","w") as cds_seq_file:
    for cds_id, cds_seq in dict_cds.items():
        cds_seq_file.write(">" + cds_id + "\n")
        cds_seq_lines = textwrap.wrap(cds_seq, width = 50)
        for cds_seq_line in cds_seq_lines:
            cds_seq_file.write(cds_seq_line + "\n")




