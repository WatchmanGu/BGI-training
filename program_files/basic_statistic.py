# -*- coding: utf-8 -*-
import re
"""
build a dictionary of scaffold names with their sequences
 -->dict_scaffolds = {"scaffold_id":"ATCG...",....}
"""
with open("Q3_plum_0630.scafSeq60.Fa") as scaffold_file:
    dict_scaffolds = {}
    scaffold_names = []
    for line in scaffold_file:
        if line[0] == ">":
            scaffold_id_line = line.strip().split()
            scaffold_id = scaffold_id_line[0].strip(">")
            scaffold_names.append(scaffold_id)
            scaffold_sequence = ""
        scaffold_sequence = scaffold_sequence + line.strip()
        dict_scaffolds[scaffold_id] = scaffold_sequence
"""
build dictionary of contigs with their sequences
--> dict_contigs = {"contig_id": "ATCG...",...}
"""
remove_N_pattern = re.compile(r"N+")
dict_contigs = {}
contig_id_number = 1
contig_names = []
for scaffold_name in scaffold_names:
    raw_sequence = dict_scaffolds[scaffold_name]
    raw_contigs = re.split(remove_N_pattern,raw_sequence)
    for contig in raw_contigs:
        contig_id = "contig" + str(contig_id_number)
        contig_names.append(contig_id)
        contig_id_number = contig_id_number + 1
        dict_contigs[contig_id] = contig
#count Total Length, Effective Lenth, N Length, GC rate(%,GC_number/Total Length)
total_length = 0
effective_lenth = 0
N_length = 0
GC_number = 0 
for scaffold_id, scaffold_sequence in dict_scaffolds.items():
    total_length = total_length + len(scaffold_sequence)
for contig_id, contig_sequence in dict_contigs:
    effective_lenth = effective_lenth + len(contig_sequence)
    for base in contig_sequence:
        if base == "G" or base == "C":
            GC_number = GC_number + 1
N_length = total_length - effective_lenth
#accurate to two decimal places
gc_rate = round(100 * GC_number/total_length, 2)
"""
count contig N50--> contig_N50_length, scaffold N50--> scaffold_N50_length
build cmp functions for sorted()-->arrange sequences by length from short to long
"""
def scaffold_length_cmp(scaffold1,scaffold2):
    if len(dict_scaffolds[scaffold1]) < len(dict_scaffolds[scaffold2]):
        return -1
    if len(dict_scaffolds[scaffold1]) > len(dict_scaffolds[scaffold2]):
        return 1
    return 0
def contig_length_cmp(contig1,contig2):
    if len(dict_contigs[contig1]) < len(dict_contigs[contig2]):
        return -1
    if len(dict_contigs[contig1]) > len(dict_contigs[contig2]):
        return 1
    return 0
sorted_scaffold_names = sorted(scaffold_names,scaffold_length_cmp)
sorted_contig_names = sorted(contig_names,contig_length_cmp)
scaffold_N50_length = 0
contig_N50_length = 0
for scaffold_id in sorted_scaffold_names:
    scaffold_N50_length = scaffold_N50_length + len(dict_scaffolds[scaffold_id])
    if scaffold_N50_length >= total_length/2:
        break
for contig_id in sorted_contig_names:
    contig_N50_length = contig_N50_length + len(dict_contigs[contig_id])
    if contig_N50_length >= effective_lenth/2:
        break
"""
show results
save dict_scaffolds & dict_contigs as file
"""
with open("scaffold_basic_infomation.txt","w") as results_file:
    result_title = "Name\tTotal Length\tEffective Length\tN Length\tScaffold N50\tContig N50\tGC rate(%)\n"
    results = "scaffold.fa" + "\t" + str(total_length) + "\t" + str(effective_lenth) + "\t" + str(N_length) + "\t" + str(scaffold_N50_length) + "\t" + str(contig_N50_length) + "\t" + str(gc_rate) + "\n"
    results_file.write(result_title)
    results_file.write(results)
with open("dict_scaffolds.txt","w") as dict_scaffold_file:
    for scaffold_id, scaffold_sequence in dict_scaffolds.items():
        dict_scaffold_file.write(scaffold_id + "\t" + scaffold_sequence + "\n")
with open("dict_contigs.txt","w") as dict_contigs_file:
    for contig_id, contig_sequence in dict_contigs.items():
        dict_contigs_file.write(contig_id + "\t" + contig_sequence + "\n")        