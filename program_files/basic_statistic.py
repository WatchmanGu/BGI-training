# -*- coding: utf-8 -*-
import re
import functools
import time

#record start time
start = time.clock()

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

#count Total Length, Effective Length, N Length, GC rate(%,GC_number/Total Length)
total_length = 0
effective_length = 0
N_length = 0
gc_number = 0 
for scaffold_id, scaffold_sequence in dict_scaffolds.items():
    total_length = total_length + len(scaffold_sequence)
for contig_id, contig_sequence in dict_contigs.items():
    effective_length = effective_length + len(contig_sequence)
    gc_number = gc_number + contig_sequence.count("G") + contig_sequence.count("C")
N_length = total_length - effective_length

gc_rate = round(100 * gc_number/effective_length, 2)


"""
count contig N50--> contig_N50_length, scaffold N50--> scaffold_N50_length
build cmp functions for sorted()-->arrange sequences by length from short to long
"""
def scaffold_length_cmp(scaffold1,scaffold2):
    if len(dict_scaffolds[scaffold1]) < len(dict_scaffolds[scaffold2]):
        return 1
    if len(dict_scaffolds[scaffold1]) > len(dict_scaffolds[scaffold2]):
        return -1
    return 0
def contig_length_cmp(contig1,contig2):
    if len(dict_contigs[contig1]) < len(dict_contigs[contig2]):
        return 1
    if len(dict_contigs[contig1]) > len(dict_contigs[contig2]):
        return -1
    return 0
sorted_scaffold_names = sorted(scaffold_names,key = functools.cmp_to_key(scaffold_length_cmp))
sorted_contig_names = sorted(contig_names,key = functools.cmp_to_key(contig_length_cmp))
scaffold_N50_length = 0
contig_N50_length = 0
for scaffold_id in sorted_scaffold_names:
    while scaffold_N50_length < total_length/2:
        scaffold_N50_length = scaffold_N50_length + len(dict_scaffolds[scaffold_id])
    scaffold_N50_length = len(dict_scaffolds[scaffold_id])
    break    
for contig_id in sorted_contig_names:
    while contig_N50_length < effective_length/2:
        contig_N50_length = contig_N50_length + len(dict_contigs[contig_id])
    contig_N50_length = len(dict_contigs[contig_id])
    break
"""
show results
save dict_scaffolds & dict_contigs as file
"""
with open("scaffold_basic_information.txt","w") as results_file:
    print("{0:16}{1:16}{2:16}{3:16}{4:16}{5:16}{6:16}".format("Name","Total Length","Effective Length","N Length","Scaffold N50","Contig N50","GC rate(%)"))
    print("{0:16}{1:16}{2:16}{3:16}{4:16}{5:16}{6:16}".format("scaffold.fa",str(total_length),str(effective_length),str(N_length),str(scaffold_N50_length),str(contig_N50_length),str(gc_rate)))
    result_title = "\t".join(["Name","Total Length","Effective Length","N Length","Scaffold N50","Contig N50","GC rate(%)"]) + "\n"
    results = "\t".join(["scaffold.fa",str(total_length),str(effective_length),str(N_length),str(scaffold_N50_length),str(contig_N50_length),str(gc_rate)]) + "\n"
    results_file.write(result_title)
    results_file.write(results)

"""
with open("dict_scaffolds.txt","w") as dict_scaffold_file:
    for scaffold_id, scaffold_sequence in dict_scaffolds.items():
        dict_scaffold_file.write(scaffold_id + "\t" + scaffold_sequence + "\n")
with open("dict_contigs.txt","w") as dict_contigs_file:
    for contig_id, contig_sequence in dict_contigs.items():
        dict_contigs_file.write(contig_id + "\t" + contig_sequence + "\n")
"""

#record end time, and show running time
end = time.clock()
print("Running time:%s seconds"%(end - start))
print("Txt file here --> scaffold_basic_information.txt")