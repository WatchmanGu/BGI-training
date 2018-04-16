# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time

#record start time
start = time.clock()

#save data
with open("chr13.scaffold.gff","r") as gff_file:
    list_mRNA_id = []
    dict_mRNA_length = {}
    dict_gene_CDS = {}
    dict_gene_exon = {}
    dict_gene_intro = {}
    for line in gff_file:
        gff_fields = line.strip().split("\t")
        if gff_fields[2] == "mRNA":
            mRNA_id = gff_fields[8].strip().split(";")[0].split("=")[1]
            list_mRNA_id.append(mRNA_id)
            dict_mRNA_length[mRNA_id] = int(gff_fields[4]) + 1 - int(gff_fields[3])
            dict_gene_CDS[mRNA_id] = 0
            dict_gene_exon[mRNA_id] = []
            dict_gene_intro[mRNA_id] = []
            gene_exon_beg = 0
            gene_exon_end = int(gff_fields[4])
        if gff_fields[2] == "CDS":
            dict_gene_CDS[mRNA_id] = dict_gene_CDS[mRNA_id] + int(gff_fields[4]) + 1 - int(gff_fields[3])
            dict_gene_exon[mRNA_id].append(str(int(gff_fields[4]) + 1 - int(gff_fields[3])))
            gene_exon_beg = int(gff_fields[3])
            gene_intro_length = gene_exon_beg - 1 - gene_exon_end
            if gene_intro_length > 0:
                dict_gene_intro[mRNA_id].append(str(gene_intro_length))
            gene_exon_end = int(gff_fields[4])

#Output as a txt file: gene_id mRNA_length CDS_length len_of_exons len_of_intros number_of_exons
with open("advanced_statistic.txt","w") as advanced_statistic_file:
    for mRNA_id in list_mRNA_id:
        advanced_statistic_file.write(mRNA_id + "\t")
        advanced_statistic_file.write(str(dict_mRNA_length[mRNA_id]) + "\t")
        advanced_statistic_file.write(str(dict_gene_CDS[mRNA_id]) + "\t")
        advanced_statistic_file.write(";".join(dict_gene_exon[mRNA_id]) + "\t")
        advanced_statistic_file.write(";".join(dict_gene_intro[mRNA_id]) + "\t")
        advanced_statistic_file.write(str(len(dict_gene_exon[mRNA_id])) + "\n")

"""
Draw a figure to show the result.
Step 1: Count the number.
"""
with open("advanced_statistic.txt","r") as advanced_statistic_file:
    list_mRNA_length = []
    list_CDS_length = []
    list_exon_length = []
    list_intro_length = []
    list_exon_num = []
    for line in advanced_statistic_file:
        gene_attrs = line.strip().split("\t")
        list_mRNA_length.append(int(gene_attrs[1]))
        list_CDS_length.append(int(gene_attrs[2]))
        for exon_length in gene_attrs[3].split(";"):
            list_exon_length.append(int(exon_length))
        for intro_length in gene_attrs[4].split(";"):
            if intro_length != "":    
                list_intro_length.append(int(intro_length))
        list_exon_num.append(int(gene_attrs[5]))

"""
Draw a figure to show the result.
Step 2: Plot the Histogram.
"""
#mRNA_length
list_mRNA_length = sorted(list_mRNA_length)
max_length = list_mRNA_length[-1] - list_mRNA_length[-1]%100 + 200
bin_list = range(0,max_length,100)
plt.hist(list_mRNA_length,bins=bin_list,facecolor="#FF4500")
plt.xlabel("bases/bp")
plt.ylabel("number")
plt.title("Distribution of mRNA Length")
plt.savefig("Histogram of mRNA Length.png",dpi=3000,bbbox_inches="tight")
plt.close()
#CDS_length
list_CDS_length = sorted(list_CDS_length)
max_length = list_CDS_length[-1] - list_CDS_length[-1]%50 + 100
bin_list = range(0,max_length,50)
plt.hist(list_CDS_length,bins=bin_list,facecolor="#006400")
plt.xlabel("bases/bp")
plt.ylabel("number")
plt.title("Distribution of CDS Length")
plt.savefig("Histogram of CDS Length.png",dpi=3000,bbbox_inches="tight")
plt.close()
#exon_length
list_exon_length = sorted(list_exon_length)
max_length = list_exon_length[-1] - list_exon_length[-1]%10 + 20
bin_list = range(0,max_length,10)
plt.hist(list_exon_length,bin_list,facecolor="#48D1CC")
plt.xlabel("bases/bp")
plt.ylabel("number")
plt.title("Distribution of exon Length")
plt.savefig("Histogram of exon Length.png",dpi=3000,bbbox_inches="tight")
plt.close()
#intro_length
list_intro_length = sorted(list_intro_length)
max_length = list_intro_length[-1] - list_intro_length[-1]%10 +20
bin_list = range(0,max_length,10)
plt.hist(list_intro_length,bin_list,facecolor="#4B0082")
plt.xlabel("bases/bp")
plt.ylabel("number")
plt.title("Distribution of intro Length")
plt.savefig("Histogram of intro Length.png",dpi=3000,bbbox_inches="tight")
plt.close()
#exon_num
list_exon_num = sorted(list_exon_num)
max_num = list_exon_num[-1] + 2
bin_list = range(0,max_num,1)
plt.hist(list_exon_num,bin_list,facecolor="#708090")
plt.xlabel("numbers per mRNA")
plt.ylabel("number")
plt.title("Distribution of exon numbers")
plt.savefig("Histogram of exon numbers.png",dpi=3000,bbbox_inches="tight")
plt.close()

#record end time, and show running time
end = time.clock()
print("Running time:%s seconds"%(end - start))
print("Save five figures successfully! ")