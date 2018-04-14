# -*- coding: utf-8 -*-
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
            dict_gene_CDS[mRNA_id] = dict_gene_CDS[mRNA_id] + int(gff_fields[4] + 1 - int(gff_fields[3]))
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
        advanced_statistic_file.write(str(dict_gene_CDS[mRNA_id]) + "\t")
        advanced_statistic_file.write(";".join(dict_gene_exon[mRNA_id]) + "\t")
        advanced_statistic_file.write(";".join(dict_gene_intro[mRNA_id]) + "\t")
        advanced_statistic_file.write(str(len(dict_gene_exon[mRNA_id]) + "\n"))