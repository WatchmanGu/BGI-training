# -*- coding: utf-8 -*-
import time 
import sys
#record start time
start = time.clock()

#build a dictionary for chr13.agp: dict_scaffold_location_on_chr13-->{scaffold_id:location_beg,...}
dict_scaffold_beg_on_chr13 = {}
dict_scaffold_end_on_chr13 = {}
dict_scaffold_direction_on_chr13 = {}
with open("chr13.agp") as chr13_agp_file:
    for line in chr13_agp_file:
        agp_fields = line.strip().split('\t')
        if agp_fields[5] != "100":
            scaffold_id = agp_fields[5]
            chromosome_beg = agp_fields[1]
            chromosome_end = agp_fields[2]
            dict_scaffold_beg_on_chr13[scaffold_id] = chromosome_beg
            dict_scaffold_end_on_chr13[scaffold_id] = chromosome_end
            dict_scaffold_direction_on_chr13[scaffold_id] = agp_fields[8]

"""
Read the chr13.scaffold.gff file and replace the location field with the location on the chromosome.
Than output the results as chromosome_gene.gff.
"""
with open("chromosome_gene.gff","w") as chromosome_gene_anotation_file:
    #Clean up the output file if it already exit.
    chromosome_gene_anotation_file.write("")
with open("chr13.scaffold.gff") as chr13_scaffold_gff_file:
    negative_negative_ids = []
    positive_positive_ids = []
    positive_negative_ids = []
    negative_positive_ids = []
    for line in chr13_scaffold_gff_file:
        gff_fields = line.strip().split('\t')
        scaffold_id = gff_fields[0]
        cds_id = gff_fields[8].strip().split(";")[0].split("=")[1]
        scaffold_beg = int(gff_fields[3])
        scaffold_end = int(gff_fields[4])
        if dict_scaffold_direction_on_chr13[scaffold_id] == "+" or dict_scaffold_direction_on_chr13[scaffold_id] == "n":
            chromosome_beg = int(dict_scaffold_beg_on_chr13[scaffold_id])         
            cds_on_chromosome_beg = chromosome_beg - 1 + scaffold_beg
            cds_on_chromosome_end = chromosome_beg - 1 + scaffold_end
            gff_fields[3] = str(cds_on_chromosome_beg)
            gff_fields[4] = str(cds_on_chromosome_end)
            if gff_fields[6] == "-":
                gff_fields[6] = "-"
                if cds_id not in positive_negative_ids:
                    positive_negative_ids.append(cds_id)
            if gff_fields[6] == "+":
                gff_fields[6] = "+"
                if cds_id not in positive_positive_ids:
                    positive_positive_ids.append(cds_id)
        if dict_scaffold_direction_on_chr13[scaffold_id] == "-":  
            chromosome_end = int(dict_scaffold_end_on_chr13[scaffold_id])          
            cds_on_chromosome_beg = chromosome_end + 1 - scaffold_end
            cds_on_chromosome_end = chromosome_end + 1 - scaffold_beg
            gff_fields[3] = str(cds_on_chromosome_beg)
            gff_fields[4] = str(cds_on_chromosome_end)
            if gff_fields[6] == "-":
                gff_fields[6] = "++"
                if cds_id not in negative_negative_ids:
                    negative_negative_ids.append(cds_id)
            if gff_fields[6] == "+":
                gff_fields[6] = "--"
                if cds_id not in negative_positive_ids:
                    negative_positive_ids.append(cds_id)
            gff_fields[6] = gff_fields[6][0]
        gff_fields[0] = "Pa13"
        
        with open("chromosome_gene.gff","a") as chromosome_gene_anotation_file:
            for field in gff_fields:
                chromosome_gene_anotation_file.write(field)
                chromosome_gene_anotation_file.write("\t")
            chromosome_gene_anotation_file.write("\n")

#record end time, and show running time
end = time.clock()
print("Running time:%s seconds"%(end - start))
print("Txt file here --> chromosome_gene.gff")