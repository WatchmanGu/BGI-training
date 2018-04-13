# -*- coding: utf-8 -*-

#build a dictionary for chr13.agp: dict_scaffold_location_on_chr13-->{scaffold_id:location_beg,...}
dict_scaffold_location_on_chr13 = {}
with open("chr13.agp") as chr13_agp_file:
    for line in chr13_agp_file:
        agp_fields = line.strip().split('\t')
        scaffold_id = agp_fields[0]
        chromosome_beg = agp_fields[1]
        dict_scaffold_location_on_chr13[scaffold_id] = chromosome_beg
"""
Read the chr13.scaffold.gff file and replace the location field with the location on the chromosome.
Than output the results as chromosome_gene.gff.
"""
with open("chr13.scaffold.gff") as chr13_scaffold_gff_file:
    for line in chr13_scaffold_gff_file:
        gff_fields = line.strip().split('\t')
        if gff_fields[2] == "CDS":
            scaffold_id = gff_fields[0]
            scaffold_beg = int(gff_fields[3])
            scaffold_end = int(gff_fields[4])
            chromosome_beg = int(dict_scaffold_location_on_chr13[scaffold_id])
            cds_on_chromosome_beg = chromosome_beg - 1 + scaffold_beg
            cds_on_chromosome_end = chromosome_beg - 1 + scaffold_end
            gff_fields[3] = str(cds_on_chromosome_beg)
            gff_fields[4] = str(cds_on_chromosome_end)
        with open("chromosome_gene.gff","w") as chromosome_gene_anotation_file:
            for field in gff_fields:
                chromosome_gene_anotation_file.write(field)
                chromosome_gene_anotation_file.write("\t")
            chromosome_gene_anotation_file.write("\n")



