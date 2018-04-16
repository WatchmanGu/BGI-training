# -*- coding: utf-8 -*-
import textwrap
#Condon table
condon_table = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L",
                "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
                "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
                "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
                "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
                "CCT":"P","CCC":"P","CCA":"p","CCG":"P",
                "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
                "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
                "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
                "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
                "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
                "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
                "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
                "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
                "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
                "GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
#Open the CDS.fa --> gene ids with their seqs.
with open("cds.fa","r") as cds_file:
    gene_ids = []
    dict_gene_seqs = {}
    for line in cds_file:
        if line[0] == ">":
            gene_id = line.strip().strip(">")
            if gene_id not in gene_ids:
                gene_ids.append(gene_id)
                gene_seq = ""
        gene_seq = gene_seq + line.strip()
        dict_gene_seqs[gene_id] = gene_seq
#Translate DNAs
with open("protein.fa","w") as protein_file:
    for gene_id, gene_seq in dict_gene_seqs.items():
        codon_rang = range(0,len(gene_seq),3)
        protein = ""
        for site in codon_rang:
            codon = gene_seq[site:site+3]
            if codon in condon_table:
                amino_acid = condon_table[codon]
            else:
                amino_acid = "X"
            protein = protein + amino_acid
        protein_file.write(">" + gene_id + "\n")
        protein_lines = textwrap.wrap(protein,width=50)
        for protein_line in protein_lines:
            protein_file.write(protein_line + "\n")