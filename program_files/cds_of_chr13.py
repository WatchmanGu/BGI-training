# -*- coding: utf-8 -*-

def sequence_reverse(string):
    """
    A function to reverse a sequence.
    e.g. ATCG --> GCTA.
    Remind: '+', positive, 5'->3';'-', negative, 3'->5'.
    """
    return string[::-1]
"""
Read chromosome_gene.gff file and map the CDS to the chromosome
Than output as a file, including the gene_id and the sequence
"""
#with open