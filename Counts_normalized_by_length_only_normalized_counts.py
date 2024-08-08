#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from decimal import *
getcontext().prec = 10

# ARGS ORDER: gff3 tab_file_1 tab_file_2 ... tab_file_n
list_files = sys.argv[2:] # python commence à compter à partir de 0
gff = sys.argv[1] # python commence à compter à partir de 0

# Get for each sample the quantification from the file
dict_of_sample = {}
for f in list_files :
	dict_of_sample[f] = {}
	for line in open(f, "r") :
		s = line.strip().split("\t") # strip() retire le \n & split("\t") crée une liste à partir des colonnes séparées par \t
		dict_of_sample[f][s[0]] = int(s[1]) # s0 = gene ID & s1 = valeur (counts of reads)

# Get CDS length per parent id
cds_length = {}
for line in open(gff, "r") :
	if line[0] == "#" :
		continue
	s = line.strip().split("\t")

	if s[2] != "CDS" :
		continue

	substraction = int(s[4]) - int(s[3])

	# Parse lines like:
	# ID=FUN_000002-T1.exon3;Parent=FUN_000002-T1;
	frmt = s[8].split(";") # split last column containing Parent and ID
	for el in frmt : # For each element in the split column
		sel = el.split("=") # Re-split by "="
		if sel[0] == "Parent" : # If first element of the column == Parent
			if sel[1] not in cds_length : # If Parent value is not already in the dict
				cds_length[sel[1]] = substraction # Create the key value pair
			else : # If it is already in: add substraction
				cds_length[sel[1]] += substraction
"""
for key, value in cds_length.items() :
	print(key, "\t", value)
"""

# counts est un dictionnaire
gene_tpm = {}
for sample, counts in dict_of_sample.items() :
	#print(sample, counts)
	for funid, count in counts.items() :
		up = Decimal(count)
		down = Decimal(cds_length[funid])
		tpm = Decimal(up/down) # CHANGER LA FORMULE ICI
		if funid not in gene_tpm :
			gene_tpm[funid] = {} # On crée un dictionnaire dont les clés sont des funid (gene) et les values sont des dictionnaires dont les clés sont des samples et les values des tpms
		gene_tpm[funid][sample] = tpm

output = "Counts_normalized_by_lentgh_63.tpms"
f = open(output, "w")

for funid, sample_dict in gene_tpm.items() :
	out_string = funid + "\t"
	for sample, tpm in sample_dict.items() :
		out_string += str(tpm) + "\t"
	f.write(out_string + "\n")

f.close()
