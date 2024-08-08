#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Modules
import sys
import argparse
import os

from decimal import *
getcontext().prec = 10 # Report 10 digits in the Decimal lib

def main(args) :
	list_files = list(set(args.input)) # List of input files provided after the --input argument

	# Get for each sample the quantification from the file
	dict_of_samples = {}
	for input_file in list_files : # For each file in the input list
		dict_of_samples[input_file] = {}
		for line in open(input_file, "r") : # For each line in the current file
			s = line.strip().split("\t") # strip() retire le \n & split("\t") crée une liste à partir des colonnes séparées par \t
			dict_of_samples[input_file][s[0]] = Decimal((s[1])) # s0 = gene ID & s1 = valeur # Fill a dictionary key:value gene_id:value

	# Total number of transcripts per sample
	total_number_of_transcripts = {}
	for sample, count_dict in dict_of_samples.items() : # sample = key & counts = value (here value is another dictionary)
		total_number_of_transcripts[sample] = sum(count_dict.values()) # sum of values of all genes in the considered sample

	for k, v in total_number_of_transcripts.items() :
		print(k, v)

	# counts est un dictionnaire
	gene_tpm = {}
	for sample, count_dict in dict_of_samples.items() :
		#print(sample, counts)
		for funid, count in count_dict.items() :
			up = Decimal(count*1000000)
			down = Decimal(total_number_of_transcripts[sample])
			tpm = Decimal(up/down) # CHANGER LA FORMULE ICI
			if funid not in gene_tpm.keys() :
				gene_tpm[funid] = {} # On crée un dictionnaire dont les clés sont des funid (gene) et les values sont des dictionnaires dont les clés sont des samples et les values des tpms
			gene_tpm[funid][sample] = tpm

	output = args.output[0]
	f = open(output, "w")
	out = "GENE\t" + "\t".join(sample for sample in dict_of_samples.keys())
	f.write(out + "\n")

	samples_list = [sample for sample in dict_of_samples.keys()]
	for funid, sample_dict in gene_tpm.items() :
		out_string = funid + "\t"
		for sample in samples_list :
			tpm = sample_dict[sample]
			out_string += str(tpm) + "\t"
		f.write(out_string + "\n")

	f.close()

def get_args() :
	"""Argument parser"""
	parser = argparse.ArgumentParser(description='Compute TPMs from a set of genes and files')
	parser.add_argument("-i", "--input", required=True, nargs="+", help="<REQUIRED> a list of files")
	parser.add_argument("-o", "--output", required=True, nargs=1, help="<REQUIRED> a valid output path")
	return parser.parse_args()


if __name__ == "__main__" :
	args = get_args()
	main(args)
