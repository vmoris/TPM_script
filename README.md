# TPM_script
Python script to calculate TPM from count file (provided by ht-seq counts)
#for each sample run the first python script to normalize the counts by the lenght of the gene
python Counts_normalized_by_length_only_normalized_counts.py Adineta_vaga.gff3 counts_file.tabular

#run the second python script to calculate the TPM for each gene and each sample
python /media/vmoris/Elements/roti/script/python/TPM_new3.py --output TPM_Rob1_genome2019.tsv --input Counts_normalized_by_lentgh_*
