#!/bin/sh

#preparing genes_of_interest file
#in a directory where u have the w_q_de_80 file from pannzer, or any file that is a list of different expressed genes.
#usage: sh preparing_interesting_genes.sh <differential_expressed_genes_file>

awk '{print $1}' "$1" > genes_of_interest
sed -i 's/TRINITY_//g' geneID2GO
