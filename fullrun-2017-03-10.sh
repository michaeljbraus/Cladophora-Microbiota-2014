#!/bin/bash
# Full analysis of 16S, mock data. 
# RUN ONE COMMAND AT A TIME TO ENSURE THE COMMANDS ARE CORRECT, UNLESS RE-RUNNING. 
###
# First enter in QIIME (1.9)
# Merge all the files and add the sample IDs to each read. 
add_qiime_labels.py -i clado-data-pear-fasta/ -m 16S-Analysis-Tutorial/map.txt -c InputFileName -n 1 -o combined_fasta
# Count the reads, just a quick check; good to do checks like this to see it's working. 
count_seqs.py -i combined_fasta/combined_seqs.fna > log-combined_seqs.txt
# Pick OTUs. View the documentation for pick_otus.py to read about its defaults. 
pick_otus.py -i combined_fasta/combined_seqs.fna -o picked_otus_default
# Picking representative sequences. 
mkdir rep_seqs
pick_rep_set.py -i picked_otus_default/combined_seqs_otus.txt -f combined_fasta/combined_seqs.fna -o rep_seqs/rep_seqs.fasta -l rep_seqs/rep_seqs.log
# Quick count. 
count_seqs.py -i rep_seqs/rep_seqs.fasta
# Align sequences. 
align_seqs.py -i rep_seqs/rep_seqs.fasta -o pynast_aligned/
# Just making some log files of counts; feel free to check them out if you like. 
count_seqs.py -i pynast_aligned/rep_seqs_failures.fasta > pynast_aligned/rep_seqs_failures.log
count_seqs.py -i pynast_aligned/rep_seqs_aligned.fasta > pynast_aligned/rep_seqs_aligned.log
# Now assign taxonomy. Note it's using uclust, which is slow but good. 
mkdir assigned_taxonomy_uclust
assign_taxonomy.py -i rep_seqs/rep_seqs.fasta -t greengenes/gg_13_5_taxonomy.txt -r greengenes/gg_13_5_otus/rep_set/97_otus.fasta -m uclust -o assigned_taxonomy_uclust/
# Another sanity check. 
head rdp_assigned_taxonomy/rep_seqs_tax_assignments.txt
# Make the biom file. 
make_otu_table.py -i picked_otus_default/combined_seqs_otus.txt -o clado_otu_table.biom -e pynast_aligned/rep_seqs_failures.fasta -t assigned_taxonomy_uclust/rep_seqs_tax_assignments.txt
# It's a mess, but it's good to check this worked. 
head clado_otu_table.biom
# Summarize the biom table. 
biom summarize-table -i clado_otu_table.biom -o clado_otu_table_summary.txt