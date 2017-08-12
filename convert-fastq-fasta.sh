#!/bin/bash
# Convert fastq files to fasta with seqtk. 
mkdir ../clado-data-pear-fasta
for name in $(cat ./names-tube.txt)
do
    seqtk seq -a ../clado-data-pear/pear-${name}.assembled.fastq > ../clado-data-pear-fasta/pear-${name}.assembled.fasta
done
