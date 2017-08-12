#!/bin/bash
# Merge raw reads with PEAR.
# PEAR v0.9.9 [May 13, 2016]
mkdir ../clado-data-pear
for name in $(cat ./names-tube.txt)
do
    pear -m 500 -f ../clado-data/*${name}*R1.fastq -r ../clado-data/*${name}*R2.fastq -o ../clado-data-pear/pear-${name} -j 10
done

