# Cladophora-Microbiota-2014
16S Amplicon Analysis of Cladophora Collected from Lake Mendota in 2014

For a full summary of methods and analysis, see the published paper in the Journal of Phycology. 

Link: http://onlinelibrary.wiley.com/doi/10.1111/jpy.12573/full

---

This paper is based on my master's thesis. It's not as good as the published paper, but it takes a more environmental-studies-ish approach to the study of the microbiota of lake macroalgae. This is much more technical and relevant to the biology of Cladophora glomerata. 

All 16S amplicon reads for this study can be found on the NCBI Short Read Archive (SRA). 
* Link: https://www.ncbi.nlm.nih.gov/sra. (Search for "cladophora glomeratat microbiota") 
* You can see how these were uploaded by looking in "Short-Reads-Archive/"

The UW-Madison Biotechnology Center has provided specific methods as to how these amplicons were produced. 
* See "Construction and Sequencing of 16s libraries - 20161017.doc". 
* Also see "Updated_Illumina_Submission_Form_1.pdf" for reference. 

These data are listed under the NCBI BioProject "Cladophora 16S 2014 Accession: PRJNA360140 ID: 360140". 
* Link: https://www.ncbi.nlm.nih.gov/bioproject/360140

Database used for analysis of amplicon reads for Braus et al. (2017): SILVA 128. 
* Link: https://www.arb-silva.de/download/archive/qiime/. (Download "Silva_128_release.tgz")

Water data from Lake Mendota can be retrieved from the UW-Madison Center for Limnology's LTER database. 
* This study used "north_temperate_lakes_lter__daily_water_temperature_-_lake_mendota.csv" retrieved in 2016 to make the supplementary figure about lake water temperatures. 
* Link: https://lter.limnology.wisc.edu/data

---

Install the following: 
* QIIME. Link: http://qiime.org/install/install.html
* PEAR. Link: https://sco.h-its.org/exelixis/web/software/pear/
* seqtk. Link: https://github.com/lh3/seqtk

PROTOCOL:
* Decompress the reads (run747.zip) and the SILVA database. 
* Run "pairedends.sh". 
* Run "convert-fastq-fasta.sh"
* Run the QIIME commands in "Braus2017-fullrun.txt". 
* Move "clado_otu_table_SILVA.biom" into the working directory. 
* Run "Clado-Re-Analysis.R" to create tables and figures (see "Figs/"). 

---

The reviewers at Journal of Phycology requested I re-analyze the data from Zulkifly et al. (2015). 
* The protocol was the same as above, just one fastq file, and the requited files all have "Zulk2012" in their title. 
* The code currently uses the Greengenes database, but SILVA can be used as well. 
* The reads did not require processing like the Illumina reads for my 2017 study. QIIME makes "otu_table-Zulk2012.biom", and "Zulk2012-Analysis.R" makes the tables and figures (see "Zulk2012-Figs/"). 

---

*Omnis cellula e cellula. -Rudolf Virchow*