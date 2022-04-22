
# Metagenomics_analysis

USAGE:

```
sh metagenomics_analysis.sh > metagenomics_analysis.sh.log 2>&1
```



Workflow from the paper: 
Current Challenges and Best Practice Protocols for Microbiome Analysis using Amplicon and Metagenomic Sequencing

Modified and adapted to the local Linux server environment. Several databases, tools, and URLs have been adopted to the local environment.



### Additional annotation scripts were developed as an extension to the original workflow 

* PFAM and CAZyDB annotation:
* USAGE:

```
bash run_pfam_CAZyDB_annotation.sh > run_pfam_CAZyDB_annotation.sh.log 2>&1

```



### Script to Pars and merge annotations, this script was integrated into the above workflow script, however, you can also run this independently

* USAGE: 
* run-Pars-annotation.py INPUT=Inputfile DIR=Directory-fullpath

```
python run-Pars-annotation.py sample.tsv /prokka_output/

```



Cite:
Richa Bharti, Dominik G Grimm, Current challenges and best-practice protocols for microbiome analysis, Briefings in Bioinformatics, Volume 22, Issue 1, January 2021, Pages 178–193, https://doi.org/10.1093/bib/bbz155
