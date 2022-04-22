#!/bin/bash -l

echo "PFAM and CAZyDB annotation"

#Usage: ./run_pfam_CAZyDB_annotation.sh > run_pfam_CAZyDB_annotation.sh.log 2>&1


###########################################################

cluster_functional_classification_main(){
   check_and_install
   run_clustering
}

# check and install missing packages for pfam pipeline
check_and_install(){
   echo "Checking and installing packages"
   install_cdhit
   install_hmmer
   install_pfam
   install_CAZyDB
   echo "DONE checking and installing packages"
}

install_cdhit(){
   echo "Checking and installing CDHIT"
   if [ -d "${TOOLS_FOLDER}/cdhit" ]; then
      echo "CDHIT already installed"
   else
      #Installing cdhit
      echo "Installing CDHIT"
      #################
      # CDHIT
      cd ${TOOLS_FOLDER}
      wget https://github.com/weizhongli/cdhit/archive/refs/tags/V4.8.1.tar.gz
      tar -zxvf V4.8.1.tar.gz
      cd cdhit-4.8.1
      make
   fi
   echo "DONE checking and installing cdhit-4!"
}

install_hmmer(){
   echo "Checking and installing hmmer"
   if [ -d "${TOOLS_FOLDER}/hmmer" ]; then
      echo "hmmer already installed"
   else
      #Installing hmmer
      echo "Installing hmmer"
      #################
      # hmmer
      cd ${TOOLS_FOLDER}
      wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz
      tar -zxvf hmmer-3.3.2.tar.gz
      cd hmmer-3.3.2
      make
      make check
      make install
   fi
   echo "DONE checking and installing hmmer!"
}

install_pfam(){
   echo "Checking and installing pfam"
   if [ -d "${TOOLS_FOLDER}/pfam" ]; then
      echo "pfam already installed"
   else
      #Installing pfam 
      echo "Installing pfam"
      #################
      # #1.4 Domain/profiles homology searches
      # #1.4.1 Pfam database
      #pfam
      cd ${TOOLS_FOLDER}
      wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz 
      gunzip Pfam-A.hmm.gz
      hmmpress Pfam-A.hmm
   fi
   echo "DONE checking and installing pfam!"
}

install_CAZyDB(){
   echo "Checking and installing CAZyDB"
   if [ -d "${TOOLS_FOLDER}/CAZyDB" ]; then
      echo "CAZyDB already installed"
   else
      #Installing CAZyDB
      echo "Installing CAZyDB"
      #CAZyDB ; build hmm data base of CAZyDB
      cd ${TOOLS_FOLDER}
      wget https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/CAZyDB.09242021.fa
      wget https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/dbCAN-fam-HMMs.txt
      #https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/readme.txt
      hmmpress dbCAN-fam-HMMs.txt
   fi
   echo "DONE checking and installing CAZyDB!"
}


run_clustering(){
   echo "Running clustering"

###################

#CDHIT=${TOOLS_FOLDER}/cdhit/cdhit-4.8.1/

PFAM=${TOOLS_FOLDER}/pfam/
CAZyDB=${TOOLS_FOLDER}/CAZyDB/

workdir=${ANALYSIS_FOLDER}/analysis/Gene_based_analysis_onContigs/functional_classification/prokka_output/
#workdir=${ANALYSIS_FOLDER}/analysis/bin_functional_annotation/prokka_out/

#cd ${workdir}
###################

find ${workdir} -name "*.fna" |sort|while read READ_FW;
do

#While using Metagenomic_all/analysis/Gene_based_analysis_onContigs/functional_classification/prokka_output/ directory use below line,
#SN=$(basename $READ_FW | cut -d. -f1);
#While using binning, Metagenomic_all/analysis/bin_functional_annotation/prokka_out/metabat/ 
#or Metagenomic_all/analysis/bin_functional_annotation/prokka_out/metabat/ directories use below line,
SN=$(basename $READ_FW | cut -d. -f1,2);

DIRNM=$(dirname $READ_FW );
FSN=${DIRNM}"/"${SN}

echo ${SN}
#echo ${DIRNM}
echo ${FSN}

###########################################################
# CD-HIT 
# set environment path to ${CDHIT}

cd-hit -i ${FSN}.faa -o ${FSN}.cdhit.95.fa -c 0.95 -T 70 -M 0

###########################################################

#### Search pfam DB ## HHSEARCH 
hmmsearch --cpu 70 -E 1e-05 --domE 0.35 --domtblout ${FSN}.cdhit.95_PFAM.out ${PFAM}Pfam-A.hmm ${FSN}.cdhit.95.fa

###################
######
# filter PFAM output file
# Cluster all PfamIds, GeneIds, Gene Description based on PROKKA geneIDs

##hmmscan --cpu 5 --domtblout TrinotatePFAM.out ~/RNAseq_assembly_annotation/assembly_annotation/database/trinotate_database/Pfam-A.hmm Trinity_200.fa.transdecoder.pep

###### get accession # Cluster 2nd column based on 1st column
#cat ${FSN}.cdhit.95_PFAM.out | sed '1,3d' | sed '/^#/d' | awk -F" " '{print $4"\t"$2}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g > ${FSN}.cdhit.95_PFAM_clust_acc.out
###### get target name # Cluster 2nd column based on 1st column
#cat ${FSN}.cdhit.95_PFAM.out | sed '1,3d' | sed '/^#/d' | awk -F" " '{print $4"\t"$1}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g > ${FSN}.cdhit.95_PFAM_clust_t-name.out
###### get target name # Cluster 2nd column based on 1st column
#cat ${FSN}.cdhit.95_PFAM.out | sed '1,3d' | sed '/^#/d' | awk '{for (i=23;i<=NF;i++) printf("%s ",$i)} {print """\t"$4}' | awk -F "\t" '{print $NF"\t"$1}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g > ${FSN}.cdhit.95_PFAM_clust_acc-desc.out

######

###### get accession # Cluster 2nd column based on 1st column
cat ${FSN}.cdhit.95_PFAM.out | sed '1,3d' | sed '/^#/d' | awk -F" " '{print $1"\t"$5}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g > ${FSN}.cdhit.95_PFAM_clust_acc.out
###### get target name # Cluster 2nd column based on 1st column
cat ${FSN}.cdhit.95_PFAM.out | sed '1,3d' | sed '/^#/d' | awk -F" " '{print $1"\t"$4}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g > ${FSN}.cdhit.95_PFAM_clust_t-name.out
###### get target name # Cluster 2nd column based on 1st column
cat ${FSN}.cdhit.95_PFAM.out | sed '1,3d' | sed '/^#/d' | awk '{for (i=23;i<=NF;i++) printf("%s ",$i)} {print """\t"$1}' | awk -F "\t" '{print $NF"\t"$1}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g > ${FSN}.cdhit.95_PFAM_clust_acc-desc.out

####################

###########################################################
paste <(cut -f1,2 ${FSN}.cdhit.95_PFAM_clust_acc.out) <(cut -f2 ${FSN}.cdhit.95_PFAM_clust_t-name.out) <(cut -f2 ${FSN}.cdhit.95_PFAM_clust_acc-desc.out) | sed 's/ ;/;/g' | sed 's/  ;/;/g' > ${FSN}.cdhit.95_PFAM_clust_join.out

###### clean files

rm ${FSN}.cdhit.95_PFAM_clust_acc.out
rm ${FSN}.cdhit.95_PFAM_clust_t-name.out
rm ${FSN}.cdhit.95_PFAM_clust_acc-desc.out

######
# build hmm data base of CAZyDB
#https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/readme.txt

#hmmpress dbCAN-fam-HMMs.txt

hmmsearch --cpu 70 -E 1e-05 --domE 0.35 --domtblout ${FSN}.cdhit.95.dbCAN.out.dm ${CAZyDB}dbCAN-fam-HMMs.txt ${FSN}.cdhit.95.fa

sh ${CAZyDB}hmmscan-parser.sh ${FSN}.cdhit.95.dbCAN.out.dm > ${FSN}.cdhit.95.dbCAN.out.dm.ps
cat ${FSN}.cdhit.95.dbCAN.out.dm.ps | awk '$5<1e-15&&$10>0.35' > ${FSN}.cdhit.95.dbCAN.out.dm.stringent 

###### 
####################
# Cluster all PfamIds, GeneIds, Gene Description based on PROKKA geneIDs

###### get accession # Cluster 2nd column based on 1st column
cat ${FSN}.cdhit.95.dbCAN.out.dm | sed '1,3d' | sed '/^#/d' | awk -F" " '{print $1"\t"$5}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g > ${FSN}.cdhit.95.dbCAN_clust_acc.out
###### get target name # Cluster 2nd column based on 1st column
cat ${FSN}.cdhit.95.dbCAN.out.dm | sed '1,3d' | sed '/^#/d' | awk -F" " '{print $1"\t"$4}' | sed s'/.hmm//'g | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g > ${FSN}.cdhit.95.dbCAN_clust_t-name.out
###### get target name # Cluster 2nd column based on 1st column
cat ${FSN}.cdhit.95.dbCAN.out.dm | sed '1,3d' | sed '/^#/d' | awk '{for (i=23;i<=NF;i++) printf("%s ",$i)} {print """\t"$1}' | awk -F "\t" '{print $NF"\t"$1}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g > ${FSN}.cdhit.95.dbCAN_clust_acc-desc.out

####################

###########################################################
paste <(cut -f1,2 ${FSN}.cdhit.95.dbCAN_clust_acc.out) <(cut -f2 ${FSN}.cdhit.95.dbCAN_clust_t-name.out) <(cut -f2 ${FSN}.cdhit.95.dbCAN_clust_acc-desc.out) | sed 's/ ;/;/g' | sed 's/  ;/;/g' > ${FSN}.cdhit.95.dbCAN_clust_join.out


###################
rm ${FSN}.cdhit.95.dbCAN_clust_acc.out
rm ${FSN}.cdhit.95.dbCAN_clust_t-name.out
rm ${FSN}.cdhit.95.dbCAN_clust_acc-desc.out

######
###################
# Stats 

######
# PFAM
#cat ${FSN}.cdhit.95_PFAM.out | sed '1,3d' | sed '/^#/d' | awk -F" " '{print $4"\t"$1}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g | awk -F"; " '{ print $0"\t"NF }' | sort -t $'\t' -n -k3 -r | cut -d$'\t' -f1,3
#
cat ${FSN}.cdhit.95_PFAM.out | sed '1,3d' | sed '/^#/d' | awk '{for (i=23;i<=NF;i++) printf("%s ",$i)} {print """\t"$1}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g | awk -F"; " '{ print $0"\t"NF }' | sort -t $'\t' -n -k3 -r > ${FSN}.cdhit.95_PFAM_clust_target-description.tsv

######
# CAZyDB - dbCAN
#cat ${FSN}.cdhit.95.dbCAN.out.dm | sed '1,3d' | sed '/^#/d' | awk -F" " '{print $4"\t"$1}' | sed s'/.hmm//'g | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g | awk -F"; " '{ print $0"\t"NF }' | sort -t $'\t' -n -k3 -r | cut -d$'\t' -f1,3
#
cat ${FSN}.cdhit.95.dbCAN.out.dm | sed '1,3d' | sed '/^#/d' | awk '{for (i=23;i<=NF;i++) printf("%s ",$i)} {print """\t"$1}' | awk ' !(a[$1]) {a[$1]=$0} a[$1] {w=$1; $1=""; a[w]=a[w] $0} END {for (i in a) print a[i]}' FS="\t" OFS="; " | awk '{$2=""; print $0}' FS=";" OFS=";" | sed s'/;;/;/'g | awk -F"; " '{ print $0"\t"NF }' | sort -t $'\t' -n -k3 -r > ${FSN}.cdhit.95.dbCAN_clust_target-description.tsv

done


###################
# USAGE
# run-Pars-annotation.py INPUT=Inputfile DIR=Directory-fullpath
#python run-Pars-annotation.py SWL.tsv /Metagenomic_all/analysis/Gene_based_analysis_onContigs/functional_classification/prokka_output/

# run on all samples for functional_classification
workdir=${ANALYSIS_FOLDER}/analysis/Gene_based_analysis_onContigs/functional_classification/prokka_output/
find ${workdir} -name "*.fna" |sort|while read READ_FW;
do
#While using Metagenomic_all/analysis/Gene_based_analysis_onContigs/functional_classification/prokka_output/ directory use below line,
SN=$(basename $READ_FW | cut -d. -f1);
DIRNM=$(dirname $READ_FW );
echo ${SN}
echo ${DIRNM}

python run-Pars-annotation.py ${SN}.tsv ${DIRNM}"/"

done

###### Optional usage #######
## run on all samples for bin_functional_annotation
#workdir=${ANALYSIS_FOLDER}/analysis/bin_functional_annotation/prokka_out/
#find ${workdir} -name "*.fna" |sort|while read READ_FW;
#do
##While using binning, Metagenomic_all/analysis/bin_functional_annotation/prokka_out/metabat/ directories use below line,
#SN=$(basename $READ_FW | cut -d. -f1,2);
#DIRNM=$(dirname $READ_FW );
#echo ${SN}
#echo ${DIRNM}
#python run-Pars-annotation.py ${SN}.tsv ${DIRNM}"/"
#done

###### 
echo "DONE running clustering!"

}

bin_functional_classification_main

