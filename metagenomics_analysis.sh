#!/usr/bin/env bash

### Workflow for whole genome shotgun metagenomics analysis

# Root folder name"
NAME=Metagenomic_all

# Raw data folder path
#SRC_RAWDATA='/data1/Active_Projects/Metagenomic_QC/rawdata/'
#LINKPATH_DB='/data1/Active_Projects/paper_scripts/reference/'

SRC_RAWDATA='/home/gala0002/proj/proj_oliyad-2022/data/F21FTSEUHT0827_MICaviM/Clean/*'
LINKPATH_DB='/home/gala0002/proj/proj_oliyad-2022/references'


metagenomics_analysis_main(){
   create_folders
   set_variables # -> Never comment this function 
   #fetch_example_data # -> Uncomment this function if you want to run it on an example data 
   #copy_rawdata

<<DONTRUN
   # Keep the below scripts here as a comment to avoid execution 
DONTRUN
   run_qc
   run_assembly
   run_coassembly
   run_reference_analysis
   run_comparative_analysis
   run_coverage_and_bining
   run_binrefinement
   run_bin_taxonomic_classification
   run_bin_functional_classification
   run_comparative_analysis
   run_reference_analysis

######

   echo $LINKPATH_DB
}

create_folders(){

   echo "Creating sub-folders..."

   # Sub-folders in the root folder
   for FOLDER in analysis tools rawdata reference bin
   do
      mkdir -p $NAME/$FOLDER
   done
   echo "DONE creating sub-folders!"
}

# setting variable path
set_variables(){
   echo "Setting variables for paths..."

   export ROOT_FOLDER_NAME=$NAME
   export TOOLS_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/tools
   export RAWDATA_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/rawdata
   export ANALYSIS_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/analysis
   export REFERENCE_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/reference
   export BIN_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/bin
   export LINKPATH_DB=$LINKPATH_DB
   echo "DONE setting variables for paths!"
}


# copy raw data from source folder to analysis folder structure 
copy_rawdata(){
   lst=$(ls -d $SRC_RAWDATA/*.fq.gz)
   for file in $lst 
   do
      echo "Copying ${file}"
      cp ${file} ${RAWDATA_FOLDER}/
   done
   echo "DONE copying rawdata!"
}


# fetch raw data from web to the analysis folder structure 
fetch_example_data(){

   mkdir -p $NAME/example_data

   cd $NAME/example_data

   wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_1.fastq.gz
   wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_2.fastq.gz

   SRC_RAWDATA=$NAME/example_data
   cd -
}

# Run QC analysis 
run_qc(){
   echo "Running Quality Control" 

   . ./run_qc.sh

   echo "DONE Quality Control!" 
}

# Run assembly 
run_assembly(){
   echo "Running Assembly" 

   . ./run_assembly.sh

   echo "DONE Running Assembly!" 
}

# Run assembly 
run_coassembly(){
   echo "Running Coassembly" 

   . ./run_coassembly.sh

   echo "DONE Running Coassembly!" 
}

run_reference_analysis(){
   echo "Running Reference Analysis" 

   . ./run_reference_analysis.sh

   echo "DONE Running Reference Analysis!" 
}

run_comparative_analysis(){
   echo "Running Comparative Analysis" 

   . ./run_comparative_analysis.sh

   echo "DONE Running Comparative Analysis!" 
}

run_coverage_and_bining(){
   echo "Running Coverage and Bining" 

   . ./run_coverage_bining.sh

   echo "DONE running Coverage and Bining!" 
}

run_binrefinement(){
   echo "Running Binrefinement" 

   . ./run_binrefinement_v2.sh

   echo "DONE running Binrefinement!" 
}

run_bin_taxonomic_classification(){
   echo "Running Bin Taxonomic Classification" 

   . ./run_bin_taxonomic_classification.sh

   echo "DONE running Bin Taxonomic Classification!" 
}

run_bin_functional_classification(){
   echo "Running Bin Functional Classification" 

   . ./run_bin_functional_classification.sh

   echo "DONE running Bin Functional Classification!" 
}

metagenomics_analysis_main
