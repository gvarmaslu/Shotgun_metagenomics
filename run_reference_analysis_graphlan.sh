#!/usr/bin/evn bash

### Refereence based Analysis


ref_analysis_main(){
   #check_and_install
   #run_kraken
   run_metaphlan
   #run_diamond
   #run_humann
}

# Check and install missing packages for the QC pipeline
check_and_install(){
   echo "Checking and installing packages for Reference Based Analysis"
   install_kraken
   download_krakendb
   install_metaphlan
   download_metaphlan
   install_megan
   download_megandb
   install_diamond
   download_nrdb
   install_humann
   download_humanndb
   install_graphlan
   echo "DONE checking and installing packages for Reference Based Analysis!"
}


run_kraken(){
   echo "Running kraken2"

   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/kraken
   cd ${ANALYSIS_FOLDER}/reference_classification/kraken

   # running kraken
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/kraken2-2.0.8-beta/kraken2 \
      --db ${REFERENCE_FOLDER}/reference_database/kraken2DB \
      --threads 70 \
      --use-names \
      --paired \
      ${i} \
      ${i%1*}2.unmerged.final.clean.fq \
      --report ${ANALYSIS_FOLDER}/reference_classification/kraken/$(basename ${i%1*}kreport) > \
      ${ANALYSIS_FOLDER}/reference_classification/kraken/$(basename ${i%1*}kraken)
   done

   echo "DONE running kraken2!"
}


run_metaphlan(){
   echo "Running metaphlan"

   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan2
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/profiled_samples
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/output_images
   cd ${ANALYSIS_FOLDER}/reference_classification/metaphlan2

   #<<COM1
   #finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   #echo "Running metaphlan"
   #for i in ${finallist}
   #do
   #   ${TOOLS_FOLDER}/metaphlan2/metaphlan2.py \
   #   --input_type multifastq \
   #   ${i},${i%1*}2.unmerged.final.clean.fq,${i%1*}u.final.clean.fq \
   #   --bowtie2db ${REFERENCE_FOLDER}/reference_database/metaphlan2/db_v20/mpa_v20_m200 \
   #   --bt2_ps very-sensitive \
   #   --nproc 70 \
   #   --mpa_pkl ${REFERENCE_FOLDER}/reference_database/metaphlan2/db_v20/mpa_v20_m200.pkl \
   #   --bowtie2out ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/$(basename ${i%1*}bt2out) > ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/profiled_samples/$(basename ${i%1*}txt)
   #done
   #COM1


   echo "Creating sample profiles using metaphlan"
   profilelist=$(ls -d ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/profiled_samples/* | awk '{print $NF}')

   #Merge the output
   ${TOOLS_FOLDER}/metaphlan2/utils/merge_metaphlan_tables.py \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/profiled_samples/*.txt > \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/merged_abundance_table.txt

   # create a species abundance table
   grep -E "(s__)|(^ID)" ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/merged_abundance_table.txt \
   | grep -v "t__" | sed 's/^.*s__//g' > \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/merged_abundance_table_species.txt

   ########
   #conda create -n graphlan
   #conda activate graphlan
   #conda install -c biobakery graphlan
   #conda install -c bioconda export2graphlan
   #conda install -c anaconda biopython
   ########

   echo "Creating heatmap"
   profilelist=$(ls -d ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/profiled_samples/* | awk '{print $NF}')
   ${TOOLS_FOLDER}/hclust2/hclust2.py \
   -i ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/merged_abundance_table_species.txt \
   -o ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/output_images/abundance_heatmap_species.png \
   --ftop 25 \
   --f_dist_f braycurtis \
   --s_dist_f braycurtis \
   --cell_aspect_ratio 0.5 \
   -l \
   --flabel_size 6 \
   --slabel_size 6 \
   --max_flabel_len 100 \
   --max_slabel_len 100 \
   --minv 0.1 \
   --dpi 300

   # creating files for graphlan
   echo "Running graphlan"
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/graphlan_output
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/graphlan_output/output_images
   export2graphlan.py \
   --skip_rows 1,2 \
   -i ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/merged_abundance_table.txt \
   --tree ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/graphlan_output/merged_abundance.tree.txt \
   --annotation ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/graphlan_output/merged_abundance.annot.txt \
   --most_abundant 100 \
   --abundance_threshold 1 \
   --least_biomarkers 10 \
   --annotations 5,6 \
   --external_annotations 7 \
   --min_clade_size 1

   graphlan_annotate.py \
   --annot ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/graphlan_output/merged_abundance.annot.txt \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/graphlan_output/merged_abundance.tree.txt \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/graphlan_output/merged_abundance.xml

   ${TOOLS_FOLDER}/graphlan/graphlan.py \
   --dpi 300 \
   --size 15 \
   --format pdf \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/graphlan_output/merged_abundance.xml \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan2/graphlan_output/output_images/merged_abundance.pdf
   #  --external_legends

   echo "DONE running metaphlan!"
}

run_diamond(){
   echo "Running diamond"

   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/diamond_output

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.12.interleaved.final.clean.fa | awk '{print $NF}')
   for i in ${finallist}
   do
	#if [ -f "${ANALYSIS_FOLDER}/reference_classification/diamond_output/$(basename ${i%12*}.rma)" ]; then
	#echo "diamond and MEGAN:daa2rma already exists"
	#else
	#${TOOLS_FOLDER}/diamond/diamond blastx \

	diamond blastx \
	--threads 70 \
	--query ${i} \
	--db ${REFERENCE_FOLDER}/reference_database/nr.dmnd \
	--daa ${ANALYSIS_FOLDER}/reference_classification/diamond_output/$(basename ${i%12*})

	echo "diamond done...."

	/bioinfo/MEGAN/megan-6/tools/daa2rma \
	--in ${ANALYSIS_FOLDER}/reference_classification/diamond_output/$(basename ${i%12*}.daa) \
	--acc2taxa ${REFERENCE_FOLDER}/reference_database/megan_ref/nucl_acc2tax-Jul2019.abin \
	--acc2interpro2go ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2interpro-Jul2019X.abin \
	--acc2seed  ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2seed-May2015XX.abin  \
	--acc2eggnog ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2eggnog-Jul2019X.abin \
	-fwa true \
	-t 70 \
	--out ${ANALYSIS_FOLDER}/reference_classification/diamond_output/$(basename ${i%.12*}.rma)

	echo "MEGAN:daa2rma done..."
     #fi

   done

   echo "DONE running diamond!"
}


run_humann(){
   echo "Running humann"
   #
   #conda activate mpa
   #/home/gala0002/anaconda3/envs/mpa/lib/python3.7/site-packages/metaphlan/
   #${TOOLS_FOLDER}/metaphlan2/
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/humann
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.interleaved.final.clean.fa | awk '{print $NF}')
   for i in ${finallist}
   do
      humann \
      --input ${i} \
      --metaphlan /home/gala0002/anaconda3/envs/mpa/lib/python3.7/site-packages/metaphlan/ \
      --output ${ANALYSIS_FOLDER}/reference_classification/humann \
      --nucleotide-database ${REFERENCE_FOLDER}/reference_database/humann_databases/chocophlan \
      --protein-database ${REFERENCE_FOLDER}/reference_database/humann_databases/uniref \
      --threads 70
   done

   genefamilieslist=$(ls -d ${ANALYSIS_FOLDER}/reference_classification/humann/*_genefamilies.tsv | awk '{print $NF}')
   for i in ${genefamilieslist}
   do
      humann_renorm_table \
      --input $i \
      --output ${i%_genefamilies*}_genefamilies_relab.tsv \
      --units relab

      humann_join_tables \
      --input ${ANALYSIS_FOLDER}/reference_classification/humann \
      --output ${i%_genefamilies*}_genefamilies.tsv \
      --file_name genefamilies_relab

      humann_join_tables \
      --input ${ANALYSIS_FOLDER}/reference_classification/humann \
      --output ${i%_genefamilies*}_pathcoverage.tsv \
      --file_name pathcoverage

      humann_join_tables \
      --input ${ANALYSIS_FOLDER}/reference_classification/humann \
      --output ${i%genefamilies*}_pathabundance.tsv \
      --file_name pathabundance_relab
   done

   echo "DONE running humann!"
}


#install krken2 and kraken2 database
install_kraken(){
   echo "Checking and installing kraken"
   if [ -d "${TOOLS_FOLDER}/kraken2-2.0.8-beta" ]; then
      echo "kraken already installed"
   else
      #Installing kraken 
      echo "Installing kraken"
      cd ${TOOLS_FOLDER}
      wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz
      tar xzf v2.0.8-beta.tar.gz
      cd kraken2-2.0.8-beta
      sh install_kraken2.sh ${TOOLS_FOLDER}/kraken2-2.0.8-beta
      cd ..
      rm v2.0.8-beta.tar.gz
      #export PATH=`pwd`/kraken2/:$PATH
      #export PATH=`pwd`/kraken2-build/:$PATH
      #export PATH=`pwd`/kraken2-inspect/:$PATH
   fi
   echo "DONE checking and installing kraken!"
}

download_krakendb(){
   echo "Checking and downloading kraken database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/kraken2DB" ]; then
      echo "kraken database already installed"
   else
      mkdir -p ${REFERENCE_FOLDER}/reference_database/
      ln -s ${LINKPATH_DB}/reference_database/kraken2DB ${REFERENCE_FOLDER}/reference_database/kraken2DB
   fi
   echo "DONE checking and downloading kraken database!"
}

#install metaphlan and dependencies
install_metaphlan(){
   echo "Checking and installing metaphlan"
   if [ -d "${TOOLS_FOLDER}/metaphlan2" ]; then
      echo "metaphlan already installed"
   else
      #Installing metaphlan 
      echo "Installing metaphlan"
      cd ${TOOLS_FOLDER}
      git clone https://github.com/biobakery/MetaPhlAn2 metaphlan2
      git clone https://github.com/biobakery/graphlan
      #install hclust
      #wget https://bitbucket.org/nsegata/hclust2/get/tip.zip
      #unzip tip.zip
      #rm tip.zip
      git clone https://github.com/SegataLab/hclust2
      #export PATH=`pwd`/graphlan/:$PATH
   fi
   echo "DONE checking and installing metaphlan!"
}

download_metaphlan(){
   echo "Checking and downloading metaphlan database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/metaphlan2" ]; then
      echo "metaphlan database already installed"
   else
      mkdir -p ${REFERENCE_FOLDER}/reference_database/
      ln -s ${LINKPATH_DB}/reference_database/metaphlan2 ${REFERENCE_FOLDER}/reference_database/metaphlan2
   fi
   echo "DONE checking and downloading metaphlan database!"
}

#install megan and download dependency files
install_megan(){
   echo "Checking and installing megan"
   if [ -d "${TOOLS_FOLDER}/megan" ]; then
      echo "megan already installed"
   else
      #Installing megan 
      echo "Installing megan"
      cd ${TOOLS_FOLDER}
      mkdir -p ${TOOLS_FOLDER}/megan 
      cd megan
      wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_6_16_4.sh
      sh MEGAN_Community_unix_6_16_4.sh
      cd ..
   fi
   echo "DONE checking and installing megan!"
}

download_megandb(){
   echo "Checking and downloading megan database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/megan_ref" ]; then
      echo "megan database already installed"
   else
      mkdir -p ${REFERENCE_FOLDER}/reference_database/
      ln -s ${LINKPATH_DB}/reference_database/megan_ref ${REFERENCE_FOLDER}/reference_database/megan_ref
   fi
   echo "DONE checking and downloading megan database!"
}

#install diamond
install_diamond(){
   echo "Checking and installing diamond"
   if [ -d "${TOOLS_FOLDER}/diamond" ]; then
      echo "diamond already installed"
   else
      #Installing diamond 
      echo "Installing diamond"
      cd ${TOOLS_FOLDER}
      mkdir -p  ${TOOLS_FOLDER}/diamond/
      cd diamond
      wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz
      tar xzf diamond-linux64.tar.gz
   fi
   echo "DONE checking and installing diamond!"
}

download_nrdb(){
   echo "Checking and downloading NR database"
   if [ -f "${REFERENCE_FOLDER}/reference_database/nr.dmnd" ]; then
      echo "NR database already installed"
   else
      mkdir -p ${REFERENCE_FOLDER}/reference_database/
      ln -s ${LINKPATH_DB}/reference_database/nr.dmnd ${REFERENCE_FOLDER}/reference_database/nr.dmnd
   fi
   echo "DONE checking and downloading NR database!"
}

# install humann and dependencies
install_humann(){
   echo "Checking and installing humann"
   if [ -d "${TOOLS_FOLDER}/humann" ]; then
      echo "humann already installed"
   else
      cd ${TOOLS_FOLDER}
      git clone https://github.com/biobakery/humann humann
      cd humann
      python3 setup.py install --user
   fi
   echo "DONE checking and installing humann!"
}

download_humanndb(){
   echo "Checking and downloading humann database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/humann_databases" ]; then
      echo "humann database already installed"
   else
      mkdir -p ${REFERENCE_FOLDER}/reference_database
      ln -s ${LINKPATH_DB}/reference_database/humann_databases ${REFERENCE_FOLDER}/reference_database/humann_databases
   fi
   echo "DONE checking and downloading humann database!"
}

install_graphlan(){
   echo "Checking and installing graphlan"
   if [ -d "${TOOLS_FOLDER}/graphlan" ]; then
      echo "graphlan already installed"
   else
      echo "Installing graphlan"
      git clone https://github.com/biobakery/graphlan
      #cd graphlan

   fi
   echo "DONE checking and installing graphlan!"

}

ref_analysis_main
