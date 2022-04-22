#!/usr/bin/evn bash

### Annotation with KEGG
#KEGG annotation procedure


#KOBAS
#Download KOBAS files and unzip

cd tools/
mkdir -p ./KOBAS/sqlite3
mkdir -p ./KOBAS/seq_pep
cd ./KOBAS/sqlite3
wget http://kobas.cbi.pku.edu.cn/download/sqlite3/ko.db.gz
gunzip ko.db.gz
cd /data/KOBAS/seq_pep
wget http://kobas.cbi.pku.edu.cn/download/seq_pep/ko.pep.fasta.gz
gunzip ko.pep.fasta.gz


#Make a blastable database

makeblastdb -in ko.pep.fasta -dbtype prot

#Make a Diamond searchable database

diamond makedb --in ko.pep.fasta -d ko.pep

#first, run a diamond search agains the KOBAS KO database

diamond blastp -d /data/DATABASES/KOBAS/seq_pep/ko -q sequence_files/SDB_ONE.faa -o SDB_ONE.faa.dmd -e 1e-10 -k 1
gunzip SDB_ONE.faa.dmd.gz

#Extract KO numbers and get gene IDs from faa file

cut -f 1,2 SDB_ONE.faa.dmd > SDB_ONE.ORF_G_IDs
grep '>' sequence_files/SDB_ONE.faa | sed 's/>//g' > SDB_ONE.ORF_IDs
cut -f 2 SDB_ONE.ORF_G_IDs > SDB_ONE.G_IDs

#Make an SQlite database and add the tables

sqlite3 annotate.sqlite

.separator \t
create table ORF_G_IDs (ORF, GID);
.import SDB_ONE.ORF_G_IDs ORF_G_IDs

.separator " "
create table KoGenes (KO, GID);
.import KoGenes KoGenes

.separator \t
create table ORFS (ORF);
.import SDB_ONE.ORF_IDs ORFS

CREATE TABLE IF NOT EXISTS out as
SELECT ORF_G_IDs.ORF, ORF_G_IDs.GID, KoGenes.KO FROM ORF_G_IDs JOIN KoGenes ON ORF_G_IDs.GID = KoGenes.GID;

CREATE TABLE IF NOT EXISTS allout as
SELECT * FROM ORFS LEFT JOIN out ON ORFS.ORF = out.ORF

.separator "\t"
.output SDB_ONE.G_ID_KO_ORF_ID
SELECT * FROM allout


