#!/usr/bin/python


"""
#Script to Pars and merge annotations 
# Author: GV Saripella
# Organization: Swedish University of Agricultural Sciences (SLU)
# 

#####------Inputs-------
# USAGE: 
# run-Pars-annotation.py INPUT=Inputfile DIR=Directory-fullpath

python run-Pars-annotation.py SWL.tsv /prokka_output/

#############################################------

####------

"""
import sys
import re
import subprocess
from subprocess import *
from subprocess import call

class SearchDB():
	def Search_FLS(self,readfl1,workdir):
		"""
		Calling Search local DB
		"""
		def srchdb_dbCAN(GName,FPATH):
			try:
				True
				cmdFls1 = "LANG=C grep -w '"+str(GName)+"\t' "+str(FPATH)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AnnoPars = cmdFls2.strip().decode().split("\t")
				grepout = "\t".join(AnnoPars[1:])
			except:
				False
				grepout = str("."+"\t"+"."+"\t"+".")
			return grepout
		def srchdb_PFAM(GName,FPATH):
			try:
				True
				cmdFls1 = "LANG=C grep '"+str(GName)+"\t' "+str(FPATH)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AnnoPars = cmdFls2.strip().decode().split("\t")
				grepout = "\t".join(AnnoPars[1:])
			except:
				False
				grepout = str("."+"\t"+"."+"\t"+".")
			return grepout
		if len(readfl1.split(".")[:-1]) ==2:
			readfl = ".".join(readfl1.split(".")[:-1])
		else:
			readfl = readfl1.split(".")[:-1]
		with open(workdir+readfl1,'r') as f1, open(workdir+readfl+"_PFAM-CAZyDB_anno.tsv",'w') as output:
			first_line0 = f1.readline().strip().split("\t")
			first_lines = f1.readlines()
			output.write(str(str("\t".join(first_line0))+"\t"+"PFAM_clust_acc-ID"+"\t"+"PFAM_clust_t-name"+"\t"+"PFAM_clust_acc-desc"+"\t"+"CAZyDB_clust_acc-ID"+"\t"+"CAZyDB_clust_t-name"+"\t"+"CAZyDB_clust_acc-desc"+"\n"))

			for lns in first_lines:
				lns_sp =  lns.strip().split("\t")
				if lns_sp[1] == "CDS":
					searchfl_dbCAN = workdir+readfl+".cdhit.95.dbCAN_clust_join.out"
					searchfl_PFAM = workdir+readfl+".cdhit.95_PFAM_clust_join.out"
					Anno_out_pars1 = srchdb_dbCAN(lns_sp[0],searchfl_dbCAN)
					Anno_out_pars2 = srchdb_PFAM(lns_sp[0],searchfl_PFAM)
					out_print = str("\t".join(lns_sp)+"\t"+Anno_out_pars2+"\t"+Anno_out_pars1)
					output.write(out_print+"\n")

		print("Done seach for ...")
		return None

clF1 = SearchDB().Search_FLS(sys.argv[1],sys.argv[2])


