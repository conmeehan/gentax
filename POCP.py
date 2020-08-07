#! /usr/bin/env python
import sys
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

"""
Author: Conor Meehan 23/01/2019

Script to calculate the percentage of conserved proteins (POCP) between two isolates
This is calculated as follows (based on Qin et al, Journal of Bacteriology 2014):
Proteins from the query genome were considered conserved when they had a 
BLAST match with an E value of less than 1e-5, 
a sequence identity of more than 40%, 
and an alignable region of the query protein sequence of more than 50%.
The percentage of conserved proteins (POCP) between two genomes was calculated as 
[(C1+C2)/(T1+T2)]*100%, 
where C1 and C2 represent the conserved number of proteins in the two genomes being compared, respectively, 
and T1 and T2 represent the total number of proteins in the two genomes being compared, respectively.


Input:
A BLASTp output in tabular format (outfmt 6) for isolate 1 as query and isolate 2 as subject
A BLASTp output in tabular format (outfmt 6) for isolate 2 as query and isolate 1 as subject
A protein file (e.g. faa file) for the isolate 1
A protein file (e.g. faa file) for the isolate 2 

Output:
Save file with the PCOP inside (no rounding)
File is named as isolate1_isolate2_PCOP.txt
Names are taken from the protein faa files with the folder names (everything before /) and the faa stripped


NOTE: requires biopython
NOTE: Names in faa files are cut at the first space, as they would be in the BLAST output

Usage:
python POCP.py --blast1 genome1QueryBLASTpFile --blast2 genome2QueryBLASTpFile --protein1 faaFileSpecies1 --protein2 faaFileSpecies2 

"""

#parse the inputs
parser = argparse.ArgumentParser(description= 'Script to calculate the percentage of conserved proteins (POCP) between two genomes. See script header for details.')
parser.add_argument('--blast1', required=True, help='A BLASTp output in tabular format (outfmt 6) for isolate 1 as query and isolate 2 as subject')
parser.add_argument('--blast2', required=True, help='A BLASTp output in tabular format (outfmt 6) for isolate 2 as query and isolate 1 as subject')
parser.add_argument('--protein1', required=True, help='A protein file (e.g. faa file) for the isolate 1')
parser.add_argument('--protein2', required=True, help='A protein file (e.g. faa file) for the isolate 2')

args = parser.parse_args()


#create the save file
save1=re.sub(".faa","", args.protein1)
save1=re.sub(".*/","",save1)
save2=re.sub(".faa","", args.protein2)
save2=re.sub(".*/","",save2)

try:
	save=open(save1+"_"+save2+"_PCOP.txt",'w')
except IOError:
	print('no room for save file')
	sys.exit()


#open the BLAST files
try:
	blast1F=open(args.blast1, 'r')
except IOError:
	print("\n BLAST 1 not found in folder.")
	sys.exit()
try:
	blast2F=open(args.blast2, 'r')
except IOError:
	print("\n BLAST 2 not found in folder.")
	sys.exit()

#read in the protein file of isolate 1 and for each protein, get the length and name (SeqIO auto cuts it at first space) and save to a dictionary
prots1all = SeqIO.to_dict(SeqIO.parse(args.protein1, "fasta"))
prots1={}
for prot in prots1all:
	prots1[prot]=len(prots1all[prot].seq)
#repeat for the protein file of isolate 2
prots2all = SeqIO.to_dict(SeqIO.parse(args.protein2, "fasta"))
prots2={}
for prot in prots2all:
	prots2[prot]=len(prots2all[prot].seq)

#get the totals
T1=len(prots1)
T2=len(prots2)

#go through the first BLAST file and for each protein that has a match meeting the minimum requirements, add 1 to the counter and the protein name to the list (to avoid counting the same protein twice)
C1=0.0
matches1=[]
while 1:
	line=blast1F.readline()
	if not line:
		break
	line=line.rstrip()
	sections=line.split()
	name=sections[0]
	id=float(sections[2])
	lenMatch=float(sections[3])
	evalue=float(sections[10])
	#if already matched, skip
	if name in matches1:
		continue
	#check evalue is less than 1e-05
	if evalue>1e-05:
		continue	
	#check id is greater than 40%
	if id<40:
		continue
	#get the length of the match as a percetnage of the total query protein length
	queryProtLen=prots1[name]
	matchPC=lenMatch/queryProtLen
	if matchPC<0.5:
		continue
	
	#if reached this stage, then have passed all checks. Increase the counter by 1 and add the name to the matched list
	matches1.append(name)
	C1+=1
blast1F.close()

#repeat for isolate 2 as query
C2=0.0
matches2=[]
while 1:
	line=blast2F.readline()
	if not line:
		break
	line=line.rstrip()
	sections=line.split()
	name=sections[0]
	id=float(sections[2])
	lenMatch=float(sections[3])
	evalue=float(sections[10])
	#if already matched, skip
	if name in matches2:
		continue
	#check evalue is less than 1e-05
	if evalue>1e-05:
		continue	
	#check id is greater than 40%
	if id<40:
		continue
	#get the length of the match as a percetnage of the total query protein length
	queryProtLen=prots2[name]
	matchPC=lenMatch/queryProtLen
	if matchPC<0.5:
		continue
	
	#if reached this stage, then have passed all checks. Increase the counter by 1 and add the name to the matched list
	matches2.append(name)
	C2+=1
blast2F.close()

#calculate the PCOP
PCOP=((C1+C2)/(T1+T2))*100


#save to file
save.write(str(PCOP)+"\n")
save.close()
sys.exit()





