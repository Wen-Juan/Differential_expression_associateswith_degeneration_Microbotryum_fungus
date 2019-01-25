## extract_ortho_seqs.py

import sys
import os
import getopt
import decimal
import numpy as np
import re


try:
	opts, args = getopt.getopt(sys.argv[1:], 'f:m:o:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

input_fasta = None
orth_matrix_filename  = None
output_dir_name = "extract_ortho_seqs_outdir"


#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** fasta_file_tidier.py | Written by DJP, 15/01/18 in Python 3.5 in Lausanne ****\n")
		print("Takes a fasta file, and ortholog matrix file")
		
		print("\n***** USAGE *****\n")		
		print("\npython3 extract_ortho_seqs.py -f [fasta file name] -m [matrix_file] -o [output filename]\n\n")
		
		sys.exit(2)
		
	elif opt in ('-f'):
		input_fasta = arg
	elif opt in ('-m'):
		orth_matrix_filename = arg
	elif opt in ('-o'):
		output_dir_name = arg
	else:
		print("i dont know")
		sys.exit(2)


##### FIRST unwrap fasta - precautionary will be necessary for some files 
### note making a temp unwrapped fasta file  - removed at end

output_fasta_name = input_fasta + ".TEMP_extract_fasta_file" 

output_file = open(output_fasta_name, "w")
print("\nUnwrapping fasta file")
count = 0
in_file = open(input_fasta)
for line in in_file:
	count = count + 1
	line = line.rstrip("\n")
	if line.startswith(">") and count == 1:
		output_file.write(line + "\n")
	elif line.startswith(">") and count > 1:
		output_file.write("\n" + line + "\n")
	else: 
		output_file.write(line)	

output_file.close()


### add seqs to dictionary
name_list = []
seq_list = []
len_list = []
seq_dict = {}

done = 0
seq_file_1 = open(output_fasta_name)
for line in seq_file_1:
	lineA = line.rstrip("\n")
	if lineA.startswith(">"):
		lineB = lineA.replace(">", "")
		name_list.append(lineB)
	else:
		seq_list.append(lineA)
		done = done + 1
		seq_len = len(lineA)
		len_list.append(seq_len)
			

for element in range(0,len(name_list)):
	name1 = name_list[element]
	seq1 = seq_list[element].replace(" ", "") ## remove gaps if seq comes from gblocks 
	seq_dict[name1] = seq1

#print(seq_dict)
## tidyup
seq_file_1.close()
os.remove(output_fasta_name)

print("Read " + str(done) + " sequences from " + input_fasta)


### make output dir


try:
	os.mkdir(output_dir_name)
except OSError as exc:
	print("\n********* WARNING ************\n" + output_dir_name + " already exists. files may be overwritten.")


#####################################################################################
#### read in otho matrix and output



line_N = 0
species_names = []
orth_matrix_file = open(orth_matrix_filename)
for line in orth_matrix_file:
	line_N = line_N + 1
	line = line.rstrip("\n")
	if line_N == 1:
		species_names = line.replace("_", "").replace("#", "").split("\t")
		print("\nGroup names:")
		print(species_names)
	else:
		seq_want = line.split("\t")
		species_names_out = ""
		for i in species_names:
			species_names_out = species_names_out + "_" + i
		
		OG_name = "OG-" + str(line_N - 1)
		s_N = 0
		curr_out_file = open(os.path.join(output_dir_name, OG_name + species_names_out + ".fa"), "w")
		for s in seq_want:
			
			
			curr_out_file.write(">" + OG_name + "_" + species_names[s_N] + "_" + s + "\n" + seq_dict.get(s) + "\n")
			
			s_N = s_N + 1
			
			#prseq_dict.get(s)
		curr_out_file.close()
		

print("\nOutputted " + str(line_N - 1) + " files to " + output_dir_name + "\n")

print("Done, Dr Ma\n\n\n")
























