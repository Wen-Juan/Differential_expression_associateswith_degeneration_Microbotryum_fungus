###This python script can extract gene sequences from a fasta file, using an index input file of gene list.

infile1 = open("Auto_NonDE_a1list.txt","r")
infile2 = open("MvSl-1064-A1-R4_prot1.fasta","r")
outfile = open("Auto_NonDE_a1list_prot.fasta","w")

table = {}
keep = "no"
nothing = 0

for line1 in infile1:
	line1 = line1.strip("\n")
	line1 = line1.split()
	seq_I_want = line1[0]
	table[seq_I_want] = seq_I_want
	print(line1)

for line2 in infile2:
	if line2.startswith(">"):
		line2 = line2.strip("\n")
		fullline2 = line2
		line2 = line2.split(">")
		seq_name = line2[1]
		if seq_name in table:
			outfile.write(fullline2 + "\n")
			keep = "yes"
		else:
			keep = "no"
	else:
		if keep == "yes":
			outfile.write(line2)
		else:
			nothing += 1
