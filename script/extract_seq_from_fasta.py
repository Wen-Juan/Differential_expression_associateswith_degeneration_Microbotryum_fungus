infile1 = open("Auto_DE_a2.txt","r")
infile2 = open("MvSl-1064-A2-R4_cds.fasta","r")
outfile = open("Auto_DE_a2_cds.fasta","w")

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
