#! /usr/bin/env python

list_pairs = open("Mvs11_MvslD_70sim50cov_singlecopy_ortholog_list.txt",'r')
MvSl = open("MvSl-1064-A1-R4_gene1.fasta",'r')
MvDv = open("MldSilD_cds1.fasta", 'r')

MvSl_dict = {}
MvDv_dict = {}
keep = "no"

temp_file_name = "outfile" + "tempy.temp"
temp_file = open(temp_file_name, "w")
for line1 in MvSl:

    if line1.startswith(">"):
        line1 = line1.strip("\n")
        fullline1 = line1
        line1 = line1.split(">")
        rt_id = line1[1]
        temp_file.write(">" + rt_id + "\t")
    else:
        seq = line1.strip("\n")
        temp_file.write(seq + "\n")

temp_file.close()

MvSl_seq_set = set()
temp_file2 = open(temp_file_name)
for line in temp_file2:
    if line.startswith(">"):
        line = line.rstrip("\n").split("\t")
        seqid = line[0]
        MvSl_seq_set.add(seqid)
        seq = line[1]
        MvSl_dict[seqid] = seq

temp_file2.close()

### file2
temp_file_name = "outfile" + "tempy2.temp"
temp_file = open(temp_file_name, "w")
for line1 in MvDv:
    if line1.startswith(">"):
        line1 = line1.strip("\n")
        fullline1 = line1
        line1 = line1.split(">")
        rt_id = line1[1]
        temp_file.write(">" + rt_id + "\t")
    else:
        seq = line1.strip("\n")
        #print(seq)
        temp_file.write(seq + "\n")

temp_file.close()

MvDv_seq_set = set()
temp_file2 = open(temp_file_name)
for line in temp_file2:
    if line.startswith(">"):
        line = line.rstrip("\n").split("\t")
        seqid = line[0]
        MvDv_seq_set.add(seqid)
        seq = line[1]
        MvDv_dict[seqid] = seq

####
file_N = 0
for line in list_pairs:
    file_N = file_N + 1
    out_seqfile_name = "Mvsl_MvDv_A1pair_" + str(file_N) +  ".fa"
    out_seqfile = open(out_seqfile_name, "w")
    line = line.rstrip("\n").split("\t")
    id_1 = ">" + line[0]
    id_2 = ">" + line[1]
    if id_2 in MvDv_seq_set:
        seq_1 = MvSl_dict.get(id_1)
        seq_2 = MvDv_dict.get(id_2)
        out_seqfile.write(str(id_1) + "\n" + str(seq_1) + "\n" + str(id_2) + "\n" + str(seq_2) + "\n")
