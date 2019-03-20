file2 = open("list.txt","r")

for line in file2:
	line = line.split(" ")
	print ('\n'.join(line))
