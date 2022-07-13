import sys

fasta_in = sys.argv[1]
fasta_out = sys.argv[2]

f = open(fasta_in)
fasta_lines = f.readlines()
f.close()

f = open(fasta_out, 'w')

for ii in range(int(len(fasta_lines)/2)):
	f.write("%s" % fasta_lines[ii * 2])
	seq = fasta_lines[ii * 2 + 1].replace('\n', '')
	f.write("%s\n" % seq[50:(-50)])

f.close()

