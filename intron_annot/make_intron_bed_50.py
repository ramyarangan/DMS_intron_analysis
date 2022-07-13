# Make intron bed file for intersecting with reads aligned to 
# window +-50 around intron
import sys

fasta_filename = sys.argv[1] # e.g. standard_allsize_extend50.fa
bed_filename = sys.argv[2]

f = open(fasta_filename)
fasta_lines = f.readlines()
f.close()

f = open(bed_filename, 'w')

for ii in range(int(len(fasta_lines)/2)):
	tag_line = fasta_lines[ii * 2][1:].replace('\n', '')
	seq = fasta_lines[ii * 2 + 1].replace('\n', '')
	f.write("%s\t%d\t%d\n" % (tag_line, 50, len(seq) - 50))

f.close()
