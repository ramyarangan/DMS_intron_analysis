import sys

dat_filename = sys.argv[1]
fasta_filename = sys.argv[2]

f = open(dat_filename)
dat_lines = f.readlines()
f.close()

f = open(fasta_filename, 'w')

for ii in range(int(len(dat_lines)/2)):
	dat_line = dat_lines[ii * 2]
	dat_items = dat_line.split(" ")[1].split("\t")
	chrnum = dat_items[0]
	start_pos = dat_items[1]
	end_pos = dat_items[2]
	strand_dir = dat_items[5].replace('\n', '')
	tag = chrnum + ":" + start_pos + "-" + end_pos + "(" + strand_dir + ")"
	f.write(">%s\n" % tag)

	seq = dat_lines[ii * 2 + 1].replace('\n', '')
	f.write("%s\n" % seq)

f.close()
