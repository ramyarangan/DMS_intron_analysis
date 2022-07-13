# Goal: Make two files based on the standard_introns set: 
# bedfile for the full pre-mRNA (intron-containing gene) 
# fasta file for the coding sequence (no intron)

f = open('standard_introns.fa')
intron_lines = f.readlines()
f.close()

introns = {}
for ii in range(int(len(intron_lines)/2)):
	intron_line = intron_lines[ii * 2]
	name = intron_line.split(' ')[1].replace('\n', '')
	intron_line = intron_line.split(' ')[0].replace('>', '')
	if name not in introns:
		introns[name] = [intron_line]
	else:
		introns[name] += [intron_line]

f = open('all_coding.fasta')
coding_annots = f.readlines()
f.close()

f = open('coding_orfs_introns.fa', 'w')
f_2 = open('intron_name_to_coding_orf.txt', 'w')

found = []

do_write = False
for coding_annot in coding_annots:
	if coding_annot[0] == '>':
		do_write = False
		gene_name = coding_annot.split(' ')[1]
		backup_gene_name = coding_annot.split(' ')[0].replace('>','')
		if gene_name in introns.keys():
			for intron in introns[gene_name]:
				f_2.write('%s\t%s\n' % (intron, introns[gene_name][0]))
			f.write('>%s %s\n' % (introns[gene_name][0], gene_name))
			found += [gene_name]
			do_write = True
		elif backup_gene_name in introns.keys():
			for intron in introns[backup_gene_name]:
				f_2.write('%s\t%s\n' % (intron, introns[backup_gene_name][0]))
			f.write('>%s %s\n' % (introns[backup_gene_name][0], gene_name))
			found += [backup_gene_name]
			do_write = True
	elif do_write:
		f.write(coding_annot)

f.close()
f_2.close()

for intron_gene in introns.keys():
	if intron_gene not in found:
		print(intron_gene)
print(len(found))
print(len(introns.keys()))