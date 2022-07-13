f = open('ares_standard_introns.csv')
lines = f.readlines()
f.close()

f = open('standard_introns.fa', 'w')

for line in lines[1:]:
	line_items = line.split(',')
	id_string = line_items[0] + ":" + line_items[1] + "-" \
		+ line_items[2] + " " + line_items[6]
	f.write(">%s\n" % id_string)
	f.write("%s" % line_items[7])

f.close()