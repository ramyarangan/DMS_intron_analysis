from pymol import cmd, util

cmd.load('4v88_bundle1_bundle2.pdb')

cmd.set('dot_solvent', 1)
cmd.set('dot_density', 1)
cmd.set('solvent_radius', 3)
stored.residues = []
cmd.iterate('name P and chain y', 'stored.residues.append(resi)')

sasa_per_residue = []
for i in stored.residues[1000:]:
	cur_area = cmd.get_area('chain y and resi %s and ((resn A and name N1) or (resn C and name N3))' % i) #  and ((resn A and name N1) or (resn C and name N3))
	print("%s %f" % (i, cur_area))
	sasa_per_residue.append(cur_area)

print(sum(sasa_per_residue))
print(cmd.get_area('all'))  # just to check that the sum of sasa per residue equals the total area
