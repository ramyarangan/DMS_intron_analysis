from os import path 
import numpy as np
from scipy.stats import pearsonr
from matplotlib import pyplot as plt 

## For the genes where we have sufficient coverage from the -pladB control to 
## predict reactivity values, compare reactivity values -pladB vs +pladB
regions = ["chrVII:555830-556307", "chrVII:311015-311526", "chrIV:579478-580017"]
seqs = {
	"chrVII:555830-556307": "GUAUGUUUGGAGGAUACGAAUAACGAUAGAAAACAUGAGUGAAUUUCCGUCCACGAAAAAAUGUUAACAUAAAAUGCAAGAGAACAAUUAAUCGAAUAAUGUUAAAUUAUUGUAAAACAAUGUGUAUGAUGAGGAGGAAUGUACCUAAGCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAACAGCUUUUGCAUAUUCAAUCCAGGCAUAGGGCGACUAUUUAGCACUCAACGAUUUUUAAGCUUGUGUAUUGCUGACAUAAAUUCCGGCUUUAGAAUCCAAUAUUGAAAAACGUGAGUACGCAGAGGAGAUAGAAGAAAAGUAGGAAGUUACCGUUUAUAUUGAUUUGUGAAAUGCAUACUCCGUUGGAUGUGGGGCAACAUAGAUUUAAGUGUGGAUGAAAAUUAUGUGCUCAUUGUGAAAAAAAAGUUUUGCUUUUACUAACAAAUUUUUUUAUUAUUUGUUUUCAAUAG", 
	"chrVII:311015-311526": "GUAUGUAGUUCCAUUUGGAAGAGGGAAUGAAAGAACCAAGACGGUGACUUUUUUUUUAGUGUUGUGCAACCAAUAUGUCGUGUGUAUAUCAUGGUACAGGAGAAUGUCAAUCAGCUAAGUGUACUCAACAUAUUUCUUUGUGUUUUGAUUGCGAACUUUGUAUUACCAUCUCACUGUUGAGACGGCUUAUUUGAGGUAAUAGCUCGAGUAAAUGUACUCUUCCAUCGCAAACUGAGCAAAAAGAAAGUGUGCAUAGCCUUUGUCAUACUUCUCCUUUAUUAUACCAUGAUAUUCAGAACAGUCAUACUGUCUACUCAUUUUACGGCUAUAAAAGGUAACUUUCAUUUAGAUUAUGGAAAGCACUAAUUAUCGCUGUAUCAAAUGGUUGUAGAGAGCGCAAUUAUGAAAAAGAGUUACCACGUUUCUUUUGUUUCGAUAAAAUGUCCAGUUGAAAACCUGUUUUACUAACGAUUUAAAAAUUGUAUUUCAUUACAAUAUUUUUUUUGUACAG",
	"chrIV:579478-580017": "GUAUGUUUAUUAACACCAUAGCGAGAUAUUAAUGCAAAAGUUGCAUUGAAUAGUUCGCUAAAUCAGAUGACACUCUAAUGUGGAAUUCAAAAGUGGAUUUCUAAUAUAAUUUGUCUCUGUCGGAUCACAAUUCUAUUACAAGUUCCGGUGUGUACACAGGUAUAGUUUAUACUGGAGAGUAGUUUCUACUCGCUGUACAUUAGCUGGGUGAUUCCAAUUUCUUUUACAAAUAUGUUGCAUUAGUUUAACAGGUUAUACUAUCUGCCGUUUCUCAGUAUAAUUUACGCCGGAAAAUUACUGAUGGCUAGCCGCCUUUAUGAAUUAGUUUUCACAAAGCUCAUAACAUAACACGUUAACCUAUCGGAGGAGAACCAAGAUUGAAGAAUCACCCGGAAUAGUUAUACUUUAAUGGAAUUGUAUGGUCUGAACGAGGAAAUAUGUCAUGAUACACUUUUCUUCAAGCCAUAUGAAUCUUCAUGUUACUAACAUUCGAUAAAUUUUUUGGAAUAUCCAAUUCCACUAAAUAUUACUUUAAACAG"
}
nd_reac_dir = 'combined_1221/reactivity/reactivity_intron_nd/'
reac_dir = 'combined_1221/reactivity/reactivity_intron/'

# Returns raw reactivities, smoothed linearly interpolated reactivities and coverage (in RPKM)
def get_reactivities(name, reac_dir):
	reac_filename = reac_dir + name + ".txt"
	
	if not path.exists(reac_filename):
		print(reac_filename)
		return np.array([]), 0

	f = open(reac_filename)
	lines = f.readlines()
	f.close()

	cov = float(lines[0].split(' ')[1])
	reacs = []
	cnt = 0
	for line in lines[1:]:
		if line.replace('\n', '') == 'nan':
			reacs += [np.nan]
		else: 
			reacs += [float(line)]
		cnt += 1
	
	reacs = np.array([-999 if np.isnan(x) else x for x in reacs])

	return reacs, cov	

# Pearson correlation between two reactivity lists (no NaN's)
def get_r_val(reacs_1, reacs_2):
	plt.scatter(reacs_1, reacs_2)
	plt.show()
	pearson = pearsonr(np.array(reacs_1), np.array(reacs_2))
	return pearson

def make_bargraphs(nd_reacs, d_reacs, region):
	nd_reacs[nd_reacs < 0] = 0
	d_reacs[d_reacs < 0] = 0
	seq = seqs[region]
	plt.figure(figsize=(10, 0.5))
	plt.bar(list(range(len(seq))), nd_reacs, color='red', width=1)
	plt.savefig('../figures/' + region + '_nd' + '.png', format='png', dpi=300)

	plt.figure(figsize=(10, 0.5))
	plt.bar(list(range(len(seq))), d_reacs, color='red', width=1)
	plt.savefig('../figures/' + region + '_d' + '.png', format='png', dpi=300)

for region in regions: 
	nd_reac, _ = get_reactivities(region, nd_reac_dir)
	d_reac, _ = get_reactivities(region, reac_dir)
	nd_reac = np.array(nd_reac)
	d_reac = np.array(d_reac)

	make_bargraphs(nd_reac, d_reac, region)
	mask = np.logical_and(nd_reac > 0, d_reac > 0)
	nd_reac = nd_reac[mask]
	d_reac = d_reac[mask]

	r_val = get_r_val(nd_reac, d_reac)
	print(region)
	print(r_val)
