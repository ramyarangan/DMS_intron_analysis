import xml.etree.ElementTree as ET
from os import path, remove, makedirs
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import numpy as np 
from matplotlib import pyplot as plt 

exp_read_counts = {'d3': 557428287, 'd45': 380745898, 'd': 938174185, 'nd': 169742497, 'rouskin': 296022844}
exp_read_lens = {'d3': 300, 'd45': 300, 'd': 300, 'nd': 300, 'rouskin': 50, 'd3_0.25': 300, 'd45_0.25': 300}

intron_fasta_file = "../intron_annot/standard_introns.fa"
coding_fasta_file = "../intron_annot/coding_orfs_introns.fa"
intron_to_coding_name_file = "../intron_annot/intron_name_to_coding_orf.txt"

d_stats_file = 'combined_1221/run_data/d_stats.txt'
d_coding_stats_file = 'combined_1221/run_data/d_coding_stats.txt'
d_premrna_stats_file = 'combined_1221/run_data/d_premrna_stats.txt'
d3_stats_file = 'combined_1221/run_data/d3_stats.txt'
d45_stats_file = 'combined_1221/run_data/d45_stats.txt' #'combined/run_data/d_stats.txt'
d3_coding_stats_file = 'combined_1221/run_data/d3_coding_stats.txt'
d45_coding_stats_file = 'combined_1221/run_data/d45_coding_stats.txt'
nd_stats_file = 'combined_1221/run_data/nd_stats.txt'
nd_coding_stats_file = 'combined_1221/run_data/nd_coding_stats.txt'

# Avoid snoRNAs since most reads for these are for the snoRNA rather than the intron
snorna_introns = ['chrI:142253-142619', 'chrVII:364964-365432', 'chrXIII:499878-500151', 'chrXI:283095-283421',\
	'chrXII:856574-857057', 'chrXIII:163308-163716', 'chrXIV:721770-722302', 'chrXVI:173665-174072', 'chrXVI:173162-174072']
stem_introns = ['chrVII:859260-859473', 'chrIV:1355236-1355550', 'chrII:592416-592768', 'chrXII:522669-523028',\
	'chrXIV:443826-444171', 'chrXIII:225891-226289', 'chrIV:491559-491898', 'chrIII:177906-178213',\
	'chrXVI:75985-76223', 'chrV:166771-166874']
special_rpg_introns = ['chrVII:439093-439323', 'chrXVI:404956-405457', 'chrII:604514-604927', 'chrX:73796-74204']

def get_rpg_introns():
	rpg_introns = []

	f = open(intron_fasta_file)
	lines = f.readlines()
	f.close()

	for line in lines: 
		if len(line) > 0 and line[0] == '>':
			tag = line.split(" ")[0].replace(">", "")
			gene_name = line.split(' ')[1].replace('\n', '')
			if (gene_name[:2] == "RP"):
				rpg_introns += [tag]

	rpg_introns += ["YML6", "MRPL44"]
	return rpg_introns
all_rpg_introns = get_rpg_introns()

# Gets sequences in FASTQ file into a couple different data structures
def get_names_seqs(fasta_file):
	f = open(fasta_file)
	lines = f.readlines()
	f.close()

	names = []
	seqs = []
	name_seq_dict = {}
	cur_seq = ''
	for line in lines:
		if line[0] == '>':
			if cur_seq != '':
				seqs += [cur_seq]
				name_seq_dict[names[-1]] = cur_seq
			cur_seq = ''
			name = line.split(' ')[0].replace('>', '').replace('\n', '')
			names += [name]
		else:
			cur_seq += line.replace('\n', '')
	seqs += [cur_seq]
	name_seq_dict[names[-1]] = cur_seq
	return names, seqs, name_seq_dict

# Get dictionary of name of construct to length of each construct from a fasta file
def fill_length_dict(filename):
	f = open(filename)
	seqs = f.readlines()
	f.close()

	len_dict = {}
	tag = ""
	cur_seq = ""
	for seq in seqs:
		if len(seq) > 0 and seq[0] == ">":
			if tag != "":
				len_dict[tag] = len(cur_seq)
			tag = seq.split(" ")[0].replace(">", "").replace('\n', '')
			cur_seq = ""
		else:
			cur_seq += seq.replace("\n", "")
	if tag != "":
		len_dict[tag] = len(cur_seq)	
	return len_dict

# Get a dictionary from the name of the intron to the name of the coding sequence
def fill_name_dict(filename):
	f = open(filename)
	tag_lines = f.readlines()
	f.close()

	name_dict = {}
	for tag_line in tag_lines:
		intron_name = tag_line.split("\t")[0].replace("\n", "")
		coding_name = tag_line.split("\t")[1].replace("\n", "")
		name_dict[intron_name] = coding_name 
	return name_dict

# RPKM mode: RPKM normalized using the construct lengths from the fasta file and experiment read counts
# Not RPKM mode: normalize to construct lengths from the fasta file, multiply by read length. This 
#               is the expected coverage per position in the construct.
def get_norm_cov(stats, lengths_dict, experiment, do_rpkm, do_rpk, inc_names=[]):
	cov_dict = {}
	for stat in stats:
		if len(stat) > 0 and stat[0] == 'c':
			name = stat.split('\t')[0]
			if inc_names != [] and name not in inc_names:
				continue
			coverage = int(stat.split('\t')[1].replace('\n',''))
			coverage = coverage * exp_read_lens[experiment]/lengths_dict[name] # coverage * 1000/lengths_dict[name]
			# if do_rpkm:
			# 	permil_factor = exp_read_counts[experiment]/1000000
			# 	coverage = coverage/permil_factor
			# if not do_rpkm and not do_rpk:
			# 	coverage = coverage * exp_read_lens[experiment]
			cov_dict[name] = coverage
	return cov_dict

# Fill coverage dictionary with values from the RNAframework stats files
# RNAframework stats files is the number of reads that map onto each sequence
def fill_coverage_dict(intron_filename, experiment, fasta_file, do_rpkm=False, do_rpk=False):
	f = open(intron_filename)
	intron_stats = f.readlines()
	f.close()

	intron_lens = fill_length_dict(fasta_file)

	intron_cov_dict = get_norm_cov(intron_stats, intron_lens, experiment, do_rpkm, do_rpk)

	return intron_cov_dict

def get_max_num_introns():
	d3_cov = fill_coverage_dict(d3_stats_file, 'd3', intron_fasta_file)
	return len(d3_cov.keys()) - len(snorna_introns)

# Retained intron fraction is the RPKM coverage of intron region (unspliced transcripts) over
# the RPKM coverage of the coding region (which is from unspliced + spliced transcripts)
def get_ri_fractions(intron_rpkms, coding_rpkms):
	name_dict = fill_name_dict(intron_to_coding_name_file)

	ri_fractions = {}
	for intron_key in intron_rpkms.keys():
		if intron_key in name_dict.keys():
			if name_dict[intron_key] in coding_rpkms.keys():
				coding_key = name_dict[intron_key]
				ri_fractions[intron_key] = intron_rpkms[intron_key]/coding_rpkms[coding_key]
	return ri_fractions

# Get reactivity values from xml files for an intron
def get_reac(cur_intron, exp, reac_dir='combined_1221/'):
	cur_reac_file = reac_dir + 'rfnorm_reactivity/rfnorm_' + exp + '/' + cur_intron + '.xml'
	if path.exists(cur_reac_file):
		tree = ET.parse(cur_reac_file)
		reac_line = tree.getroot()[0][1].text.replace('\n', '')
		reac_line = reac_line.replace('\t', '')
		return reac_line.split(',')
	else:
		return []

def get_float_reac(cur_intron, exp):
	reac = get_reac(cur_intron, exp)
	return np.array([float(x) if x != 'NaN' else np.nan for x in reac])
	
# Checks if at least 10% of the values are > 0 and not NaN.
def has_data(reactivities, thresh=0.1):
	num_values = 0
	for reac_value in reactivities:
		if reac_value == 'NaN':
			continue
		if float(reac_value) > 0.001:
			num_values += 1
	if num_values == 0:
		return False
	if num_values/len(reactivities) > thresh:
		return True
	return False

# Filter two reactivity lists such that any position that's NaN in either
# list gets filtered out in both. 
def filter_lists(reacs_1, reacs_2):
	reacs_1_filt = []
	reacs_2_filt = []
	if len(reacs_1) != len(reacs_2):
		print("ERROR: Two reactivity lists being compared are not the same length.")
	for ii, reac_1 in enumerate(reacs_1):
		reac_2 = reacs_2[ii]
		if reac_1 != 'NaN' and reac_2 != 'NaN':
			reacs_1_filt += [float(reac_1)]
			reacs_2_filt += [float(reac_2)]
	return [reacs_1_filt, reacs_2_filt]

# Pearson correlation between two reactivity lists (no NaN's)
def get_r_val(reacs_1, reacs_2):
	[reacs_1_filt, reacs_2_filt] = filter_lists(reacs_1, reacs_2)
	#plt.scatter(reacs_1_filt, reacs_2_filt)
	#plt.show()
	pearson = pearsonr(np.array(reacs_1_filt), np.array(reacs_2_filt))
	#print(pearson)
	return pearson

# Spearman rank correlation between two reactivity lists (no NaN's)
def get_p_val(reacs_1, reacs_2):
	[reacs_1_filt, reacs_2_filt] = filter_lists(reacs_1, reacs_2)
	return spearmanr(np.array(reacs_1_filt), np.array(reacs_2_filt))

def write_reactivity(name, exp, cov, outdir):
	reac = get_float_reac(name, exp)

	if not path.exists(outdir):
		makedirs(outdir)

	if len(reac) > 0:
		f = open(outdir + '/' + name + '.txt', 'w')
		f.write('Coverage(RPKM): %f\n' % cov[name])
		print("Coverage for %s in %s: %f" % (name, exp, cov[name]))
		for reac in reac:
			f.write('%f\n' % reac)
		f.close()

# This is for writing final reactivities that are from an alignment combining all experiments
def write_all_reactivities(outdir='reactivity'):
	intron_names, intron_seqs, intron_name_seq_dict = get_names_seqs(intron_fasta_file)

	d_cov = fill_coverage_dict(d_stats_file, 'd', intron_fasta_file, do_rpkm=False)

	d3_cov = fill_coverage_dict(d3_stats_file, 'd3', intron_fasta_file, do_rpkm=False)
	d45_cov = fill_coverage_dict(d45_stats_file, 'd45', intron_fasta_file, do_rpkm=False)

	avg_cov = {}
	for cur_key in d_cov.keys():
		avg_cov[cur_key] = (d3_cov[cur_key] + d45_cov[cur_key])/2

	nd_cov = fill_coverage_dict(nd_stats_file, 'nd', intron_fasta_file, do_rpkm=False)
	nd_coding_cov = fill_coverage_dict(nd_coding_stats_file, 'nd', coding_fasta_file, do_rpkm=False)

	name_dict = fill_name_dict(intron_to_coding_name_file)

	for cur_intron in intron_names:
		if cur_intron not in snorna_introns:
			intron_outdir = outdir + '_intron'
			write_reactivity(cur_intron, 'd', avg_cov, intron_outdir)
			write_reactivity(cur_intron, 'nd', nd_cov, intron_outdir + '_nd')
			if cur_intron in name_dict:
				cur_coding = name_dict[cur_intron]
				coding_outdir = outdir + '_coding'
				write_reactivity(cur_coding, 'nd_coding', nd_coding_cov, coding_outdir + '_nd')
			intron_seq = intron_name_seq_dict[cur_intron]

def plot_nodrug_drug_compare(write_ri=True):
	d3_cov = fill_coverage_dict(d3_stats_file, 'd3', intron_fasta_file, do_rpkm=False)
	d45_cov = fill_coverage_dict(d45_stats_file, 'd45', intron_fasta_file, do_rpkm=False)
	nd_cov = fill_coverage_dict(nd_stats_file, 'nd', intron_fasta_file, do_rpkm=False)
	d3_cov_rpkm = fill_coverage_dict(d3_stats_file, 'd3', intron_fasta_file, do_rpkm=True)

	d3_coding_cov = fill_coverage_dict(d3_coding_stats_file, 'd3', coding_fasta_file, do_rpkm=False)
	d45_coding_cov = fill_coverage_dict(d45_coding_stats_file, 'd45', coding_fasta_file, do_rpkm=False)
	nd_coding_cov = fill_coverage_dict(nd_coding_stats_file, 'nd', coding_fasta_file, do_rpkm=False)

	d3_ri_dict = get_ri_fractions(d3_cov, d3_coding_cov)
	d45_ri_dict = get_ri_fractions(d45_cov, d45_coding_cov)
	nd_ri_dict = get_ri_fractions(nd_cov, nd_coding_cov)
	avg_ri = []
	nd_ri = []
	
	if write_ri:
		f = open("ri_fractions.txt", 'w')

	for intron_key in d3_ri_dict.keys():
		if intron_key not in snorna_introns:
			avg_ri += [(d3_ri_dict[intron_key] + d45_ri_dict[intron_key])/2]
			nd_ri += [nd_ri_dict[intron_key]]
			if write_ri:	
				f.write("%s: %f %f %f %f\n" % (intron_key, avg_ri[-1], nd_ri[-1], 
					(d3_cov[intron_key] + d45_cov[intron_key])/2, nd_cov[intron_key]))
	if write_ri:
		f.close()
	
	plt.scatter(nd_ri, avg_ri, color="black", s=10)
	(xmin, xmax) = plt.xlim()
	plt.plot([xmin, xmax], [xmin, xmax], '--', color="black")
	plt.title('Retained intron fraction for each intron')
	plt.xlabel('-drug intron RPKM / coding RPKM')
	plt.ylabel('+drug intron RPKM / coding RPKM')
	plt.yscale('log', base=10)
	plt.xscale('log', base=10)
	plt.savefig('../figures/nodrug_drug_compare.png', format='png', dpi=300)
	# plt.show()

# Compare the two replicates using data from both sequencing runs
def get_rvals_pvals_covs():
	d3_cov = fill_coverage_dict(d3_stats_file, 'd3', intron_fasta_file, do_rpk=True)
	d3_coding_cov = fill_coverage_dict(d3_coding_stats_file, 'd3', coding_fasta_file, do_rpk=True)

	d45_cov = fill_coverage_dict(d45_stats_file, 'd45', intron_fasta_file, do_rpk=True)
	d45_coding_cov = fill_coverage_dict(d45_coding_stats_file, 'd45', coding_fasta_file, do_rpk=True)

	no_reac_d3 = 0
	no_reac_d45 = 0

	rvals = [] # Pearson correlation (linear)
	pvals = [] # Spearman rank correlation
	covs = []
	covs_coding = []
	names = []

	all_keys = list(d3_cov.keys())
	for cur_intron in all_keys:
		# Skip introns that do not have a coding sequence
		name_dict = fill_name_dict(intron_to_coding_name_file)
		if cur_intron not in name_dict.keys():
			continue
		coding_key = name_dict[cur_intron]
		if coding_key not in d3_coding_cov.keys():
			continue

		if cur_intron not in snorna_introns:
			d3_reac = get_reac(cur_intron, 'd3')
			d45_reac = get_reac(cur_intron, 'd45')

			skip_compare = False
			if not has_data(d3_reac):
				skip_compare = True
				no_reac_d3 += 1
			if not has_data(d45_reac):
				skip_compare = True
				no_reac_d45 += 1

			if not skip_compare:
				r_val = get_r_val(d3_reac, d45_reac)
				p_val = get_p_val(d3_reac, d45_reac)
				avg_cov = (d3_cov[cur_intron] + d45_cov[cur_intron])/2
				avg_cov_coding = (d3_coding_cov[coding_key] + d45_coding_cov[coding_key])/2
				if not np.isnan(r_val[0]) and not np.isnan(p_val.correlation):
					names += [cur_intron]
					rvals += [r_val[0]]
					pvals += [p_val.correlation]
					covs += [avg_cov]
					covs_coding += [avg_cov_coding]

	print('D3 reactivities not passing threshold: ' + str(no_reac_d3) + ' of ' + str(len(d3_cov.keys())-7))
	print('D45 reactivities not passing threshold: ' + str(no_reac_d45) + ' of ' + str(len(d45_cov.keys())-7))
	return rvals, pvals, covs, covs_coding, names

def get_covs_rvals_names_list(log_covs, rvals, names, names_list):
	new_log_covs = []
	new_rvals = []
	for ii, name in enumerate(names): 
		if name in names_list: 
			new_log_covs += [log_covs[ii]]
			new_rvals += [rvals[ii]]
	return (new_log_covs, new_rvals)

def plot_coverage_rval(rvals, covs, coding_covs, names, add_rpgs=True, add_special_rpgs=True):
	num_05 = sum(np.array(rvals) > 0.5)
	# num_05 = sum(np.array(covs) > 1971)
	print("Number of constructs with over 0.5 reactivity correlation: ")
	print(num_05)

	num_05 = sum(np.array(rvals) > 0.6)
	print("Number of constructs with over 0.6 reactivity correlation: ")
	print(num_05)

	num_05 = sum(np.array(rvals) > 0.7)
	print("Number of constructs with over 0.7 reactivity correlation: ")
	print(num_05)

	plt.hist(rvals, color='forestgreen', alpha=0.85)
	plt.show()

	log_covs = np.log(np.array(covs))
	plt.scatter(log_covs, rvals, color="black", label="All introns", s=10)
	if add_rpgs:
		(rpg_log_covs, rpg_rvals) = get_covs_rvals_names_list(log_covs, rvals, names, all_rpg_introns)
		plt.scatter(rpg_log_covs, rpg_rvals, color="red", label="RPG introns", s=10)
	if add_special_rpgs:
		(rpg_log_covs, rpg_rvals) = get_covs_rvals_names_list(log_covs, rvals, names, special_rpg_introns)
		plt.scatter(rpg_log_covs, rpg_rvals, color="gold", edgecolors='black', label="RPG introns putative structure", s=10)

	locs, labels = plt.xticks()
	new_labels = np.exp(np.array(locs))
	new_labels = [str(int(label)) for label in new_labels]
	plt.xticks(ticks=locs, labels=new_labels)
	plt.title('R value correlation between replicates vs total read coverage')
	plt.xlabel('Average coverage (reads per base)')
	plt.ylabel('R value')
	if add_rpgs or add_special_rpgs:
		plt.legend()
	plt.savefig('../figures/coverage_rval.png', format='png', dpi=300)
	# plt.show()

	log_coding_covs = np.log(np.array(coding_covs))
	plt.scatter(log_coding_covs, rvals)
	locs, labels = plt.xticks()
	new_labels = np.exp(np.array(locs))
	new_labels = [str(int(label)) for label in new_labels]
	plt.xticks(ticks=locs, labels=new_labels)
	plt.title('R value correlation between replicates vs coding region read coverage')
	plt.xlabel('Coverage (RPKM)')
	plt.ylabel('R value')
	plt.show()

def plot_ri_fraction():
	d3_cov = fill_coverage_dict(d3_stats_file, 'd3', intron_fasta_file, do_rpkm=True)
	d3_coding_cov = fill_coverage_dict(d3_coding_stats_file, 'd3', coding_fasta_file, do_rpkm=True)

	d45_cov = fill_coverage_dict(d45_stats_file, 'd45', intron_fasta_file, do_rpkm=True)
	d45_coding_cov = fill_coverage_dict(d45_coding_stats_file, 'd45', coding_fasta_file, do_rpkm=True)

	d3_ri = get_ri_fractions(d3_cov, d3_coding_cov)
	d45_ri = get_ri_fractions(d45_cov, d45_coding_cov)

	ri_fractions= []
	avg_cov = []
	for intron_key in d3_ri.keys():
		ri_fractions += [(d3_ri[intron_key] + d45_ri[intron_key])/2]
		avg_cov += [(d3_cov[intron_key] + d45_cov[intron_key])/2]

	plt.hist(ri_fractions, color='forestgreen', bins=200, alpha=0.85)
	plt.xlim((0, 2))
	plt.xticks(ticks=np.arange(0, 2, 0.1))
	plt.title('Histogram of RI fractions for all pre-mRNA')
	plt.xlabel('RI fraction')
	plt.ylabel('Number of introns')
	plt.show()

	print("Cov RI correlation:")
	print(pearsonr(np.array(avg_cov), np.array(ri_fractions)))

def plot_intron_vs_coding_coverage(intron_covs, coding_covs):
	r_val = pearsonr(np.array(intron_covs), np.array(coding_covs))
	print("Correlation between intron coverage and spliced mRNA coverage: %s\n" % str(r_val))

	plt.scatter(intron_covs, coding_covs, color="blue")
	plt.xlabel("Intron coverage (RPKM)")
	plt.ylabel("Spliced mRNA coverage (RPKM)")
	plt.title("Coverage: Intron vs Spliced mRNA")
	plt.show()

num_introns = get_max_num_introns()

plot_nodrug_drug_compare()
 
plot_ri_fraction()

[rvals, pvals, covs, covs_coding, names] = get_rvals_pvals_covs()
plot_intron_vs_coding_coverage(covs, covs_coding)
plot_coverage_rval(rvals, covs, covs_coding, names, add_rpgs=True, add_special_rpgs=False)

write_all_reactivities(outdir='combined_1221/reactivity/reactivity')
