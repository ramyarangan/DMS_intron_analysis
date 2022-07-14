from matplotlib import pyplot as plt 
from scipy import stats 
import pandas as pd 
import seaborn as sns 
import numpy as np

intron_fasta_file = "../intron_annot/standard_introns.fa"
ri_fractions_file = "pladb/ri_fractions.txt"
intron_stats_file = "pladb/intron_stats.txt"

COVERAGE_CUTOFF = 1000

def get_ri_fractions():
	f = open(ri_fractions_file)
	ri_fraction_lines = f.readlines()
	f.close()

	pladb_ri = {}
	nopladb_ri = {}
	coverage = {}
	for ri_fraction_line in ri_fraction_lines:
		ri_items = ri_fraction_line.split(' ')
		tag = ri_items[0][:-1]
		pladb_ri[tag] = float(ri_items[1])
		nopladb_ri[tag] = float(ri_items[2])
		coverage[tag] = float(ri_items[3])

	return pladb_ri, nopladb_ri, coverage

def get_pwm_scores_lengths():
	f = open(intron_stats_file)
	intron_stats_lines = f.readlines()
	f.close()

	pwms = {}
	lengths = {}

	for intron_stats_line in intron_stats_lines:
		intron_stats = intron_stats_line.split(" ")
		pwms[intron_stats[0]] = [float(x) for x in intron_stats[1:4]]
		lengths[intron_stats[0]] = [int(x) for x in intron_stats[4:]]

	return pwms, lengths

def get_rpg_introns():
	rpg_introns = {}

	f = open(intron_fasta_file)
	lines = f.readlines()
	f.close()

	for line in lines: 
		if len(line) > 0 and line[0] == '>':
			tag = line.split(" ")[0].replace(">", "")
			gene_name = line.split(' ')[1].replace('\n', '')
			is_rpg = False
			if (gene_name[:2] == "RP"):
				is_rpg = True
			elif (gene_name == "YML6") or (gene_name == "MRPL44"):
				is_rpg = True
			rpg_introns[tag] = is_rpg
	
	return rpg_introns

def compare_classes(vals_1, vals_2, val1_tag, val2_tag, ylabel, plt_title="", 
	savefig=False, filename=""):
	vals_1 = np.array(vals_1)
	vals_2 = np.array(vals_2)
	stat, pval = stats.mannwhitneyu(vals_1, vals_2)
	print("Mann Whitney U rank test p value:")
	print(pval)

	df = pd.DataFrame(columns=[val1_tag, val2_tag])

	vals1_added = [np.nan] * max(len(vals_1), len(vals_2))
	vals1_added[:len(vals_1)] = vals_1 
	df[val1_tag] = vals1_added
	vals2_added = [np.nan] * max(len(vals_1), len(vals_2))
	vals2_added[:len(vals_2)] = vals_2
	df[val2_tag] = vals2_added
	plt.figure(figsize=(4, 3))
	my_pal = {val1_tag: "purple", val2_tag: "darkorange"}
	ax = sns.violinplot(data=df, palette=my_pal, cut=0)
	plt.ylabel(ylabel)
	plt.title(plt_title)
	if savefig:
		plt.savefig(filename, format='png', dpi=300)
	else: 
		plt.show()

def plot_scatter(vals_1, vals_2, tag, xlabel, ylabel):
	r = stats.pearsonr(np.array(vals_1), np.array(vals_2))
	print("Correlation for %s:" % tag)
	print(r)

	plt.scatter(vals_1, vals_2, color='black')
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(tag)
	plt.show()

# RP genes RI Fraction increase vs the rest
def compare_rpg_ris():
	rpg_introns = get_rpg_introns()
	pladb_ri, nopladb_ri, coverage = get_ri_fractions()

	rpg_ris = []
	rpg_ri_ratios = []
	non_rpg_ris = []
	non_rpg_ri_ratios = []

	for cur_key in pladb_ri.keys():
		if coverage[cur_key] < COVERAGE_CUTOFF: 
			continue

		ri_ratio = pladb_ri[cur_key]/nopladb_ri[cur_key]
		if rpg_introns[cur_key]:
			rpg_ris += [pladb_ri[cur_key]]
			rpg_ri_ratios += [ri_ratio]
		else:
			non_rpg_ris += [pladb_ri[cur_key]]
			non_rpg_ri_ratios += [ri_ratio]

	compare_classes(rpg_ris, non_rpg_ris, "RPG introns", \
		"Other introns", "RI +pladB", "RPG vs non-RPG retained intron fractions +pladB")
	compare_classes(rpg_ri_ratios, non_rpg_ri_ratios, "RPG introns", \
		"Other introns", "(RI +pladB)/(RI -pladB)", \
		"RPG vs non-RPG retained intron fraction ratios", \
		savefig=True, filename="../figures/rpg_ri_fractions.png")

def compare_metrics_to_ri(metric_dict, metric_labels): 
	pladb_ri, nopladb_ri, coverage = get_ri_fractions()

	for ii, metric_label in enumerate(metric_labels):
		metric_vals = []
		ris = []
		ri_ratios = []

		for cur_key in pladb_ri.keys():
			if coverage[cur_key] < COVERAGE_CUTOFF: 
				continue

			ris += [pladb_ri[cur_key]]
			ri_ratios += [pladb_ri[cur_key]/nopladb_ri[cur_key]]

			metric_vals += [metric_dict[cur_key][ii]]

		plot_scatter(metric_vals, ris, metric_label + "vs RI", \
			metric_label, "RI")
		plot_scatter(metric_vals, ri_ratios, \
			metric_label + "vs RI ratio", metric_label, "RI ratio")

# PWM agreement with splice site sequences
def compare_pwm_ri_fraction():
	pwms, _ = get_pwm_scores_lengths()
	pwm_labels = ["5'SS", "BP", "3'SS"]
	compare_metrics_to_ri(pwms, pwm_labels)

# Length distributions
def compare_len_ri_fraction():
	_, lengths = get_pwm_scores_lengths()
	len_labels = ["Total Length", "5'SS-BP Length", "BP-3'SS Length"]
	compare_metrics_to_ri(lengths, len_labels)

def compare_long_short_ris():
	_, lengths = get_pwm_scores_lengths()
	pladb_ri, nopladb_ri, coverage = get_ri_fractions()

	long_ris = []
	long_ri_ratios = []
	short_ris = []
	short_ri_ratios = []

	for cur_key in pladb_ri.keys():
		if coverage[cur_key] < COVERAGE_CUTOFF: 
			continue

		ri_ratio = pladb_ri[cur_key]/nopladb_ri[cur_key]
		# Long introns have total length > 200
		if lengths[cur_key][0] > 200:
			long_ris += [pladb_ri[cur_key]]
			long_ri_ratios += [ri_ratio]
		else:
			short_ris += [pladb_ri[cur_key]]
			short_ri_ratios += [ri_ratio]

	compare_classes(long_ris, short_ris, "Long introns", \
		"Short introns", "RI +pladB", "Long vs short retained intron fractions +pladB")
	compare_classes(long_ri_ratios, short_ri_ratios, "Long introns", \
		"Short introns", "(RI +pladB)/(RI -pladB)", \
		"Long vs short retained intron fraction ratios", \
		savefig=True, filename="../figures/long_short_ri_fractions.png")

compare_rpg_ris()
# compare_pwm_ri_fraction()
# compare_len_ri_fraction()
compare_long_short_ris()
