import sys
import argparse
import pandas as pd 

from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import seaborn as sns 
from matplotlib import pyplot as plt

from secstruct_util import *
from get_secstruct_features import *
from zipper_stem_feature import *
from accessibility_feature import *
from get_gini_features import *

parser = argparse.ArgumentParser(description='Parameters for secondary structure feature summary data')
parser.add_argument('--write_table', default=False, action='store_true', \
	help='Write secondary structure feature table to csv file')
parser.add_argument('--csv_table_filename', type=str, help='Output file for secondary structure feature table')
parser.add_argument('--heatmap', default=False, action='store_true', \
	help='Hierarchical clustering and display heatmap for secondary structure features')
parser.add_argument('--pca', default=False, action='store_true', help='PCA cluster and display')
parser.add_argument('--tsne', default=False, action='store_true', help='tSNE cluster and display')
args = parser.parse_args()

do_heatmap = args.heatmap
do_pca = args.pca
do_tsne = args.tsne
write_table = args.write_table
csv_filename = ""
if write_table:
	csv_filename = args.csv_table_filename

fasta_file = "../intron_annot/standard_introns.fa"
secstruct_dir = "rnastructure_sherlock_1221/intron/"
dat_file = "../intron_annot/standard_introns.dat"
reac_cov_dir = "../analyze_run/combined_1221/reactivity/reactivity_intron/"
mpl_stats_file = "mpl_cache/intron_1221.txt"
# mpl_stats_file = ""
extended_dat_file = "../intron_annot/standard_allsize_extend50_baseinfo.dat"
extended_reac_dir = "../analyze_run/combined_1221/reactivity/reactivity_allsize_extend50/"
mut_freq_file = '../analyze_run/combined_1221/rfcount/d_all_dedup_view.txt'

COVERAGE_CUTOFF = 1971
OVERHANG_5PRIME = 50
OVERHANG_3PRIME = 50

def get_features_df():
	columns = ["Length", "hasZipper", "ZipperdG", "hasEnd", "EndStemdG", "normMLD", "longestStem", \
		"Average BPP", "Max Gini", "5'SSaccess", "BPaccess", "3'SSaccess"] 
	feature_df = pd.DataFrame(columns=columns)
	
	names, _, _, name_symbol_dict = get_names_seqs_from_fasta(fasta_file)

	data = []
	name_list = []
	for name in names:
		feature_vals = []

		intron_len = get_length(name, fasta_file)

		mpl = get_mpl_from_file(name, mpl_stats_file, reac_cov_dir)
		
		feature_vals += [mpl * intron_len]
		
		longest_stem_len = get_longest_stem_len(name, fasta_file, secstruct_dir, reac_cov_dir)

		feature_vals += [longest_stem_len]
		bpps = get_bpps_tag(name, fasta_file, secstruct_dir, reac_cov_dir)

		avg_bpp = -1
		if bpps == None or len(bpps) == 0:
			avg_bpp = 0
		else:
			avg_bpp = sum(bpps)/len(bpps)
		feature_vals += [avg_bpp]

		gini_vals = get_gini_windows_tag(name, fasta_file, mut_freq_file)
		if len(gini_vals) == 0:
			gini_vals = [0]
		feature_vals += [max(gini_vals)]

		ss_reac = get_reac_bp_from_tag(name, extended_dat_file, extended_reac_dir)
		feature_vals += list(ss_reac)

		feature_vals = np.array(feature_vals)
		if sum(np.isnan(feature_vals)) > 0:
			continue
		if sum(feature_vals < 0) > 0:
			continue
		feature_vals = list(feature_vals)

		_, _, _, _, _, _, dG = get_best_zipper_stem_tag(name, fasta_file, secstruct_dir, \
			dat_file, do_second=True)
		if dG == "" or float(dG) > 0:
			feature_vals = [0.0, 0] + feature_vals
		else:
			feature_vals = [1.0, -float(dG)] + feature_vals

		_, _, _, _, _, _, dG = get_best_zipper_stem_tag(name, fasta_file, secstruct_dir, dat_file)
		if dG == "" or float(dG) > 0:
			feature_vals = [0.0, 0] + feature_vals
		else:
			feature_vals = [1.0, -float(dG)]+ feature_vals

		feature_vals = [float(intron_len)] + feature_vals

		val_dict = dict(zip(columns, feature_vals))
		data.append(val_dict)
		name_list += [name]

	# Assemble normalized feature values for columns
	feature_df = feature_df.append(data, True)
	feature_df = (feature_df - feature_df.min())/(feature_df.max() - feature_df.min())
	return feature_df, name_list, name_symbol_dict

def add_stem_class_feature(feature_df):
	stem_class_num = feature_df["hasZipper"].values + \
		feature_df["hasEnd"].values * 2
	stem_class = np.array([""] * len(stem_class_num), dtype=object)
	stem_class[stem_class_num == 0] = "No zipper or end stem"
	stem_class[stem_class_num == 1] = "Zipper stem"
	stem_class[stem_class_num == 2] = "End stem"
	stem_class[stem_class_num == 3] = "Zipper and end stems"

	feature_df["StemClass"] = stem_class

	return feature_df

def plot_components(feature_df, col1, col2, tag1, tag2, plt_name):
	fig, (ax1,ax2) = plt.subplots(1, 2, constrained_layout=True, figsize=(15,6))


	feature_df = add_stem_class_feature(feature_df)
	p2 = sns.scatterplot(x=col1, y=col2, data=feature_df, \
		hue='StemClass', palette=sns.color_palette("Set2")[:4], s=15, ax=ax1, edgecolor='gray')
	p2.set_xlabel(tag1)
	p2.set_ylabel(tag2)

	p1 = sns.scatterplot(x=col1, y=col2, data=feature_df, \
		hue='Length', palette='BuPu', s=15, ax=ax2, edgecolor='gray')
	# alpha=0.8, - add back in to change transparency
	p1.set_xlabel(tag1)
	p1.set_ylabel(tag2)

	norm = plt.Normalize(feature_df["Length"].min(), \
		feature_df["Length"].max())
	sm = plt.cm.ScalarMappable(cmap="BuPu", norm=norm)
	sm.set_array([])

	ax2.get_legend().remove()
	cbar = ax2.figure.colorbar(sm, ticks=[0, 1])
	cbar.set_label('Normalized length', rotation=270)

	plt.savefig("../figures/" + plt_name + ".png", format='png', dpi=300)

def make_tsne():
	feature_df, _, _ = get_features_df()
	feature_array = np.asarray(feature_df.values)

	tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
	tsne_results = tsne.fit_transform(feature_df)

	feature_df['tsne-2d-one'] = tsne_results[:,0]
	feature_df['tsne-2d-two'] = tsne_results[:,1]

	plot_components(feature_df, 'tsne-2d-one', 'tsne-2d-two', \
		'tSNE1', 'tSNE2', 'tSNE_introns')

def make_pca():
	feature_df, _, _ = get_features_df()
	feature_array = np.asarray(feature_df)

	pca = PCA(n_components=3)
	pca_result = pca.fit_transform(feature_df.values)

	feature_df['pca-one'] = pca_result[:,0]
	feature_df['pca-two'] = pca_result[:,1]
	feature_df['pca-three'] = pca_result[:,2]
	
	print_str = 'Explained variation per principal component: {}'
	print(print_str.format(pca.explained_variance_ratio_))

	plot_components(feature_df, 'pca-one', 'pca-two', \
		'PCA 1', 'PCA 2', 'PCA_introns')

def make_heatmap():
	feature_df, names, name_symbol_dict = get_features_df()
	feature_array = np.asarray(feature_df)

	# Get class labels and colors based on hierarchical clustering
	row_linkage = hierarchy.linkage(distance.pdist(feature_array), method='ward')
	labels = fcluster(row_linkage, 5, criterion='maxclust')
	labels = np.array(labels)
	print(labels)
	palette = sns.color_palette("deep")
	palette = [palette[0], palette[1], palette[8] ,palette[6], palette[3]]
	row_colors = [palette[label - 1] for label in labels]
	print("Class 1 size: %d, %f" % (np.sum(labels == 1), np.sum(labels == 1)/len(list(labels))))
	print("Class 2 size: %d, %f" % (np.sum(labels == 2), np.sum(labels == 2)/len(list(labels))))
	print("Class 3 size: %d, %f" % (np.sum(labels == 3), np.sum(labels == 3)/len(list(labels))))
	print("Class 4 size: %d, %f" % (np.sum(labels == 4), np.sum(labels == 4)/len(list(labels))))
	print("Class 5 size: %d, %f" % (np.sum(labels == 5), np.sum(labels == 5)/len(list(labels))))

	# Plot dendogram with clusters shown
	# pd.set_option("display.max_rows", None, "display.max_columns", None)
	g = sns.clustermap(feature_df, row_linkage=row_linkage, method='ward', \
		col_cluster=False, cmap="BuGn", row_colors=row_colors, cbar_kws={"ticks": [], "shrink": 0.5})

	f = open("dendrogram_order.txt", 'w')
	for ii in g.dendrogram_row.reordered_ind:
		f.write("%s\t%s\n" % (names[ii], name_symbol_dict[names[ii]]))
	f.close()

	plt.show()


def write_table_csv(csv_filename):
	names, _, _, name_symbol_dict = get_names_seqs_from_fasta(fasta_file)

	f = open(csv_filename, 'w')
	
	num_start_zipper = 0
	num_end_zipper = 0
	num_either_stem = 0
	num_long = 0

	heading_str = "Gene Symbol,Coordinates,Length,Coverage, Zipper stem strand 1 start, Zipper stem strand 1 end, "
	heading_str += "Zipper stem strand 2 start, Zipper stem strand 2 end, Zipper stem dG (kcal/mol), "
	heading_str += "BP-3'SS stem strand 1 start, BP-3'SS stem strand 1 end, "
	heading_str += "BP-3'SS stem strand 2 start, BP-3'SS stem strand 2 end, BP-3'SS stem dG (kcal/mol), "
	heading_str += "Normalized MLD, Longest stem, 5'SS accessibility, BP accessibility, 3'SS accessibility"

	f.write("%s\n" % heading_str)

	for name in names:
		csv_str = name_symbol_dict[name] + "," + name
		
		intron_len = get_length(name, fasta_file)
		csv_str += "," + str(intron_len)
		
		cov = get_cov(name, reac_cov_dir)
		csv_str += "," + str(cov)
		
		if cov < COVERAGE_CUTOFF:
			f.write("%s\n" % csv_str)
			continue

		if intron_len > 200:
			num_long += 1

		has_either_stem = False
		zipper_items = get_best_zipper_stem_tag(name, fasta_file, secstruct_dir, dat_file)
		csv_str += "," + ",".join(zipper_items[2:])

		if zipper_items[2] != "":
			num_start_zipper += 1
			has_either_stem = True

		zipper_items = get_best_zipper_stem_tag(name, fasta_file, secstruct_dir, dat_file, do_second=True)
		csv_str += "," + ",".join(zipper_items[2:])
		if zipper_items[2] != "":
			num_end_zipper += 1
			has_either_stem = True

		if has_either_stem:
			num_either_stem += 1

		mpl = get_mpl_from_file(name, mpl_stats_file, reac_cov_dir)
		csv_str += ","
		if mpl > 0:
			csv_str += '{0:.3f}'.format(mpl)

		longest_stem_len = get_longest_stem_len(name, fasta_file, secstruct_dir, reac_cov_dir)
		csv_str += ","
		if longest_stem_len > 0:
			csv_str += str(longest_stem_len)

		ss_reac = get_reac_bp_from_tag(name, extended_dat_file, extended_reac_dir)
		if ss_reac[0] > -1:
			ss_reac = ['{0:.3f}'.format(x) for x in ss_reac]
			csv_str	+= "," + ",".join(ss_reac)
		else:
			csv_str += ',,,'

		f.write("%s\n" % csv_str)

	f.close()
	print("Number of long introns with DMS data: %d" % num_long)
	print("Number of zipper stems: %d" % num_start_zipper)
	print("Number of stems between BP and 3'SS: %d" % num_end_zipper)
	print("Number of introns with either zipper or end stem: %d" % num_either_stem)

if write_table:
	write_table_csv(csv_filename)
if do_heatmap:
	make_heatmap()
if do_pca:
	make_pca()
if do_tsne:
	make_tsne()
