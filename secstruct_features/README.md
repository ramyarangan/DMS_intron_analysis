### **Fig. 4E-H**

Comparing intron and control sequences' reactivity profiles and secondary structures. The maximum extrusion from ends metric is cached for all introns in the folder `mee_cache/`; all other metrics are quick to evaluate and computed each time the plot is run.

Run: `python compare_secstruct_features.py ../intron_annot/standard_introns.fa rnastructure_sherlock_1221/intron/ ../analyze_run/combined_1221/reactivity/reactivity_intron/ ../analyze_run/combined_1221/rfcount/d_all_dedup_view.txt mee_cache/intron_mee_1221.txt ../intron_annot/standard_introns.dat ../intron_annot/coding_orfs_introns.fa rnastructure_sherlock_1221/coding_nd/ ../analyze_run/combined_1221/reactivity/reactivity_coding_nd/ ../analyze_run/combined_1221/rfcount/nd_coding_dedup_view.txt mee_cache/coding_mee_1221.txt`

The parameters here point to the following files for intron and coding sequences: sequence files, secondary structure prediction directories, reactivity values, per-base sequencing coverage, cached maximum extrusion from ends. For introns, we also pass in an intron database file including branchpoint locations for each intron.

Example figures:
Fig. 4E (output to `../figures/Gini.png`): 

Fig. 4F (output to `../figures/MEE_normalized_violinplot.png`):

Fig. 4G (output to `../figures/longest_stem_len_compare.png`):

Fig 4H (output to `../figures/bpp_compare.png`):


Example output statistics comparing intron vs coding values. Also includes the number of zipper stems and end stems:
```
Mann Whitney U rank test p value comparing intron / coding longest stem lengths:
0.0005278018928243133
Number of intron stems passing cutoff: 83
Number of coding stems passing cutoff: 83
Mann Whitney U rank test p value comparing intron / coding base-pair probability in stems:
1.6980929445123513e-14
Number of intron stems passing cutoff: 380
Number of coding stems passing cutoff: 430
Mann Whitney U rank test p value comparing intron / coding Gini coeffs:
1.0948955617784742e-20
Number of intron stems passing cutoff: 437
Number of coding stems passing cutoff: 1721

Number of zipper stems: 42

Number of end stems: 30

Mann Whitney U rank test p value comparing normalized MEE:
4.506242697742745e-54
Number of intron stems passing cutoff: 276
Number of coding stems passing cutoff: 178
```


### **Fig. 5A**

Heatmap summarizing intron secondary structure features, with hierarchical clustering and dendrogram visualization.

Run: `python plot_write_features.py --heatmap`
This function also outputs the introns in the order of the dendrogram to `dendrogram/dendrogram_order.txt`.

Example figure:

Example output class popoulations: 
```
Class 1 size: 12, 0.077922
Class 2 size: 16, 0.103896
Class 3 size: 14, 0.090909
Class 4 size: 34, 0.220779
Class 5 size: 16, 0.103896
Class 6 size: 9, 0.058442
Class 7 size: 12, 0.077922
Class 8 size: 23, 0.149351
Class 9 size: 18, 0.116883
Total number of rows: 154
```


### **Fig. 5B**

tSNE plot clustering introns by secondary structure features. Note that the clustering is expected to be different for each run.

Run: `python plot_write_features.py --tsne`

Example figure (output to `../figures/tSNE_introns.png`): 


### **Table S2**

Table summarizing secondary structures and structural featuers for introns.

Run: `python plot_write_features.py --write_table --csv_table_filename secstruct_features.csv`

Table is output in `secstruct_features.csv`. 

Output data: 
```
Number of long introns with DMS data: 113
Number of zipper stems: 42
Number of downstream stems: 30
Number of introns with either zipper or downstream stem: 60
```

