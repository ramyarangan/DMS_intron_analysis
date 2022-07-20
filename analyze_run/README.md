

Fig. 1F: Compare reactivity profiles (-pladB vs +pladB) for three introns with sufficient coverage in the -pladB experiment.
python compare_pladb_plus_minus_reactivity.py

Example output image (in figures/):

Example output r-vals: 
chrVII:555830-556307
(0.7567325317877045, 4.463571135865212e-40)
chrVII:311015-311526
(0.9273482810294933, 1.919830745169201e-101)
chrIV:579478-580017
(0.8354252119460018, 8.387530606708506e-70)


Fig. S1: Compare RI fraction +pladB vs -pladB for different classes: 
python analyze_pladb_ri_fraction_genes.py

Example output image (in figures/long_short_ri_fractions.png and figures/rpg_ri_fractions.png):


Example output stats: 
Mann Whitney U rank test p value:
2.4122316867070156e-21
Mann Whitney U rank test p value:
9.085459258551757e-14


Fig. S2A: DMS mutational frequency. 

Generate statistics used to make plots with:
python analyze_signal_noise.py d3



Fig. S2B: DMS validation using rRNA structures and reactivity
python plot_rrna_roc_curve.py

Example output images:

Example output stats:
Pearson correlation: 0.912067
R^2 value: 0.831865
AUC for ['18s', '25s']: 0.911435560847249



Fig. S2D: R-value between replicate reactivity vs intron coverage
Fig. 1E: +pladB vs -pladB RI fraction
python reactivity_correlation.py

Example output images:
In figures/coverage_rval.png: 


In figures/nodrug_drug_compare.png:



