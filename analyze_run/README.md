### **Fig. 1F**
Compare reactivity profiles (-pladB vs +pladB) for three introns with sufficient coverage in the -pladB experiment.

Run: `python compare_pladb_plus_minus_reactivity.py`

Example output image (in `figures/`):
![chrIV:579478-580017_d](https://user-images.githubusercontent.com/2606810/179919448-ce475eac-d9ce-4694-bd23-96e820267cce.png)
![chrIV:579478-580017_nd](https://user-images.githubusercontent.com/2606810/179919451-4988da59-dbb9-450f-a742-ec9aaf27f435.png)

Example output r-vals: 
```
chrVII:555830-556307
(0.7567325317877045, 4.463571135865212e-40)
chrVII:311015-311526
(0.9273482810294933, 1.919830745169201e-101)
chrIV:579478-580017
(0.8354252119460018, 8.387530606708506e-70)
```

### **Fig. S1** 
Compare retained intron (RI) fraction +pladB vs -pladB for different intron classes.

Run: `python analyze_pladb_ri_fraction_genes.py`

Example output image (`figures/long_short_ri_fractions.png` and `figures/rpg_ri_fractions.png`):

![long_short_ri_fractions](https://user-images.githubusercontent.com/2606810/179919589-93214e82-5779-44f8-881f-c27eef0cd01a.png)

Example output stats: 
```
Mann Whitney U rank test p value:
2.4122316867070156e-21
Mann Whitney U rank test p value:
9.085459258551757e-14
```

### **Fig. S2A** 
Analysis of per-base DMS mutational frequency. 

To generate statistics used for plot, run: `python analyze_signal_noise.py d3`

### **Fig. S2B** 
DMS validation using rRNA structures and reactivity.

Run: `python plot_rrna_roc_curve.py`

Example output images:

<img width="556" alt="Screen Shot 2022-07-20 at 12 01 23 AM" src="https://user-images.githubusercontent.com/2606810/179919946-daaca60c-620f-4b51-82a8-f0ddbfe58719.png">

<img width="551" alt="Screen Shot 2022-07-20 at 12 01 33 AM" src="https://user-images.githubusercontent.com/2606810/179919896-34852e30-35b1-4c54-b0b8-d1bb93ed8cee.png">

Example output stats:
```
Pearson correlation: 0.912067
R^2 value: 0.831865
AUC for ['18s', '25s']: 0.911435560847249
```


### **Fig. S2D, Fig. 1E** 
R-value between replicate reactivity vs intron coverage, and +pladB vs -pladB retained intron (RI) fraction

Run: `python reactivity_correlation.py`

Example output images:
In `figures/coverage_rval.png`: 

![coverage_rval](https://user-images.githubusercontent.com/2606810/179920121-0c2d8345-7e54-4bab-a97a-a984850becad.png)

In `figures/nodrug_drug_compare.png`:

![nodrug_drug_compare](https://user-images.githubusercontent.com/2606810/179920145-af9b2891-886f-41a1-80a8-3729d9ee2ef7.png)


