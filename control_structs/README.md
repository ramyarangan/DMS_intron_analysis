### **Fig. 2A-C, Table S1**

Evaluation of DMS-guided structure predictions for control RNA structures.

Run: `evaluate_confidence_estimates.py`

Predicted secondary structures and base-pairing probabilities are included in `control_secstruct_bpps/`, and native secondary structures are in `native_secstructs`. 

Example figure (output to `../figures/control_ppv_sens_f1.png`): 

![control_ppv_sens_f1](https://user-images.githubusercontent.com/2606810/179921863-ddbea04a-ef05-4909-a047-1e6b7f717270.png)

Example output:
```
Vienna MFE prediction
Total confusion matrix
[57 50 19  0  0]
Precision: 0.532710
Recall: 0.750000
F1 Score: 0.622951

DMS-guided RNAstructure, no helix confidence estimate cutoff
Total confusion matrix
[57 34 21  0  0]
Precision: 0.626374
Recall: 0.730769
F1 Score: 0.674556

DMS-guided RNAstructure, helix confidence estimate cutoff 0.7
1800
RDN18-1
(14, 5, 15, 0, 4)
121
RDN5-1
(4, 0, 0, 0, 0)
158
RDN58-1
(1, 0, 0, 0, 0)
214
SNR7-L
(4, 3, 2, 0, 0)
568
SNR19
(18, 2, 4, 0, 0)
382
HAC1
(1, 1, 0, 0, 0)
356
ASH1
(2, 0, 0, 0, 0)
351
RPS28B
(2, 0, 0, 0, 0)
301
SFT2
(0, 0, 1, 0, 0)
72
TRR4
(1, 0, 2, 0, 0)
72
TRT2
(2, 0, 1, 0, 1)
72
IMT4
(2, 0, 1, 0, 0)
Total confusion matrix
[51 11 26  0  5]
Total confusion matrix mRNA
[5 1 1 0 0]
Precision: 0.822581
Recall: 0.662338
F1 score: 0.733813
```
