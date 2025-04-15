awk -F "," '$2~/Epithelium_Malig/ && ($3~/Marco/ || $3~/Treg/) && /MDK/' ../Lymphoid3_ICAR_trim0_LR.csv > mdk.3.csv
awk -F "," '$2~/Epithelium_Malig/ && ($3~/Marco/ || $3~/Treg/) && /MDK/' ../Lymphoid_1_2_ICAR_trim0_LR.csv > mdk.12.csv

perl get.stat.pl mdk.12.csv mdk.3.csv Lymphoid3up_vs_Lymphoid_1_2down_all_gene.csv > mdk.cci.csv

