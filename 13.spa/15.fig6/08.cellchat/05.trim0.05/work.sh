perl get.stat.all.pl Lymphoid_1_2_ICAR_sub_select19_trim0.05_LR.csv Lymphoid3_ICAR_sub_select19_trim0.05_LR.csv Lymphoid_1_2_ICAR_sub_select19_trim0.05_LR.csv > ICAR_sub_select_trim0.05.stat.xls

awk -F "," '($2~/Myeloid_Marco/ || $2~/Fib/) && ($3~/Treg/ || $3~/Lymphoid_B/)' Lymphoid_1_2_ICAR_sub_select19_trim0.05_LR.csv > mac.fib.treg.b.12.xls
awk -F "," '/Treg/ && /Lymphoid_B/' Lymphoid_1_2_ICAR_sub_select19_trim0.05_LR.csv >> mac.fib.treg.b.12.xls

perl get.stat.all.pl mac.fib.treg.b.12.xls mac.fib.treg.b.3.xls Lymphoid_1_2_ICAR_sub_select19_trim0.05_LR.csv > mac.fib.treg.b.123.stat.xls

