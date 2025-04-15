awk '/Myeloid_pDC/ && /Epithelium_Malig_Migration/' ../Lymphoid3_ICAR_trim0_LR.csv  > lym3.csv
awk '/Myeloid_pDC/ && /Epithelium_Malig_Migration/' ../Lymphoid_1_2_ICAR_trim0_LR.csv > lym12.csv

perl get.stat.pl lym12.csv lym3.csv ../Lymphoid3up_vs_Lymphoid_1_2down_all_gene.csv > Myeloid_pDC_Epithelium_Malig_Migration.stat.xls

perl get.stat.pl lym12.csv lym3.csv diff.Myeloid_pDC_Lymphoid3.csv diff.Myeloid_pDC_Lymphoid_1_2.csv |les

