sc.pl.umap(adata_concat, color=[
'CCL8','GJA4','MYH11','MCAM','RGS5','IL6',
'COL5A1','COL5A2','COL6A3','DCN','FN1','LUM','POSTN','VCAN',
'CD74','HLA-DRA','HLA-DRB1','CCL21','CXCL12',
'KRT19','KRT8','SAA1','SLPI',
'APOA2','FABP1','FABP4','FRZB','GPX3',
'CXCL1','C3','C7','FBLN1','IGFBP6',
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".aaaa.bbbb_cccc_dddd_eeee.Fibroblast.umap.markergene.png")

sc.pl.umap(adata_concat, color=[
"POSTN","KRT19","IL6","FBLN2","CD74","KCNN3","SCN7A","P2RY1",
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".aaaa.bbbb_cccc_dddd_eeee.less.Fibroblast.umap.markergene.png")


