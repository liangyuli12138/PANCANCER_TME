import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/EC_Alveolar");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.EC.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["FENDRR","VIPR1","MYZAP","SLC6A4","HPGD","CD82","TMEM100","ACE","TNFSF10","ADGRL2","PAPSS2","STXBP6","EPAS1","TBX2","PRX","ARHGAP18","RALA","EDNRB","EMP2","ADRB1","GPX3","CAV2","CR381670.1","CA4","S100A4","AFF3","ACVRL1","FOXF1","PRKCE","SULT1A1","IL7R","SFTPC","TSPAN12","EPB41L2","DPYSL3","ZNF331","CARD16","RCSD1","PDE1C","PRICKLE2","KHDRBS2","OCLN","LMO7","SOSTDC1","TBX3","AC104984.4","S100A3","ITGA1","CAV1","JCAD",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Alveolar",
                  save="EC_Alveolar.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Alveolar.markergene.0.xxx.pdf")

marker_genes_dict=["PDE3B","MYLK","VAT1","CYB5A","ALDH2","NDRG4","CASZ1","PLXNC1","MEIS1","CYB561","NEDD4L","MARVELD2","HLA-E","AC019257.8","B3GALNT1","RAMP2","IL32","PHLDB1","AC138776.1","ANXA3","HECW2","CYP3A5","DISP1","ITGA4","XKR6","CMTM3","PCDH17","MAOA","STK4","NCALD","GCOM1","ITGA3","MAP4K2","SFTPB","SLPI","ITM2B","MFSD6","TSPAN15","DACH1","PASK","FRY","IFNGR1","CAVIN1","GBP4","ANKRD20A17P","CRMP1","TUSC3","SPOCK2","DAAM1","EXOSC7",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Alveolar",
                  save="EC_Alveolar.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Alveolar.markergene.1.xxx.pdf")
