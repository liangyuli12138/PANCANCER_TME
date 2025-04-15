import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Mural_cell.0905.h5ad")
marker_genes_dict = {
"Mural_cell_Pericyte1":["PDGFRB","MYO1B","ARHGDIB","HIGD1B","NOTCH3","RGS5","CYGB","SPARC","NDUFA4L2","COL4A1"],
"Mural_cell_Pericyte2":["CD44","IGFBP5","GGT5","CCL2","CLSTN2","CXCL12"],
"Mural_cell_SMC1":["ADIRF","MUSTN1","RERGL","PLN","SORBS2","CRIP1","NET1"],
"Mural_cell_SMC2":["RAMP1","PALLD","MYLK","TCEAL4","NCAM1","COLEC12","PSAP","SFRP4"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["Mural_cell_Pericyte1","Mural_cell_Pericyte2","Mural_cell_SMC1","Mural_cell_SMC2",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Mural_cell.pdf")


