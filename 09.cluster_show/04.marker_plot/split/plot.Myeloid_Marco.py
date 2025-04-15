import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Myeloid_Marco.0905.h5ad")
marker_genes_dict = {
"Myeloid_Marco_C1QC":["C1QC","C1QB","APOC1"],
"Myeloid_Marco_SPP1":["SPP1","CCL4","CCL3L1","CXCR4"],
"Myeloid_Marco_LYVE1":["LYVE1","FOLR2","STAB1","F13A1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["Myeloid_Marco_C1QC","Myeloid_Marco_SPP1","Myeloid_Marco_LYVE1",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Myeloid_Marco.pdf")


