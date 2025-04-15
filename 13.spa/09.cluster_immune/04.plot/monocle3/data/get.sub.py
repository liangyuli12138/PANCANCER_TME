import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import harmonypy as hm
import pandas as pd


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/D01972C4/D01972C4_cellbin.final.h5ad")

cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/monocle3/data/index.Lymphoid_B_memory")
bdata_concat = adata_concat[cellist["cell"],:]
bdata_concat.write_h5ad("index.Lymphoid_B_memory.h5ad")
bdata_concat.obs.to_csv("index.Lymphoid_B_memory.obs")

cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/monocle3/data/index.Lymphoid_CD4_Treg")
bdata_concat = adata_concat[cellist["cell"],:]
bdata_concat.write_h5ad("index.Lymphoid_CD4_Treg.h5ad")
bdata_concat.obs.to_csv("index.Lymphoid_CD4_Treg.obs")

cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/monocle3/data/index.Lymphoid_CD8_Tex")
bdata_concat = adata_concat[cellist["cell"],:]
bdata_concat.write_h5ad("index.Lymphoid_CD8_Tex.h5ad")
bdata_concat.obs.to_csv("index.Lymphoid_CD8_Tex.obs")

cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/monocle3/data/index.Myeloid_cDC2")
bdata_concat = adata_concat[cellist["cell"],:]
bdata_concat.write_h5ad("index.Myeloid_cDC2.h5ad")
bdata_concat.obs.to_csv("index.Myeloid_cDC2.obs")

cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/monocle3/data/index.Myeloid_Marco_C1QC")
bdata_concat = adata_concat[cellist["cell"],:]
bdata_concat.write_h5ad("index.Myeloid_Marco_C1QC.h5ad")
bdata_concat.obs.to_csv("index.Myeloid_Marco_C1QC.obs")


