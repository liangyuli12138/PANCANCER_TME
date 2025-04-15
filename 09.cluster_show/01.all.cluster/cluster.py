import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/01.all.cluster")
adata = sc.read_h5ad("pancancer.ref.0723.final.h5ad")

cellist = pd.read_csv("merge.all.all.level.filterB.sev.input")
adata = adata[cellist["cell"],:]

atlist = pd.read_csv("merge.all.all.level.filterB.sev.at",index_col=0)
adata.obs = adata.obs.join(atlist)

adata.write_h5ad("pancancer.ref.0905.final.h5ad")
adata.obs.to_csv("pancancer.ref.0905.final.obs.csv")
adata.var.to_csv("pancancer.ref.0905.final.var.csv")

palettedir3={
"Lymphoid_Plamsa":"#1F77B4",
"Myeloid":"#AEC7E8",
"Lymphoid_T_NK":"#FF7F0E",
"EC":"#FFBB78",
"Fibroblast":"#2CA02C",
"Lymphoid_B":"#98DF8A",
"Mast":"#D62728",
"Epithelium":"#FF9896",
"Mural_cell":"#C5B0D5",
}

palettedir1={
"Lymphoid_CD4":"#2CA02C",
"Lymphoid_CD8":"#FF7F0E",
"Lymphoid_ILC":"#1F77B4",
"Lymphoid_MAIT":"#AEC7E8",
"Lymphoid_NK":"#FFBB78",
"Lymphoid_NKT":"#98DF8A",
"Lymphoid_B":"#C49C94",
"Lymphoid_Plamsa":"#B85E17",
"Myeloid_Mono":"#C5B0D5",
"Myeloid_DC":"#9467BD",
"Myeloid_Marco":"#06A0FC",
"Myeloid_Mast":"#BCBD22",
"Fibroblast_mCAF":"#E377C2",
"Fibroblast_iCAF":"#DBDB8D",
"Fibroblast_apCAF":"#32FE04",
"Mural_Pericyte":"#FF9896",
"Mural_SMC":"#17BECF",
"EC":"#F7B6D2",
"Epithelium_Malig":"#D62728",
"Epithelium_Normal":"#9EDAE5",
}

palettedir2={
"Lymphoid_CD4_Tn":"#AEC7E8",
"Lymphoid_CD4_Tcm":"#808000",
"Lymphoid_CD4_Tfh":"#88DBDB",
"Lymphoid_CD4_Th17":"#C0FB04",
"Lymphoid_CD4_CTL":"#FD7003",
"Lymphoid_CD4_Treg":"#E377C2",
"Lymphoid_CD4_Tstr":"#C49C94",
"Lymphoid_CD8_Tn":"#37B01C",
"Lymphoid_CD8_Tm":"#89ACC4",
"Lymphoid_CD8_Teff":"#E8FE9D",
"Lymphoid_CD8_Tex":"#A3E70C",
"Lymphoid_CD8_Tstr":"#F66B6B",
"Lymphoid_CD8_Tisg":"#05FEB2",
"Lymphoid_ILC":"#237311",
"Lymphoid_MAIT":"#8B64F1",
"Lymphoid_NK_CD16+":"#F9A767",
"Lymphoid_NK_CD56+":"#4904FE",
"Lymphoid_NKT":"#B85E17",
"Lymphoid_B_naive":"#A9F4AE",
"Lymphoid_B_memory":"#075103",
"Lymphoid_Plamsa_IGLC":"#BCBD22",
"Lymphoid_Plamsa_IGKC":"#065A8C",
"Myeloid_Mono":"#DEA65E",
"Myeloid_cDC1":"#691768",
"Myeloid_cDC2":"#97B732",
"Myeloid_cDC3":"#F305FD",
"Myeloid_pDC":"#000075",
"Myeloid_Marco_C1QC":"#B62473",
"Myeloid_Marco_SPP1":"#05FEFE",
"Myeloid_Marco_LYVE1":"#FFD8B1",
"Myeloid_Mast":"#D62728",
"Fibroblast_mCAF_POSTN":"#FCC405",
"Fibroblast_mCAF_KRT19":"#06A0FC",
"Fibroblast_mCAF_WNT5A":"#53E27D",
"Fibroblast_iCAF_IGFBP6":"#F59DA5",
"Fibroblast_iCAF_IL6":"#864EAC",
"Fibroblast_iCAF_KCNN3":"#04A597",
"Fibroblast_apCAF_CD74":"#8B1E3E",
"Mural_cell_Pericyte1":"#6DC3F7",
"Mural_cell_Pericyte2":"#7DF8D2",
"Mural_cell_SMC1":"#BC68F4",
"Mural_cell_SMC2":"#FABED4",
"EC_Angiogenic":"#32FE04",
"EC_Artery":"#F80387",
"EC_Capillary":"#17BECF",
"EC_Vein":"#F96A94",
"EC_Lymph":"#97989A",
"EC_Alveolar":"#F4D260",
"EC_Glomerular":"#9905FD",
"EC_Sinusoidal":"#DCBEFF",
"Epithelium_Malig_Migration":"#FD0750",
"Epithelium_Malig_Cycle":"#DBDB8D",
"Epithelium_Malig_cEMT":"#A516AC",
"Epithelium_Malig_Interferon":"#F4F804",
"Epithelium_Malig_Stress":"#BF9C25",
"Epithelium_Malig_Basal":"#F791FC",
"Epithelium_Malig_Glandular":"#93FD7B",
"Epithelium_Malig_Multi":"#526DFA",
"Epithelium_Normal":"#9EDAE5",
}

sc.pl.umap(adata, color="groups_pri",palette = palettedir2,size=0.3, sort_order= "False", save= ".groups_pri.cluster.0.5.png")
sc.pl.umap(adata, color="groups_sce",size=0.3, sort_order= "False", save= ".groups_sce.cluster.0.5.png")
sc.pl.umap(adata, color="groups_thi",palette = palettedir1,size=0.3, sort_order= "False", save= ".groups_thi.cluster.0.5.png")
sc.pl.umap(adata, color="groups_fif",size=0.3, sort_order= "False", save= ".groups_fif.cluster.0.5.png")
sc.pl.umap(adata, color="groups_six",size=0.3, sort_order= "False", save= ".groups_six.cluster.0.5.png")
sc.pl.umap(adata, color="groups_sev",size=0.3, sort_order= "False", save= ".groups_sev.cluster.0.5.png")

sc.pl.umap(adata, color="groups_pri",palette = palettedir2,size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_pri.cluster.0.5.png")
sc.pl.umap(adata, color="groups_sce",size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_sce.cluster.0.5.png")
sc.pl.umap(adata, color="groups_thi",palette = palettedir1,size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_thi.cluster.0.5.png")
sc.pl.umap(adata, color="groups_fif",palette = palettedir3, size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_fif.cluster.0.5.png")
sc.pl.umap(adata, color="groups_six",size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_six.cluster.0.5.png")
sc.pl.umap(adata, color="groups_sev",size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_sev.cluster.0.5.png")


sc.pl.umap(adata, color="Tissue", size=0.3, sort_order= "False",save=".Tissue.cluster.png")
sc.pl.umap(adata, color="Phenotype", size=0.3,sort_order= "False",save=".Phenotype.cluster.png")

sc.pl.umap(adata, color="groups_pri",palette = palettedir2,size=0.3, sort_order= "False", save= ".groups_pri.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_sce",size=0.3, sort_order= "False", save= ".groups_sce.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_thi",palette = palettedir1,size=0.3, sort_order= "False", save= ".groups_thi.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_fif",palette = palettedir3, size=0.3, sort_order= "False", save= ".groups_fif.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_six",size=0.3, sort_order= "False", save= ".groups_six.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_sev",size=0.3, sort_order= "False", save= ".groups_sev.cluster.0.5.pdf")
5
sc.pl.umap(adata, color="groups_pri",palette = palettedir2,size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_pri.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_sce",size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_sce.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_thi",palette = palettedir1,size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_thi.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_fif",palette = palettedir3, size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_fif.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_six",size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_six.cluster.0.5.pdf")
sc.pl.umap(adata, color="groups_sev",size=0.3, legend_loc='on data',legend_fontweight="light",legend_fontsize=5,legend_fontoutline=0,sort_order= "False", save= ".on.groups_sev.cluster.0.5.pdf")


sc.pl.umap(adata, color="Tissue", size=0.3, sort_order= "False",save=".Tissue.cluster.pdf")
sc.pl.umap(adata, color="Phenotype", size=0.3,sort_order= "False",save=".Phenotype.cluster.pdf")

fig, (
[axLymphoid_CD4_Tn,axLymphoid_CD4_Tcm,axLymphoid_CD4_Tfh,axLymphoid_CD4_Th17],[axLymphoid_CD4_CTL,axLymphoid_CD4_Treg,axLymphoid_CD4_Tstr,axLymphoid_CD8_Tn],[axLymphoid_CD8_Tm,axLymphoid_CD8_Teff,axLymphoid_CD8_Tex,axLymphoid_CD8_Tstr],[axLymphoid_CD8_Tisg,axLymphoid_ILC,axLymphoid_MAIT,axLymphoid_NK_CD16__],[axLymphoid_NK_CD56__,axLymphoid_NKT,axLymphoid_B_naive,axLymphoid_B_memory],[axLymphoid_Plamsa_IGLC,axLymphoid_Plamsa_IGKC,axMyeloid_Mono,axMyeloid_cDC1],[axMyeloid_cDC2,axMyeloid_cDC3,axMyeloid_pDC,axMyeloid_Marco_C1QC],[axMyeloid_Marco_SPP1,axMyeloid_Marco_LYVE1,axMyeloid_Mast,axFibroblast_mCAF_POSTN],[axFibroblast_mCAF_KRT19,axFibroblast_mCAF_WNT5A,axFibroblast_iCAF_IGFBP6,axFibroblast_iCAF_IL6],[axFibroblast_iCAF_KCNN3,axFibroblast_apCAF_CD74,axMural_cell_Pericyte1,axMural_cell_Pericyte2],[axMural_cell_SMC1,axMural_cell_SMC2,axEC_Angiogenic,axEC_Artery],[axEC_Capillary,axEC_Vein,axEC_Lymph,axEC_Alveolar],[axEC_Glomerular,axEC_Sinusoidal,axEpithelium_Malig_Migration,axEpithelium_Malig_Cycle],[axEpithelium_Malig_cEMT,axEpithelium_Malig_Interferon,axEpithelium_Malig_Stress,axEpithelium_Malig_Basal],[axEpithelium_Malig_Glandular,axEpithelium_Malig_Multi,axEpithelium_Normal,xxx]
) = plt.subplots(15, 4, figsize=(16,45),
     )

sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD4_CTL'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD4_CTL)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD4_Tcm'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD4_Tcm)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD4_Tfh'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD4_Tfh)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD4_Th17'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD4_Th17)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD4_Tn'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD4_Tn)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD4_Treg'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD4_Treg)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD4_Tstr'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD4_Tstr)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD8_Teff'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD8_Teff)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD8_Tex'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD8_Tex)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD8_Tisg'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD8_Tisg)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD8_Tm'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD8_Tm)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD8_Tn'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD8_Tn)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_CD8_Tstr'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD8_Tstr)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_ILC'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_ILC)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_MAIT'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_MAIT)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_NK_CD16+'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_NK_CD16__)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_NK_CD56+'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_NK_CD56__)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_NKT'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_NKT)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_B_memory'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_B_memory)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_B_naive'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_B_naive)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_Plamsa_IGKC'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_Plamsa_IGKC)


sc.pl.umap(adata, color='groups_pri', groups=['Lymphoid_Plamsa_IGLC'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_Plamsa_IGLC)


sc.pl.umap(adata, color='groups_pri', groups=['Myeloid_Mono'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_Mono)


sc.pl.umap(adata, color='groups_pri', groups=['Myeloid_cDC1'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_cDC1)


sc.pl.umap(adata, color='groups_pri', groups=['Myeloid_cDC2'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_cDC2)


sc.pl.umap(adata, color='groups_pri', groups=['Myeloid_cDC3'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_cDC3)


sc.pl.umap(adata, color='groups_pri', groups=['Myeloid_pDC'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_pDC)


sc.pl.umap(adata, color='groups_pri', groups=['Myeloid_Marco_C1QC'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_Marco_C1QC)


sc.pl.umap(adata, color='groups_pri', groups=['Myeloid_Marco_SPP1'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_Marco_SPP1)


sc.pl.umap(adata, color='groups_pri', groups=['Myeloid_Marco_LYVE1'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_Marco_LYVE1)


sc.pl.umap(adata, color='groups_pri', groups=['Myeloid_Mast'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_Mast)


sc.pl.umap(adata, color='groups_pri', groups=['Fibroblast_apCAF_CD74'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_apCAF_CD74)


sc.pl.umap(adata, color='groups_pri', groups=['Fibroblast_iCAF_IGFBP6'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_iCAF_IGFBP6)


sc.pl.umap(adata, color='groups_pri', groups=['Fibroblast_iCAF_IL6'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_iCAF_IL6)


sc.pl.umap(adata, color='groups_pri', groups=['Fibroblast_iCAF_KCNN3'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_iCAF_KCNN3)


sc.pl.umap(adata, color='groups_pri', groups=['Fibroblast_mCAF_KRT19'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_mCAF_KRT19)


sc.pl.umap(adata, color='groups_pri', groups=['Fibroblast_mCAF_POSTN'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_mCAF_POSTN)


sc.pl.umap(adata, color='groups_pri', groups=['Fibroblast_mCAF_WNT5A'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_mCAF_WNT5A)


sc.pl.umap(adata, color='groups_pri', groups=['Mural_cell_Pericyte1'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMural_cell_Pericyte1)


sc.pl.umap(adata, color='groups_pri', groups=['Mural_cell_Pericyte2'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMural_cell_Pericyte2)


sc.pl.umap(adata, color='groups_pri', groups=['Mural_cell_SMC1'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMural_cell_SMC1)


sc.pl.umap(adata, color='groups_pri', groups=['Mural_cell_SMC2'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMural_cell_SMC2)


sc.pl.umap(adata, color='groups_pri', groups=['EC_Artery'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC_Artery)


sc.pl.umap(adata, color='groups_pri', groups=['EC_Capillary'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC_Capillary)


sc.pl.umap(adata, color='groups_pri', groups=['EC_Angiogenic'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC_Angiogenic)


sc.pl.umap(adata, color='groups_pri', groups=['EC_Vein'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC_Vein)


sc.pl.umap(adata, color='groups_pri', groups=['EC_Lymph'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC_Lymph)


sc.pl.umap(adata, color='groups_pri', groups=['EC_Alveolar'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC_Alveolar)


sc.pl.umap(adata, color='groups_pri', groups=['EC_Glomerular'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC_Glomerular)


sc.pl.umap(adata, color='groups_pri', groups=['EC_Sinusoidal'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC_Sinusoidal)


sc.pl.umap(adata, color='groups_pri', groups=['Epithelium_Malig_Basal'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig_Basal)


sc.pl.umap(adata, color='groups_pri', groups=['Epithelium_Malig_cEMT'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig_cEMT)


sc.pl.umap(adata, color='groups_pri', groups=['Epithelium_Malig_Cycle'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig_Cycle)


sc.pl.umap(adata, color='groups_pri', groups=['Epithelium_Malig_Glandular'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig_Glandular)


sc.pl.umap(adata, color='groups_pri', groups=['Epithelium_Malig_Interferon'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig_Interferon)


sc.pl.umap(adata, color='groups_pri', groups=['Epithelium_Malig_Migration'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig_Migration)


sc.pl.umap(adata, color='groups_pri', groups=['Epithelium_Malig_Multi'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig_Multi)


sc.pl.umap(adata, color='groups_pri', groups=['Epithelium_Malig_Stress'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig_Stress)


sc.pl.umap(adata, color='groups_pri', groups=['Epithelium_Normal'],  legend_loc='on data', palette = palettedir2,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Normal)




#plt.savefig("figures/sub.group_pri.celltype.0809.pdf")
plt.savefig("figures/sub.group_pri.celltype.0809.png")

fig, (
[axLymphoid_CD4,axLymphoid_CD8,axLymphoid_ILC,axLymphoid_MAIT],[axLymphoid_NK,axLymphoid_B,axLymphoid_Plamsa,axMyeloid],[axFibroblast,axMural_cell,axEC,axEpithelium_Malig],[axEpithelium_Normal,xxx,yyy,zzz]
) = plt.subplots(4, 4, figsize=(16,12),
     )

sc.pl.umap(adata, color='groups_sce', groups=['Lymphoid_CD4'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD4)


sc.pl.umap(adata, color='groups_sce', groups=['Lymphoid_CD8'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD8)


sc.pl.umap(adata, color='groups_sce', groups=['Lymphoid_ILC'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_ILC)


sc.pl.umap(adata, color='groups_sce', groups=['Lymphoid_MAIT'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_MAIT)


sc.pl.umap(adata, color='groups_sce', groups=['Lymphoid_NK'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_NK)


sc.pl.umap(adata, color='groups_sce', groups=['Lymphoid_B'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_B)


sc.pl.umap(adata, color='groups_sce', groups=['Lymphoid_Plamsa'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_Plamsa)


sc.pl.umap(adata, color='groups_sce', groups=['Myeloid'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid)


sc.pl.umap(adata, color='groups_sce', groups=['Fibroblast'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast)


sc.pl.umap(adata, color='groups_sce', groups=['Mural_cell'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMural_cell)


sc.pl.umap(adata, color='groups_sce', groups=['EC'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC)


sc.pl.umap(adata, color='groups_sce', groups=['Epithelium_Malig'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig)


sc.pl.umap(adata, color='groups_sce', groups=['Epithelium_Normal'],  legend_loc='on data', 
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Normal)


#plt.savefig("figures/sub.group_sce.celltype.0809.pdf")
plt.savefig("figures/sub.group_sce.celltype.0809.png")


fig, (
[axLymphoid_CD4,axLymphoid_CD8,axLymphoid_ILC,axLymphoid_MAIT],[axLymphoid_NK,axLymphoid_NKT,axLymphoid_B,axLymphoid_Plamsa],[axMyeloid_Mono,axMyeloid_DC,axMyeloid_Marco,axMyeloid_Mast],[axFibroblast_mCAF,axFibroblast_iCAF,axFibroblast_apCAF,axMural_Pericyte],[axMural_SMC,axEC,axEpithelium_Malig,axEpithelium_Normal]
) = plt.subplots(5, 4, figsize=(16,15),
     )

sc.pl.umap(adata, color='groups_thi', groups=['Lymphoid_CD4'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD4)


sc.pl.umap(adata, color='groups_thi', groups=['Lymphoid_CD8'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_CD8)


sc.pl.umap(adata, color='groups_thi', groups=['Lymphoid_ILC'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_ILC)


sc.pl.umap(adata, color='groups_thi', groups=['Lymphoid_MAIT'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_MAIT)


sc.pl.umap(adata, color='groups_thi', groups=['Lymphoid_NK'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_NK)


sc.pl.umap(adata, color='groups_thi', groups=['Lymphoid_NKT'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_NKT)


sc.pl.umap(adata, color='groups_thi', groups=['Lymphoid_B'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_B)


sc.pl.umap(adata, color='groups_thi', groups=['Lymphoid_Plamsa'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axLymphoid_Plamsa)


sc.pl.umap(adata, color='groups_thi', groups=['Myeloid_Mono'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_Mono)


sc.pl.umap(adata, color='groups_thi', groups=['Myeloid_DC'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_DC)


sc.pl.umap(adata, color='groups_thi', groups=['Myeloid_Marco'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_Marco)


sc.pl.umap(adata, color='groups_thi', groups=['Myeloid_Mast'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMyeloid_Mast)


sc.pl.umap(adata, color='groups_thi', groups=['Fibroblast_mCAF'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_mCAF)


sc.pl.umap(adata, color='groups_thi', groups=['Fibroblast_iCAF'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_iCAF)


sc.pl.umap(adata, color='groups_thi', groups=['Fibroblast_apCAF'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axFibroblast_apCAF)


sc.pl.umap(adata, color='groups_thi', groups=['Mural_Pericyte'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMural_Pericyte)


sc.pl.umap(adata, color='groups_thi', groups=['Mural_SMC'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axMural_SMC)


sc.pl.umap(adata, color='groups_thi', groups=['EC'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEC)


sc.pl.umap(adata, color='groups_thi', groups=['Epithelium_Malig'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Malig)


sc.pl.umap(adata, color='groups_thi', groups=['Epithelium_Normal'],  legend_loc='on data', palette = palettedir1,
           legend_fontweight="light",legend_fontsize=8,size=0.3, legend_fontoutline=0, ax=axEpithelium_Normal)

#plt.savefig("figures/sub.group_thi.celltype.0809.pdf")
plt.savefig("figures/sub.group_thi.celltype.0809.png")


