from pathlib import Path
import os
import anndata
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import infercnvpy
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancancer/08.filter_by_gene/09.cancer_cnv/05.find_cancer_donor/aaaa/bbbb")
adata_ori = sc.read_h5ad("plot_clustering.aaaa.bbbb.remarker.h5ad")
adata_raw = adata_ori.raw.to_adata()

#atlist = pd.read_csv("pancancer.final.0317.obs.csv.at",index_col=0)
#adata_raw = adata_raw[cellist["cell"],:]
#adata_raw.obs = adata_raw.obs.join(atlist)
#adata_raw = adata_ori.raw.to_adata()
#adata_raw = adata_ori

infercnvpy.io.genomic_position_from_gtf("genes.ref.gtf", adata=adata_raw)

infercnvpy.tl.infercnv(
    adata_raw,
    reference_key="groups_tl",
    reference_cat=[
        "T_NK",
        "B",
        "Plasmacyte",
        "Myeloid",
        "mastocyte",
        "Fibroblast",
        "SMC",
        "EC",
     ],
    window_size=250,
#    n_jobs=20
)

infercnvpy.pl.chromosome_heatmap(adata_raw, groupby="groups_tl", save=".aaaa.bbbb.heatmap_cell_type.png")
infercnvpy.pl.chromosome_heatmap(adata_raw, groupby="groups_tl", save=".aaaa.bbbb.heatmap_cell_type.pdf")

infercnvpy.tl.pca(adata_raw)
infercnvpy.pp.neighbors(adata_raw)
infercnvpy.tl.leiden(adata_raw)

infercnvpy.tl.umap(adata_raw)
infercnvpy.tl.cnv_score(adata_raw)

adata_raw.write_h5ad("plot_clustering.remarker.aaaa.bbbb.cnv.cluster.score.h5ad")
adata_raw.obs.to_csv("plot_clustering.remarker.aaaa.bbbb.cnv.cluster.score.obs.csv")
adata_raw.var.to_csv("plot_clustering.remarker.aaaa.bbbb.cnv.cluster.score.var.csv")

#_, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(2, 3, figsize=(24, 16))
infercnvpy.pl.umap(adata_raw, color="cnv_leiden",legend_loc='on data', legend_fontweight='normal', show=False, save=".aaaa.bbbb.cnv.umap.cnv_leiden.png")
infercnvpy.pl.umap(adata_raw, color="cnv_score", legend_fontweight='normal', show=False, save=".aaaa.bbbb.cnv.umap.cnv_score.png")
infercnvpy.pl.umap(adata_raw, color="groups_tl", legend_fontweight='normal', show=False, save=".aaaa.bbbb.cnv.umap.Cell_type.png")
infercnvpy.pl.umap(adata_raw, color="cnv_leiden",legend_loc='on data', legend_fontweight='normal', show=False, save=".aaaa.bbbb.cnv.umap.cnv_leiden.pdf")
infercnvpy.pl.umap(adata_raw, color="cnv_score", legend_fontweight='normal', show=False, save=".aaaa.bbbb.cnv.umap.cnv_score.pdf")
infercnvpy.pl.umap(adata_raw, color="groups_tl", legend_fontweight='normal', show=False, save=".aaaa.bbbb.cnv.umap.Cell_type.pdf")

sc.pl.umap(adata_raw, color="cnv_leiden",  legend_fontweight='normal', show=False, save=".aaaa.bbbb.sc.umap.cnv_leiden.png")
sc.pl.umap(adata_raw, color="cnv_score", legend_fontweight='normal', show=False, save=".aaaa.bbbb.sc.umap.cnv_score.png")
sc.pl.umap(adata_raw, color="groups_tl", legend_fontweight='normal', show=False, save=".aaaa.bbbb.sc.umap.Cell_type.png")
sc.pl.umap(adata_raw, color="cnv_leiden",  legend_fontweight='normal', show=False, save=".aaaa.bbbb.sc.umap.cnv_leiden.pdf")
sc.pl.umap(adata_raw, color="cnv_score", legend_fontweight='normal', show=False, save=".aaaa.bbbb.sc.umap.cnv_score.pdf")
sc.pl.umap(adata_raw, color="groups_tl", legend_fontweight='normal', show=False, save=".aaaa.bbbb.sc.umap.Cell_type.pdf")

adata_raw = adata_raw[(adata_raw.obs["groups_tl"] == "Epithelium"),:]

infercnvpy.pl.umap(adata_raw, color="cnv_leiden",legend_loc='on data', legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.cnv.umap.cnv_leiden.png")
infercnvpy.pl.umap(adata_raw, color="cnv_score", legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.cnv.umap.cnv_score.png")
infercnvpy.pl.umap(adata_raw, color="groups_tl", legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.cnv.umap.Cell_type.png")
infercnvpy.pl.umap(adata_raw, color="cnv_leiden",legend_loc='on data', legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.cnv.umap.cnv_leiden.pdf")
infercnvpy.pl.umap(adata_raw, color="cnv_score", legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.cnv.umap.cnv_score.pdf")
infercnvpy.pl.umap(adata_raw, color="groups_tl", legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.cnv.umap.Cell_type.pdf")

sc.pl.umap(adata_raw, color="cnv_leiden",  legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.sc.umap.cnv_leiden.png")
sc.pl.umap(adata_raw, color="cnv_score", legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.sc.umap.cnv_score.png")
sc.pl.umap(adata_raw, color="groups_tl", legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.sc.umap.Cell_type.png")
sc.pl.umap(adata_raw, color="cnv_leiden",  legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.sc.umap.cnv_leiden.pdf")
sc.pl.umap(adata_raw, color="cnv_score", legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.sc.umap.cnv_score.pdf")
sc.pl.umap(adata_raw, color="groups_tl", legend_fontweight='normal', show=False, save=".ep.aaaa.bbbb.sc.umap.Cell_type.pdf")
