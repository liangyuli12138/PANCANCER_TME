import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/08.marker/filter_list")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872D2_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01872D2_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01872D2_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/SS200000929BL_D2_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("SS200000929BL_D2_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("SS200000929BL_D2_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="SS200000929BL_D2")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01615B5_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01615B5_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01615B5_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01615B5")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872D1_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01872D1_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01872D1_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01872D1")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D6_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01972D6_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01972D6_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01972D6")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D1_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01972D1_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01972D1_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01972D1")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C5_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01872C5_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01872C5_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01872C5")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C4_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01872C4_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01872C4_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01872C4")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D2_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01972D2_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01972D2_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01972D2")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972C1_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01972C1_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01972C1_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01972C1")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972B6_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01972B6_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01972B6_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01972B6")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01615B1_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01615B1_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01615B1_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01615B1")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872D3_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01872D3_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01872D3_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01872D3")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324E3_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B02324E3_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B02324E3_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B02324E3")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324B5_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B02324B5_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B02324B5_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B02324B5")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01613C6_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01613C6_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01613C6_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01613C6")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324E4_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B02324E4_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B02324E4_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B02324E4")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324A5_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B02324A5_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B02324A5_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B02324A5")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324A1_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B02324A1_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B02324A1_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B02324A1")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01613A5_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01613A5_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01613A5_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01613A5")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01615B2_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01615B2_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01615B2_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01615B2")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01613B1_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01613B1_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01613B1_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01613B1")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324E5_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B02324E5_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B02324E5_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B02324E5")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01613D1_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01613D1_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01613D1_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01613D1")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872D4_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01872D4_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01872D4_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01872D4")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972C4_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01972C4_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01972C4_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01972C4")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C6_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("D01872C6_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("D01872C6_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="D01872C6")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01613B2_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01613B2_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01613B2_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01613B2")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01317E6_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01317E6_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01317E6_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01317E6")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01615B3_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("B01615B3_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("B01615B3_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata, batch_categories="B01615B3")

sc.pp.normalize_total(adata_concat, target_sum=1e4)
sc.pp.log1p(adata_concat)
sc.pp.scale(adata_concat, max_value=10)

adata_concat.write_h5ad("pancancer.merge.all.st.h5ad")
adata_concat.obs.to_csv("pancancer.merge.all.st.obs")

marker_genes_dict = {
        'Lymphoid_T_NK':['CD3D','TRAC','GZMA','CCL5'],
        'Lymphoid_B':['MS4A1','BANK1','IGHM'],
        'Lymphoid_Plamsa':['MZB1','JCHAIN','IGLL5'],
        'Myeloid':['LYZ','C1QB','LST1'],
        'Mast':['TPSB2','MS4A2','CPA3'],
        'EC':['VWF','EMCN','PECAM1'],
        'Fibroblast':['LUM','FBLN1','DCN'],
        'Mural_cell':['RGS5','MYH11','MCAM'],
        'Epithelium': ['EPCAM', 'KRT8','ELF3'],
}

sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="celltype_sc", dot_max=0.9, categories_order=['Lymphoid_T_NK','Lymphoid_B','Lymphoid_Plamsa','Myeloid','Mast','EC','Fibroblast','Mural_cell','Epithelium'],vmax=3,
                  save=".cluster.markergene.pdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="celltype_sc", dot_max=0.9, categories_order=['Lymphoid_T_NK','Lymphoid_B','Lymphoid_Plamsa','Myeloid','Mast','EC','Fibroblast','Mural_cell','Epithelium'],use_raw=False, colorbar_title="mean z-score", vmin=-2, vmax=2, cmap="RdBu_r",
                  save=".cluster.markergene.zscore.pdf")

