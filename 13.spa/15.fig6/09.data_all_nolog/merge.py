import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/09.data_all_nolog")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/SS200000929BL_D2_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/SS200000929BL_D2.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "SS200000929BL_D2_" + adata.obs.index
adata_concat = adata
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01615B5_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B01615B5.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B01615B5_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D6_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01972D6.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01972D6_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D1_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01972D1.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01972D1_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C5_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01872C5.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01872C5_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C4_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01872C4.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01872C4_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D2_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01972D2.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01972D2_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972C1_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01972C1.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01972C1_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972B6_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01972B6.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01972B6_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01615B1_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B01615B1.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B01615B1_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872D3_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01872D3.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01872D3_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324E3_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B02324E3.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B02324E3_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324B5_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B02324B5.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B02324B5_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01613C6_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B01613C6.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B01613C6_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324E4_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B02324E4.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B02324E4_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324A5_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B02324A5.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B02324A5_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01615B2_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B01615B2.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B01615B2_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01613B1_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B01613B1.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B01613B1_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01613D1_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B01613D1.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B01613D1_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872D4_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01872D4.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01872D4_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972C4_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01972C4.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01972C4_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C6_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/D01872C6.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "D01872C6_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01613B2_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B01613B2.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B01613B2_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01317E6_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B01317E6.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B01317E6_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B01615B3_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/B01615B3.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "B01615B3_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

import re
# 获取原始索引
original_index = adata_concat.obs_names.to_list()

# 使用正则表达式替换后缀
pattern = r"-\d+(-\d+)*$"
new_index = [re.sub(pattern, '', idx) for idx in original_index]

# 更新索引
adata_concat.obs_names = new_index

sc.pp.normalize_total(adata_concat,target_sum=10000)

adata_concat.write_h5ad("pancancer.icar.background.nolog.cell.h5ad")
adata_concat.obs.to_csv("pancancer.icar.background.nolog.cell.obs")

