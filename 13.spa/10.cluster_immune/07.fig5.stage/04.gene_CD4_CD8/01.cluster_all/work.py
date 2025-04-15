import s.canpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
import numpy as np

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/04.gene_CD4_CD8/01.cluster_all")
adata = sc.read('cluster_expression.ori.h5ad')
csv_data = pd.read_csv('cluster2group.csv')
adata.obs['new_groups'] = csv_data.set_index('cell')['new_groups']
adata.obs['stage'] = np.where(adata.obs['new_groups'].isin(['Lymphoid1', 'Lymphoid2']), 'Early',
                             np.where(adata.obs['new_groups'] == 'Lymphoid3', 'Late', 'other'))

# 提取stage为Early和Late的两大类
early_adata = adata[adata.obs['stage'] == 'Early']
late_adata = adata[adata.obs['stage'] == 'Late']

# 对每个cell的基因表达进行加和处理
early_sum = np.sum(early_adata.X, axis=0)
late_sum = np.sum(late_adata.X, axis=0)

# 对stage的基因表达再做normalize处理
early_sum = np.array(early_sum / np.sum(early_sum))
late_sum = np.array(late_sum / np.sum(late_sum))

# 计算差异分析
fc = np.log2(late_sum / early_sum)
_, pvals = stats.ttest_ind(early_adata.X, late_adata.X)

# 计算FDR的-log10值
pvals = np.nan_to_num(pvals)
fdr = -np.log10(multipletests(pvals, method='fdr_bh')[1])

# 创建结果DataFrame
diff_df = pd.DataFrame({'fc': fc, 'fdr': fdr}, index=adata.var_names)

# 输出结果为CSV文件
diff_df.to_csv('diff.out.csv')


