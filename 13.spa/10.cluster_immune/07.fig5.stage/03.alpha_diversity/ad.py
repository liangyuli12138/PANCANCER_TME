import pandas as pd
import numpy as np

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/03.alpha_diversity")
otu_table = pd.read_csv('cell.merge.at', sep=',', header=0, index_col=0)

# 计算每个OTU的相对丰度的负对数
negative_log = -np.log(otu_table)

# 计算每个样本的Shannon多样性指数
shannon_index = negative_log * otu_table
shannon_index['Shannon'] = shannon_index.sum(axis=1)

shannon_index.to_csv('otu_table.csv')
