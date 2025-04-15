import pandas as pd
import numpy as np
from scipy.stats import zscore

data = pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden2/cluster.immune.cell.csv', index_col=0)
columns_to_zscore = data.columns[0:]
zscored_data = data[columns_to_zscore].apply(zscore)
zscored_data = zscored_data.clip(-3, 3)
zscored_data.to_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden3/zscored_data.csv')

