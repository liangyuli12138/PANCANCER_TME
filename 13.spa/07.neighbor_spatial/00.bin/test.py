from spatial_neighbors import dot_plot

pd = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/07.neighbor_spatial/00.bin/out.nbor_pvals.csv",index_col=0)
df_to_merge = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/07.neighbor_spatial/00.bin/out.zscore.csv",index_col=0)

df = pd.melt(df, id_vars ='celltype', value_name='FDR')
df_to_merge = pd.melt(df_to_merge, id_vars ='celltype', value_name='z_score')
data = pd.merge(df, df_to_merge)

data['-logFDR'] = -np.log10(data['FDR'])
data = data[(data.celltype != 'UN') & (data.variable != 'UN')]
mask = (data['FDR'] > 0.05)
data.loc[mask, ['-logFDR', 'z_score']] = np.nan
#data.to_csv(f'{group}_{region}_subtype_table.csv')
data.to_csv('subtype_table.csv')
p = dot_plot(data, 'celltype', 'variable' ,'-logFDR', 'z_score',
colorbar_title='z score', size_title='${FDR}$',
fontstyle='normal', cmap='viridis', figsize=(5,5), rotation=90
)

#ggsave(p, f'{group}_{region}_subtype.pdf')
ggsave(p, 'subtype.pdf')

