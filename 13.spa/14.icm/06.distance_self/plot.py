import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/06.distance_self")

data = pd.read_csv("Lymphoid123.distance.stat.csv")

matrix = pd.pivot_table(data, index='lym', columns='type1', values='dis', aggfunc='mean')
genes = ["Bn","Bm","Plamsa","Teff","Tex","Tfh","Tm","Tn","Treg","Marco","DC","Fibroblast","EC","Epithelium","Stroma"]
colors = sns.color_palette("viridis", len(genes))
mean_plot = matrix.plot(kind="line", y=genes, marker="o")
mean_plot.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.savefig('test.png', dpi=300)

