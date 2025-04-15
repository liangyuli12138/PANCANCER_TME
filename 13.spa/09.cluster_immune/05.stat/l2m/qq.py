import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

data = pd.read_csv('l2m.stat.merge.csv', index_col=0)

x = data['Lymphoid_percentage']
y = data['area']

plt.figure(figsize=(8, 6))
#stats.probplot(x, dist="norm", plot=plt)
stats.probplot(x, plot=plt)
plt.title('QQ_Lymphoid_percentage_vs_area')
plt.xlabel('Lymphoid_percentage')
plt.ylabel('area')
plt.savefig('QQ_Lymphoid_percentage_vs_area.png', dpi=300)

plt.figure(figsize=(8, 6))
sns.pplot(data=data,x='Lymphoid_percentage',y='area',kind='qq', height=4, aspect=2, display_kws={"identity":False, "fit":True, "reg":True, "ci":0.025})
plt.title('QQ_Lymphoid_percentage_vs_area')
plt.xlabel('Lymphoid_percentage')
plt.ylabel('area')
plt.savefig('QQ_Lymphoid_percentage_vs_area.png', dpi=300)


