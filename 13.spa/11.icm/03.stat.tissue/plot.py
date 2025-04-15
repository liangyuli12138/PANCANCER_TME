import pandas as pd
import matplotlib.pyplot as plt

# 读取文件
data = pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/03.stat.tissue/input.icm.tissue.csv', delimiter='\t')

# 按照tumore进行分组，并将icm1、icm2和icm3加和
grouped_data = data.groupby('tumore').sum()

# 提取所需数据
tumors = grouped_data.index
cm1 = grouped_data['ICM1'].tolist()
cm2 = grouped_data['ICM2'].tolist()
cm3 = grouped_data['ICM3'].tolist()

colors = ['#37B01C', '#F66B6B', '#BC68F4']

# 绘制堆叠柱状图并设置颜色
plt.bar(tumors, cm1, label='ICM1', color=colors[0])
plt.bar(tumors, cm2, bottom=cm1, label='ICM2', color=colors[1])
plt.bar(tumors, cm3, bottom=[i+j for i,j in zip(cm1, cm2)], label='ICM3', color=colors[2])


# 添加图例和标签
plt.legend(facecolor='white')  
plt.xlabel('Tumor')
plt.ylabel('Counts')
#plt.title('Stacked Bar Chart')

# 显示横坐标标签
plt.xticks(tumors, rotation=45)

# 显示纵坐标标签
plt.yticks(range(0, max(cm1+cm2+cm3)+1, 20000))

plt.tight_layout()

plt.savefig('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/03.stat.tissue/output.icm.tissue.png', format='png')

# 显示图形
plt.close()

