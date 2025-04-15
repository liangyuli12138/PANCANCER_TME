import pandas as pd
import plotly.graph_objects as go
import plotly.offline as pyo
import os

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/06.stat/sankey")

# 读取csv文件并提取数据
data = pd.read_csv("tissue.for.sankey.csv", header=None, names=["source", "target", "value"])
color_data = pd.read_csv("tissue.color.csv", header=None, names=["node", "color"])

# 创建节点和链接
nodes = list(set(data["source"]).union(set(data["target"])))
nodes = sorted(nodes)
node_dict = {node: i for i, node in enumerate(nodes)}
source_indices = [node_dict[source] for source in data["source"]]
target_indices = [node_dict[target] for target in data["target"]]
values = data["value"]

# 设置节点颜色和链接颜色
node_colors = []
link_colors = []
for node in nodes:
    color = color_data[color_data["node"] == node]["color"].values[0]
    node_colors.append(color)
    
for target in data["target"]:
    color = color_data[color_data["node"] == target]["color"].values[0]
    link_colors.append(color)

sorted_labels = sorted(nodes)

# 绘制桑基图
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=600,
        thickness=40,
        label=sorted_labels,
        color=node_colors
    ),
    link=dict(
        source=source_indices,
        target=target_indices,
        value=values,
        color=link_colors,
    )
)])

fig.update_layout(
    height=1000,
    width=800
)

# 保存输出文件为test.html
pyo.plot(fig, filename="icar.all.index.sankey.html")

fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=600,
        thickness=40,
        label=[""]*len(nodes),
        color=node_colors
    ),
    link=dict(
        source=source_indices,
        target=target_indices,
        value=values,
        color=link_colors,
    )
)])

fig.update_layout(
    height=1000,
    width=800
)

# 保存输出文件为test.html
pyo.plot(fig, filename="icar.all.sankey.html")

