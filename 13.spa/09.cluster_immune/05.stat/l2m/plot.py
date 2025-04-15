import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
data = pd.read_csv('l2m.stat.merge.csv', index_col=0)

fig, ax = plt.subplots(figsize=(8,6))
ax = sns.violinplot(data, x='group', y='Lymphoid_percentage', palette='tab10', scale='width',linewidth=0.8, height=1,order=['group4','group5','group0','group1','group2','group3','group6'])
sns.stripplot(data, x='group', y='Lymphoid_percentage', jitter=0.2, size=2, alpha=0.8, color='black', order=['group4','group5','group0','group1','group2','group3','group6'])
plt.savefig('output.Lymphoid_percentage.png', dpi=300) 

fig, ax = plt.subplots(figsize=(8,6))
ax = sns.violinplot(data, x='group', y='Myeloid_percentage', palette='tab10', scale='width',linewidth=0.8, height=1,order=['group4','group5','group0','group1','group2','group3','group6'])
sns.stripplot(data, x='group', y='Myeloid_percentage', jitter=0.2, size=2, alpha=0.8, color='black', order=['group4','group5','group0','group1','group2','group3','group6'])
plt.savefig('output.Myeloid_percentage.png', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))
ax = sns.violinplot(data, x='group', y='log10(Lym/Mye_ratio)', palette='tab10', scale='width',linewidth=0.8, height=1,order=['group4','group5','group0','group1','group2','group3','group6'])
sns.stripplot(data, x='group', y='log10(Lym/Mye_ratio)',jitter=0.2, size=2, alpha=0.8, color='black', order=['group4','group5','group0','group1','group2','group3','group6'])
plt.savefig('output.log10_LymMye_ratio.png', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))
ax = sns.violinplot(data, x='group', y='density', palette='tab10', scale='width',linewidth=0.8, height=1,order=['group4','group5','group0','group1','group2','group3','group6'])
sns.stripplot(data, x='group', y='density',jitter=0.2, size=2, alpha=0.8, color='black', order=['group4','group5','group0','group1','group2','group3','group6'])
plt.savefig('output.density.png', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))
ax = sns.violinplot(data, x='group', y='area', palette='tab10', scale='width',linewidth=0.8, height=1,order=['group4','group5','group0','group1','group2','group3','group6'])
ax.set_ylim(0, 1)
sns.stripplot(data, x='group', y='area',jitter=0.2, size=2, alpha=0.8, color='black', order=['group4','group5','group0','group1','group2','group3','group6'])
plt.savefig('output.area.png', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))
ax = sns.violinplot(data, x='group', y='elongation', palette='tab10', scale='width',linewidth=0.8, height=1,order=['group4','group5','group0','group1','group2','group3','group6'])
sns.stripplot(data, x='group', y='elongation',jitter=0.2, size=2, alpha=0.8, color='black', order=['group4','group5','group0','group1','group2','group3','group6'])
plt.savefig('output.elongation.png', dpi=300)

