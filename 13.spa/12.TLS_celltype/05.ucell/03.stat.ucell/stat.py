import pandas as pd

# 读取文件路径
with open('at.list', 'r') as file:
    file_paths = [line.strip() for line in file]

# 读取csv文件并筛选
for file_path in file_paths:
    # 提取样品id
    sample_id = file_path.split('/')[-1].split('.')[0]
    # 读取csv文件
    df = pd.read_csv(file_path)
    # 筛选LM列为Lymphoid的细胞
    lymphoid_cells = df[df['LM'] == 'Lymphoid']
    # 保存筛选结果到文件
    output_file_path = sample_id + '_lymphoid.csv'
    lymphoid_cells.to_csv(output_file_path, index=False)


