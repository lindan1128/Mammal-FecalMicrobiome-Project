import pandas as pd
import os


def generate_manifest(arc_file, output_file):
    
    current_path = os.getcwd()

    # 读取arc.txt文件
    with open(arc_file, 'r') as f:
        srr_ids = [line.strip() for line in f]

    # 创建空列表存储数据
    sample_ids = []
    filepaths = []
    directions = []  # 添加direction列表

    # 生成数据
    for i, srr_id in enumerate(srr_ids, 1):
        # 每个样本只添加一行
        sample_ids.append(f'sample{i}')

        # 生成文件路径，只使用_1.fastq
        base_path = os.path.join(current_path, srr_id)
        filepaths.append(f'{base_path}_1.fastq')

        # 添加direction，都是forward
        directions.append('forward')

    # 创建DataFrame
    manifest_df = pd.DataFrame({
        'sample-id': sample_ids,
        'absolute-filepath': filepaths,
        'direction': directions
    })

    # 保存为CSV文件
    manifest_df.to_csv(output_file, index=False)
    print(f"Manifest file has been generated: {output_file}")
    print(f"Total samples: {len(srr_ids)}")
    print(f"Total rows: {len(manifest_df)}")


# 使用函数
arc_file = 'arc.txt'
output_file = 'manifest.csv'
generate_manifest(arc_file, output_file)