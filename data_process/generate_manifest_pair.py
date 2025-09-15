import pandas as pd
import os


def generate_manifest(arc_file, output_file):

    current_path = os.getcwd()

    # Read arc.txt file
    with open(arc_file, 'r') as f:
        srr_ids = [line.strip() for line in f]

    # Create empty lists to store data
    sample_ids = []
    filepaths = []
    directions = []

    # Generate data
    for i, srr_id in enumerate(srr_ids, 1):
        # Add two rows for each sample
        sample_ids.extend([f'sample{i}'] * 2)

        # Generate file paths
        base_path = os.path.join(current_path, srr_id)
        filepaths.extend([f'{base_path}_1.fastq', f'{base_path}_2.fastq'])

        # Add direction
        directions.extend(['forward', 'reverse'])

    # Create DataFrame
    manifest_df = pd.DataFrame({
        'sample-id': sample_ids,
        'absolute-filepath': filepaths,
        'direction': directions
    })

    # Save as CSV file
    manifest_df.to_csv(output_file, index=False)
    print(f"Manifest file has been generated: {output_file}")
    print(f"Total samples: {len(srr_ids)}")
    print(f"Total rows: {len(manifest_df)}")


# Use the function
arc_file = 'arc.txt'
output_file = 'manifest.csv'
generate_manifest(arc_file, output_file)