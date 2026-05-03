#!/usr/bin/env python3
import os
import pandas as pd

repo_root = os.getcwd()
features_dir = os.getenv('FEATURES_DIR', os.path.join(repo_root, 'features_new'))
cases_info_path = os.getenv('CASES_INFO', os.path.join(repo_root, 'tcga_info', 'cases_info_somatic.csv'))

if os.path.exists(cases_info_path):
    cases_info = pd.read_csv(cases_info_path)
else:
    cases_info = pd.DataFrame(columns=['submitter_id', 'primary_site'])

files = [f for f in os.listdir(features_dir) if f.endswith('_features_filtered_GT_RD_MQ_AF_vcf_exons_jointAF_vafs.csv')]
frames = []
for filename in files:
    df = pd.read_csv(os.path.join(features_dir, filename), index_col=0)
    df['count'] = len(df)
    sample_id = filename.split('_')[0]
    tcga_id = sample_id.replace('-01A', '').replace('-01B', '').replace('-01C', '')
    df['REF_ALT'] = df['REF'] + '_' + df['ALT']
    if tcga_id in cases_info['submitter_id'].values:
        df['primary_site'] = cases_info.loc[cases_info['submitter_id'] == tcga_id, 'primary_site'].iloc[0]
    else:
        df['primary_site'] = 'Unknown'
    df['sample_id'] = sample_id
    df_combined = df.drop(columns=[col for col in ['REF', 'ALT', 'gt'] if col in df.columns])
    frames.append(df_combined)

combined_df = pd.concat(frames, axis=0)
combined_df = combined_df[~combined_df.index.duplicated(keep=False)]
combined_df.to_csv(os.path.join(features_dir, 'combined_nodup.csv'))
