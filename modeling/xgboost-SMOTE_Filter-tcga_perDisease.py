## Filtering: GT + ReadDepth + MappingQual + AllelicRatio + VCF hard filtering + exons + joint AF(1e-4) + no duplicates
## More features: conversion type: REF_ALT, AF, origin, count, VAF (neighbor, CN)
## saved combined datatable: combined_nodup.csv

## XGBoost + SMOTE ratio=1 (oversample minority to have balanced classes)

## Select samples with TCGA somatic:  0.1(14) ~ 0.9 (354)

## Train and test on each disease type


import os
import pandas as pd
import numpy as np
import joblib
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.metrics import accuracy_score, classification_report, f1_score, precision_recall_curve, auc, recall_score, roc_auc_score
from xgboost import XGBClassifier
from imblearn.over_sampling import SMOTENC
import matplotlib.pyplot as plt


model_name = "xgboost-SMOTE_Filter-tcga_perDisease"

# Set up paths
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
data_dir = os.environ.get('DATA_DIR', repo_root)
tcga_info_dir = os.environ.get('TCGA_INFO_DIR', os.path.join(data_dir, 'tcga_info'))
features_new_dir = os.environ.get('FEATURES_NEW_DIR', os.path.join(data_dir, 'features_new'))

### Read processed data

combined_df = pd.read_csv(os.path.join(features_new_dir, "combined_nodup.csv"), index_col=0)

disease_df = pd.read_csv(os.path.join(tcga_info_dir, "cases_info.txt"), sep='\t')
disease_df = disease_df[['submitter_id', 'disease_type']]
disease_df.rename(columns={'submitter_id': 'tcga_id'}, inplace=True)

combined_df['tcga_id'] = combined_df['sample_id'].replace(to_replace=r'-01[A-C]', value='', regex=True)
combined_df = pd.merge(combined_df, disease_df, on='tcga_id', how='left')

cases_info = pd.read_csv("/gpfs/commons/groups/gursoy_lab/xli/somatic_variant_prediction/tcga_info/cases_info_somatic.csv", index_col=0)
select_tcga = cases_info.loc[(cases_info['num_true_somatic'] > 14) & (cases_info['num_true_somatic'] < 354), 'submitter_id'].values
combined_df = combined_df[combined_df['tcga_id'].isin(select_tcga)]
combined_df = combined_df.drop(['tcga_id'], axis=1)

print('\n#Samples: ', len(combined_df['sample_id'].unique()))
print('Filtered Dataset: ', combined_df.shape)
print('#True somatic mutations: ', combined_df[['if_true_somatic']].values.sum())



results = []

for disease, disease_data in combined_df.groupby('disease_type'):
    sample_size = disease_data['sample_id'].nunique()
    if sample_size < 20:
        print(f"Skipping disease type: {disease} (sample size: {sample_size})")
        continue

    print(f"Processing disease type: {disease}")

    seed = 42
    unique_sample_ids = disease_data['sample_id'].unique()
    train_sample_ids, test_sample_ids = train_test_split(unique_sample_ids, test_size=0.2, random_state=seed)

    train_data = disease_data[disease_data['sample_id'].isin(train_sample_ids)]
    test_data = disease_data[disease_data['sample_id'].isin(test_sample_ids)]

    X_train = train_data.drop(['if_true_somatic', 'sample_id', 'disease_type'], axis=1)
    y_train = train_data['if_true_somatic']
    X_test = test_data.drop(['if_true_somatic', 'sample_id', 'disease_type'], axis=1)
    y_test = test_data['if_true_somatic']

    categorical_cols = X_train.select_dtypes(include=['object']).columns.tolist()
    numerical_cols = X_train.select_dtypes(exclude=['object']).columns.tolist()

    numerical_transformer = StandardScaler()
    categorical_transformer = OneHotEncoder(handle_unknown='ignore')

    preprocessor = ColumnTransformer(
        transformers=[
            ('num', numerical_transformer, numerical_cols),
            ('cat', categorical_transformer, categorical_cols)
        ]
    )

    # Resample using SMOTE
    smote = SMOTENC(random_state=seed, categorical_features=categorical_cols, sampling_strategy=1.0)
    X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)

    # Build the pipeline
    xgb_classifier = XGBClassifier(random_state=seed)

    pipeline = Pipeline(steps=[('preprocessor', preprocessor),
                               ('classifier', xgb_classifier)])

    # Train the model
    pipeline.fit(X_train_resampled, y_train_resampled)

    # Evaluate the performance
    y_pred = pipeline.predict(X_test)
    y_prob = pipeline.predict_proba(X_test)[:, 1]

    f1 = f1_score(y_test, y_pred)
    f1_weighted = f1_score(y_test, y_pred, average='weighted')
    recall_class_1 = recall_score(y_test, y_pred, pos_label=1)
    auroc_class_1 = roc_auc_score(y_test, y_prob)

    precision, recall, _ = precision_recall_curve(y_test, y_prob)
    auprc = auc(recall, precision)

    num_class_1 = disease_data['if_true_somatic'].sum()
    total_cases = disease_data.shape[0]

    results.append({
        'disease_type': disease,
        'sample_size': sample_size,
        'total_cases': total_cases,
        'num_class_1': num_class_1,
        'f1_score': f1,
        'weighted_f1': f1_weighted,
        'recall_class_1': recall_class_1,
        'auroc': auroc_class_1,
        'auprc': auprc
    })


results_df = pd.DataFrame(results)
results_df.to_csv(f'{model_name}_metrics-by-disease.csv', index=False)


