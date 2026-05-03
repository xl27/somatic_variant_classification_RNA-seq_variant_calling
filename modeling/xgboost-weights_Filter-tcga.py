## Filtering: GT + ReadDepth + MappingQual + AllelicRatio + VCF hard filtering + exons + joint AF(1e-4) + no duplicates
## More features: conversion type: REF_ALT, AF, origin, count, VAF (neighbor, CN)
## saved combined datatable: combined_nodup.csv

## XGBoost: Optimal Class Weights = 4
## Test on each disease type
## Select samples with TCGA somatic:  0.1(14) ~ 0.9 (354)




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
import matplotlib.pyplot as plt


model_name = "xgboost-weights_Filter-tcga"

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


### Split

seed = 42
unique_sample_ids = combined_df['sample_id'].unique()
train_sample_ids, test_sample_ids = train_test_split(unique_sample_ids, test_size=0.2, random_state=seed)

train_data = combined_df[combined_df['sample_id'].isin(train_sample_ids)]
test_data = combined_df[combined_df['sample_id'].isin(test_sample_ids)]


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


### Build the model

class_weights = 4
print('\nscale_pos_weight=', class_weights)


xgb_classifier = XGBClassifier(scale_pos_weight=class_weights, random_state=seed)

pipeline = Pipeline(steps=[('preprocessor', preprocessor),
                           ('classifier', xgb_classifier)])


pipeline.fit(X_train, y_train)


# Save the model
joblib.dump(pipeline, f'{model_name}_model.pkl')


### Evaluate the performance 

y_pred = pipeline.predict(X_test)
print('\n#predicted somatic: ', y_pred.sum())

f1 = f1_score(y_test, y_pred)
f1_weighted = f1_score(y_test, y_pred, average='weighted')
recall_class_1 = recall_score(y_test, y_pred, pos_label=1)
classification_rep = classification_report(y_test, y_pred)

print(f'\nF1 Score (class 1): {f1}')
print(f'F1 Score (weighted average): {f1_weighted}')
print(f'Accuracy (class 1 recall): {recall_class_1}')

print('\nClassification Report:\n', classification_rep)


# Compute precision-recall 
y_prob = pipeline.predict_proba(X_test)[:, 1]

precision, recall, thresholds = precision_recall_curve(y_test, y_prob)

auprc = auc(recall, precision)
print(f"AUPRC: {auprc:.4f}")

auroc_class_1 = roc_auc_score(y_test, y_prob)
print(f"AUROC for class 1: {auroc_class_1:.2f}")


# Extract feature importances
feature_importances = pipeline.named_steps['classifier'].feature_importances_
feature_names = pipeline.named_steps['preprocessor'].get_feature_names_out()
feature_importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': feature_importances})
feature_importance_df = feature_importance_df.sort_values(by='Importance', ascending=False)

pd.set_option('display.max_rows', 200)
print("Feature Importances:\n", feature_importance_df)



### Evaluate the model for each disease type in the test set

grouped = test_data.groupby('disease_type')

results = []

for disease, group_data in grouped:
    print(f"\nProcessing disease type: {disease}")

    X_group_test = group_data.drop(['if_true_somatic', 'sample_id', 'disease_type'], axis=1)
    y_group_test = group_data['if_true_somatic']

    y_pred = pipeline.predict(X_group_test)
    y_prob = pipeline.predict_proba(X_group_test)[:, 1]

    f1 = f1_score(y_group_test, y_pred)
    f1_weighted = f1_score(y_group_test, y_pred, average='weighted')
    recall_class_1 = recall_score(y_group_test, y_pred, pos_label=1)
    auroc_class_1 = roc_auc_score(y_group_test, y_prob)

    precision, recall, _ = precision_recall_curve(y_group_test, y_prob)
    auprc = auc(recall, precision)

    sample_size = group_data['sample_id'].nunique()
    num_class_1 = y_group_test.sum()
    total_cases = group_data.shape[0]

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



