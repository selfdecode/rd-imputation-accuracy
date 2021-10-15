import pandas as pd
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(usage='%(prog)s')
parser.add_argument('--s', dest='sample', type=str, required=True, default="", nargs=1, help='')

try:
    args = parser.parse_args()
except:
    sys.exit(0)
    
sample_file = args.sample[0]

df_sample = pd.read_table(sample_file)
s_gp = df_sample.groupby("imputed_ids")
for i,x in s_gp:
    output_name = f"{i}_ImputationAccuracy.txt"
    imputed_sample_ids=x.imputed_ids.unique()[0]
    WGS_sample_ids=x.WGS_ids.unique()[0]
    sTP=np.sum(x.TP)
    sTN=np.sum(x.TN)
    sFP=np.sum(x.FP)
    sFN=np.sum(x.FN)
    sTPr=np.round(sTP/(sTP+sFP), decimals=3)
    sFPr=np.round(sFP/(sTP+sFP), decimals=3)
    sFNr=np.round(sFN/(sTN+sFN), decimals=3)
    sTNr=np.round(sTN/(sTN+sFN), decimals=3)
    precision = np.round(np.divide(sTP, np.add(sTP, sFP)), decimals=3)
    recall = np.round(np.divide(sTP, np.add(sTP, sFN)), decimals=3)
    top = np.multiply(precision, recall)
    bottom = np.add(precision, recall)
    sRMSE=np.round(np.mean(x.RMSE), decimals=3)
    P0_per_sample=np.round(np.mean(x.concordance_P0), decimals=3)
    r2_per_sample=np.round(np.mean(x.r2), decimals=3)
    per_sample_f1 = np.round(np.multiply(2.0, np.divide(top,bottom)), decimals=3)
    labels_s=['imputed_ids','WGS_ids','F-score','concordance_P0','r2', 'precision', 'recall', 'TP', 'TN', 'FP', 'FN', 'TP_ratio', 'TN_ratio', 'FP_ratio', 'FN_ratio', 'RMSE']
    merged_results_per_sample=np.column_stack((imputed_sample_ids,WGS_sample_ids,per_sample_f1,P0_per_sample,r2_per_sample, precision, recall, sTP, sTN, sFP, sFN, sTPr, sTNr, sFPr, sFNr, sRMSE))
    merged_results_per_sample=np.vstack((np.asarray(labels_s),merged_results_per_sample))
    np.savetxt(output_name, merged_results_per_sample, fmt="%s", delimiter="\t")
print("Process Completed.")
