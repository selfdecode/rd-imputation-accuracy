import multiprocessing as mp
import timeit
import pandas as pd
import glob
import os
from functools import partial 
import subprocess as sp 
import numpy as np
from cyvcf2 import VCF
#from scipy.stats import linregress #r2
import gzip
import argparse
import sys

max_total_rows = 100000 #maximun number of rows to be loading together, summing all chunks loaded by all cores

n_max=100000 #maximun number of variants per chunk for each core, avoid RAM overload
n_min=10000 #minimun number of variants per chunk for each core, avoid interprocess communication overload

GA_file=""
WGS_file=""
IMPUTED_file=""
REF_file=""
IMPUTED_file_bcf=""
WGS_file_bcf=""

X_mode="False"
DEBUG=False

vout=""
sout=""

disable_DS=False

#just for autoencoder benchmarking purposes
multimask_mode=False

if(DEBUG==True):
    import csv

#new chunk, works on ND arrays instead of taking indexes
def chunk(data, ncores=mp.cpu_count()):

    n=int(round(len(data)/ncores-0.5))
    if(n>n_max):
        n=n_max
    if(n<n_min):
        n=n_min
    chunks=[data[i:i + n] for i in range(0, len(data), n)]

    return chunks

def convert_and_reshape(in_dict, x_in, y_in):

    l=len(in_dict[0])

    if(DEBUG==True):
        print("INPUTS",x_in,y_in)


    x=np.round(np.copy(x_in))
    y=np.round(np.copy(y_in))

    if(DEBUG==True):
        print("ROUND",x,y)


    x_o=[[in_dict[i] for i in row] for row in x]
    y_o=[[in_dict[i] for i in row] for row in y]

    x_o=np.asarray(x_o)
    y_o=np.asarray(y_o)

    if(DEBUG==True):
        print("CONVERTED", x_o,y_o)

    x_o=x_o.reshape(x_o.shape[0],x_o.shape[1]*l)
    y_o=y_o.reshape(y_o.shape[0],y_o.shape[1]*l)

    return x_o, y_o

def rmse(x, y):

    results=dict()

    #squared error
    se=np.subtract(x,y)**2

    #root mean squared error per variant
    var_rmse=np.sqrt(np.mean(se, axis=1))
    
    #squared error summmed per sample, average after merging chunks
    sample_ses=np.sum(se, axis=0)

    results['var_rmse'] = var_rmse
    results['sample_ses'] = sample_ses
    results['N'] = len(y)

    return results

def f1_score(x_in, y_in):

    #avoid division by zero
    eps = 1e-8

    #go along all axes since this function takes one variant per run, or one sample per run
    #axis: 0 (vertical, per sample), 1 (horizontal, per SNP), None (all, per sample+variant)

    dose_to_binary = {0: [1,0], 1: [1,1], 2: [0,1], -1:[0,0]}
    #dose_to_prob = {0: [1,0,0], 1: [0,1,0], 2: [0,0,1]}

    x, y = convert_and_reshape(dose_to_binary, x_in, y_in)
    #x, y = convert_and_reshape(dose_to_prob, x_in, y_in)

    results=dict()

    for axis in [0,1]:

        #TP = np.count_nonzero(np.multiply(x, y), axis=axis)
        #FP = np.count_nonzero(np.multiply(x, np.subtract(y,1.0)), axis=axis)
        #FN = np.count_nonzero(np.multiply(np.subtract(x,1.0),y), axis=axis)

        # True Positive (TP): we predict a label of 1 (positive), and the true label is 1.       
        TP = np.sum(np.logical_and(x == 1, y == 1), axis=axis)
 
        # True Negative (TN): we predict a label of 0 (negative), and the true label is 0.
        TN = np.sum(np.logical_and(x == 0, y == 0), axis=axis)
 
        # False Positive (FP): we predict a label of 1 (positive), but the true label is 0.
        FP = np.sum(np.logical_and(x == 1, y == 0), axis=axis)
 
        # False Negative (FN): we predict a label of 0 (negative), but the true label is 1.
        FN = np.sum(np.logical_and(x == 0, y == 1),  axis=axis)

        #go ahead and calculate macro F-score if on horizontal axis
        if(axis==1):
            TP0=TP
            TP = np.add(TP, eps)
            precision = np.divide(TP, np.add(TP, FP))
            recall = np.divide(TP, np.add(TP, FN))
            top = np.multiply(precision, recall)
            bottom = np.add(precision, recall)

            #macro F-score, per feature, no weights
            macro_f1 = np.multiply(2.0, np.divide(top,bottom))

            results['macro_f1']=macro_f1
            results['var_TPr']=np.round(TP/(TP+FP), decimals=3)
            results['var_FPr']=np.round(FP/(TP+FP), decimals=3)
            results['var_FNr']=np.round(FN/(TN+FN), decimals=3)
            results['var_TNr']=np.round(TN/(TN+FN), decimals=3)
            results['var_TP']=TP0
            results['var_FP']=FP
            results['var_FN']=FN
            results['var_TN']=TN
            results['var_precision']=np.round(precision, decimals=3)
            results['var_recall']=np.round(recall, decimals=3)

            if(X_mode!="False"):
                weights = np.sum(y_in, axis=axis)
                results['top'] = top
                results['precision'] = precision
                results['recall'] = recall
                results['weights'] = weights

            #weights = np.divide(weights,np.sum(weights))
            #weights = np.interp(weights, (weights.min(), weights.max()), (0.0, 1.0))

            #top = np.multiply(np.add(weights,1.0), top)
            #bottom = np.add(np.multiply(weights,precision), recall)
            #w_macro_f1 = np.divide(top,bottom)
            #results['w_macro_f1'] = w_macro_f1
            #w_macro_f1 = np.round(w_macro_f1, decimals=3)

            #macro_f1 = np.round(macro_f1, decimals=3)

        else:
            TP = TP.reshape(int(len(TP)/2),2)
            FP = FP.reshape(int(len(FP)/2),2)
            FN = FN.reshape(int(len(FN)/2),2)
            TN = TN.reshape(int(len(TN)/2),2)
            #TP = TP.reshape(int(len(TP)/3),3)
            #FP = FP.reshape(int(len(FP)/3),3)
            #FN = FN.reshape(int(len(FN)/3),3)
            #TN = TN.reshape(int(len(TN)/3),3)
            TP = np.sum(TP, axis=1)
            FP = np.sum(FP, axis=1)
            FN = np.sum(FN, axis=1)
            TN = np.sum(TN, axis=1)
            #true positives per sample
            results['TP']=TP
            #false positives per sample
            results['FP']=FP
            #true negatives per sample
            results['FN']=FN
            results['TN']=TN
            #sample size
            results['N']=len(y)*2
            #results['N']=len(y)*3

    return results

def accuracy_ratio(x_in,y_in):

    #Some people from genomics call this "concordance rate"
    #but that is just a fancy name for the old accuracy ratio used in machine learning
    #which is basically the number of correct predictions divided by the total number of predictions

    eps=1e-15

    dose_to_probability = {0: [1,0,0], 1: [0,1,0], 2: [0,0,1], -1: [0,0,0]}

    x, y = convert_and_reshape(dose_to_probability, x_in,y_in)

    results=dict()

    correct_prediction = np.equal( x, y )

    accuracy = np.mean(correct_prediction.astype(float), axis=1)
    correct_pred_per_sample = np.sum(correct_prediction.astype(float), axis=0)

    correct_pred_per_sample = correct_pred_per_sample.reshape(int(len(correct_pred_per_sample)/3),3)
    correct_pred_per_sample = np.sum(correct_pred_per_sample, axis=1)
    correct_pred_per_sample = np.divide(correct_pred_per_sample,3.0)

    #xy=np.multiply(x,y)

    #xy_per_sample=np.sum(xy,axis=0)

    #xy_per_sample=xy_per_sample.reshape(int(xy_per_sample.shape[0]/3),-1).sum(1)

    #P0 = np.divide(np.sum(xy,axis=1), len(y[0])/3)

    p11 = np.multiply(y[:,0::3],x[:,0::3])
    p11s = np.sum(p11,axis=0)
    p11 = np.sum(p11,axis=1)
    p22 = np.multiply(y[:,1::3],x[:,1::3])
    p22s = np.sum(p22,axis=0)
    p22 = np.sum(p22,axis=1)
    p33 = np.multiply(y[:,2::3],x[:,2::3])
    p33s = np.sum(p33,axis=0)
    p33 = np.sum(p33,axis=1)


    xy_per_sample=np.add(p11s,np.add(p22s,p33s))

    xy = np.add(p11,np.add(p22,p33))

    P0 = np.divide(xy, np.divide(len(y[0]),3.0))

    p12 = np.multiply(y[:,0::3],x[:,1::3])
    p12 = np.sum(p12,axis=1)

    p21 = np.multiply(y[:,1::3],x[:,0::3])
    p21 = np.sum(p21,axis=1)

    p13 = np.multiply(y[:,0::3],x[:,2::3])
    p13 = np.sum(p13,axis=1)

    p31 = np.multiply(y[:,2::3],x[:,0::3])
    p31 = np.sum(p31,axis=1)

    p23 = np.multiply(y[:,1::3],x[:,2::3])
    p23 = np.sum(p23,axis=1)

    p32 = np.multiply(y[:,2::3],x[:,1::3])
    p32 = np.sum(p32,axis=1)

    Ns=len(y[0])/3

    #p11 = np.sum(xy[:,0::3], axis=1)
    #p22 = np.sum(xy[:,1::3], axis=1)
    #p33 = np.sum(xy[:,2::3], axis=1)

    N1_i = np.add(np.add(p11, p21), p31)
    N2_i = np.add(np.add(p12, p22), p32)
    N3_i = np.add(np.add(p13, p23), p33)

    N1_j = np.add(np.add(p11, p12), p13)
    N2_j = np.add(np.add(p21, p22), p23)
    N3_j = np.add(np.add(p31, p32), p33)

    N1_ij = np.multiply(N1_i, N1_j)
    N2_ij = np.multiply(N2_i, N2_j)
    N3_ij = np.multiply(N3_i, N3_j)

    num = np.add(N1_ij, np.add(N2_ij, N3_ij))
    den = np.power(Ns,2)

    #num = np.add(eps,num)
    #den = np.add(eps,den)

    Pc = np.divide(num,den)

    num = np.subtract(P0,Pc)
    den = np.subtract(1.0,Pc)

    IQS = np.divide(num,den)
    IQS = np.nan_to_num(IQS, 0)

    results['accuracy_per_var'] = accuracy
    results['correct_pred_per_sample'] = correct_pred_per_sample
    results['xy_per_sample'] = xy_per_sample
    results['P0_per_var'] = P0
    results['IQS'] = IQS
    results['N'] = len(y)

    if(DEBUG==True):
        print(results)

    return results

def pearson_r2(x_in, y_in):

    results=dict()
    r2_results_l = []
    p_results = []

    '''
    for pos in x_in.keys():
        r2,p = linregress(x_in[pos], y_in[pos])[2:4]
        r2=r2**2
        r2_results_l.append(np.round(r2, decimals=3))
        p_results.append(np.round(p, decimals=3))
    '''

    #x = np.asarray(list(x_in.values()), dtype=float)
    #y = np.asarray(list(y_in.values()), dtype=float)
    x = np.copy(x_in)
    y = np.copy(y_in)
    if(DEBUG==True):
        np.savetxt("x_before", x)    
        np.savetxt("y_before", y)    

    #per variant
    x_sum = np.sum(x, axis=1)
    y_sum = np.sum(y, axis=1)
    xy_sum = np.sum(np.multiply(x,y), axis=1)
    x_squared_sum = np.sum(np.power(x,2), axis=1)
    y_squared_sum = np.sum(np.power(y,2), axis=1)
    N=len(y[0])

    num=np.subtract(np.multiply(xy_sum, N), np.multiply(x_sum, y_sum) )
    den=np.multiply(x_squared_sum, N)
    den=np.subtract(den, np.power(x_sum,2))
    den2=np.multiply(y_squared_sum, N)
    den2=np.subtract(den2, np.power(y_sum,2))
    den=np.sqrt(np.multiply(den, den2))

    eps=1e-15

    #fixed bug that would print correl = 1 when WGS_MAF is zero
    #num = np.add(eps,num)
    #den = np.add(eps,den)

    r2_per_variant=np.divide(num,den)
    r2_per_variant=np.nan_to_num(r2_per_variant,0)
    r2_per_variant=np.power(r2_per_variant,2)
    r2_per_variant=np.round(r2_per_variant, decimals=3)

    #per sample
    x_sum = np.sum(x, axis=0)
    y_sum = np.sum(y, axis=0)
    xy=np.multiply(x,y)
    xy_sum = np.sum(xy, axis=0)
    x_squared_sum = np.sum(np.power(x,2), axis=0)
    y_squared_sum = np.sum(np.power(y,2), axis=0)

    #results['r2_per_variant_linregress']=r2_results_l
    results['r2_per_variant_manual']=r2_per_variant
    #results['p_per_variant']=p_results

    results['x_sum']=x_sum
    results['y_sum']=y_sum
    results['xy_sum']=xy_sum
    results['x_squared_sum']=x_squared_sum
    results['y_squared_sum']=y_squared_sum

       
    if(DEBUG==True):
        results['xy']=xy
        np.savetxt("x_after", x)    
        np.savetxt("y_after", y) 
    #results['mean_x_per_sample']=mean_x
    #results['mean_y_per_sample']=mean_y
    #results['var_y_per_sample']=var_y
    #results['var_x_per_sample']=var_x
    #results['covar']=covar
    results['N']=len(y)

    return results

def extract_ga_positions(coordinates):

    GA = VCF(GA_file)
    result=[]

    for variant in GA(coordinates):
        result.append(str(variant.CHROM)+':'+str(variant.end))

    if(len(result)==0):
        return False

    return result

def extract_vcf_lines(coordinates):

    if WGS_file_bcf!="":
        WGS = VCF(WGS_file_bcf)
    else:
        WGS = VCF(WGS_file)

    lines=[]

    for variant in WGS(coordinates):
       lines.append(str(variant))

    if(len(lines))>0:
        return lines
    else:
        return False

def extract_sample_ids_from_vcf(input_file):

    sample_ids=False

    if input_file.endswith(".gz"):
        file = gzip.open(input_file,'rb')
    else:
        file = open(input_file)

    while True:
        line=file.readline()

        if input_file.endswith(".gz"):
            line=line.decode('utf-8')
            #line=line.split('\n')
        if(line[0]=='#'):
            if(line.startswith("#CHROM")):
                sample_ids=line.split('\t')[9:]
                sample_ids[-1]=sample_ids[-1].replace('\n','')
                break
        else:
            break

    file.close()

    return sample_ids

def convert_gt_to_int(gt):

    genotype_to_int={'0/0': 0, '0|0': 0.0, '0/1': 1, '0|1':1, '1/0':1, '1|0':1, '1/1':2, '1|1':2, './0':-1, './1':-1, './.':-1, '0/.':-1, '1/.':-1}

    result=genotype_to_int[gt[0:3]]

    return result

def extract_dose_from_line(vcf_line, disable_DS=False):

    vcf_line=vcf_line.split('\t')

    vcf_line[-1]=vcf_line[-1].replace('\n','')

    result_line=[]
    pos=vcf_line[0]+':'+vcf_line[1]
    #pos = chr12:13412322
    result_line.append(pos)

    snp_id=pos+'_'+vcf_line[3]+'_'+vcf_line[4]
    #snp_id = chr12:13412322_A_T
    for column in vcf_line[9:]:
        if(':' in column and disable_DS==False):
            gen_dose=column.split(':')
            result=gen_dose[1]
            result=float(result)
        else:
            result=convert_gt_to_int(column[0:3])
        result_line.append(result)

    return result_line, snp_id

def remove_positions(pos_list, dict):

    for pos in pos_list:
        dict.pop(pos, None)

def intersect_positions(dict_1,dict_2):

   s1 = set(dict_1)
   s2 = set(dict_2)
   common_keys = s1 & s2
   subdict_1 = {x: dict_1[x] for x in common_keys if x in dict_1}
   subdict_2 = {x: dict_2[x] for x in common_keys if x in dict_2}

   sorted_subdict_1 = dict(sorted(subdict_1.items(), key=lambda item: item[0]))
   sorted_subdict_2 = dict(sorted(subdict_2.items(), key=lambda item: item[0]))

   #return subdict_1, subdict_2
   return sorted_subdict_1, sorted_subdict_2

def calculate_MAF(input_file, coordinates):
    #plink --vcf HRC.r1-1.EGA.GRCh37.chr9.haplotypes.9p21.3.vcf.clean4 --freq
    #out_name=input_file.split('.gz')[0]

    cmd=f"bcftools +fill-tags {input_file} -r {coordinates} -- -t AF | bcftools query -f '%CHROM:%POS\t%AF\n'"

    #result = sp.check_output(cmd, encoding='UTF-8', shell=True)
    result = sp.check_output(cmd, encoding='UTF-8', shell=True)
    #result = sp.check_output("plink --vcf "+input_file+" --freq --out "+out_name, encoding='UTF-8', shell=True)
    maf_result={}
    result=result.split('\n')

    for line in result:
        line=line.split('\t')
        if(len(line)>1):
            maf=float(line[1])
            if(maf>0.5):
                maf=str(1-maf)
                line[1]=maf
            maf_result[line[0]]=line[1]

    #result = sp.check_output("plink --vcf "+input_file+" --freq --out "+out_name, encoding='UTF-8', shell=True)

    #maf_result = read_MAF_file(out_name+".frq")

    return maf_result

def read_MAF_file(frq_file):

    #   CHR         SNP   A1   A2          MAF  NCHROBS
    #   9   rs1333039    G    C       0.3821    54330
    #   9 rs138885769    T    C    0.0008099    54330
    #   9 rs548022918    T    G     0.000589    54330
    pd.set_option('display.float_format', '{:.6f}'.format)

    maftable = pd.read_csv(frq_file, sep='\s+', comment='#')
    maftable['MAF'] = maftable['MAF'].astype(float)

    maf = maftable['MAF'].values.tolist()
    snp = maftable['SNP'].values.tolist()

    seen = set()
    for i in snp:
        if(i in seen):
            print("WARNING: repeated SNP ID \"",i,"\". MAF results may not be reliable, please give unique names to your WGS vcf.")
        seen.add(i)

    result = dict(zip(snp, maf))

    return result

def process_lines(lines):

    results=[]
    result_line=[]
    pos=0
    chr=0
    start_pos='0'
    first=True

    imputed_dosages={}
    wgs_dosages={}
    ga_positions=[]
    snp_ids={}

    for line in lines:
        #skip comments
        if(line[0]=='#'):
            continue
        if(first==True):
            chr, start_pos = line.split('\t')[0:2]
            first=False

        pos_dosages_tmp, snp = extract_dose_from_line(line,disable_DS)
        imputed_dosages[pos_dosages_tmp[0]] = pos_dosages_tmp[1:]

    end_pos = lines[-1].split('\t')[1]

    coordinates = chr+':'+start_pos+'-'+end_pos

    if(multimask_mode==False):
        ga_pos = extract_ga_positions(coordinates)
    else:
        ga_pos = False

    wgs_lines = extract_vcf_lines(coordinates)

    if(wgs_lines==False):
        return [False]

    for line in wgs_lines:
        pos_dosages_tmp, snp = extract_dose_from_line(line,disable_DS)
        wgs_dosages[pos_dosages_tmp[0]] = pos_dosages_tmp[1:]
        snp_ids[pos_dosages_tmp[0]]=snp
    

    if(ga_pos!=False and multimask_mode==False):
        remove_positions(ga_pos, imputed_dosages)

        if(wgs_lines!=False):
            remove_positions(ga_pos, wgs_dosages)

    if(len(imputed_dosages)==0 or len(wgs_dosages)==0):
        return [False]

    imputed_dosages, wgs_dosages = intersect_positions(imputed_dosages, wgs_dosages)
    snp_ids, wgs_dosages = intersect_positions(snp_ids, wgs_dosages)
    snp_ids, imputed_dosages = intersect_positions(snp_ids,imputed_dosages)

    if(len(imputed_dosages)==0 or len(wgs_dosages)==0):
        #print(ga_pos)
        return [False]

    maf_dict={}
    wgs_maf_dict={}
    imp_maf_dict={}
    #pos->maf
    if(REF_file!=""):
        maf_dict=calculate_MAF(REF_file, coordinates)
        wgs_maf_dict=calculate_MAF(WGS_file, coordinates)
        wgs_maf_dict, snp_ids = intersect_positions(wgs_maf_dict, snp_ids)

    else:
        if(WGS_file_bcf!=""):
            maf_dict=calculate_MAF(WGS_file_bcf, coordinates)
        else:
            maf_dict=calculate_MAF(WGS_file, coordinates)

    maf_dict, snp_ids = intersect_positions(maf_dict, snp_ids)

    if(IMPUTED_file_bcf!=""):
        imp_maf_dict=calculate_MAF(IMPUTED_file_bcf, coordinates)
    else:
        imp_maf_dict=calculate_MAF(IMPUTED_file, coordinates)

    imp_maf_dict, snp_ids = intersect_positions(imp_maf_dict, snp_ids)

    f1_dict = f1_score(list(imputed_dosages.values()), list(wgs_dosages.values()))
    acc_dict = accuracy_ratio(list(imputed_dosages.values()), list(wgs_dosages.values()))
    r2_dict = pearson_r2(list(imputed_dosages.values()), list(wgs_dosages.values()))
    rmse_dict = rmse(list(imputed_dosages.values()), list(wgs_dosages.values()))

    if(DEBUG==True):
        w = csv.writer(open("imputed_dosages.csv", "w"))
        for key, val in imputed_dosages.items():
            w.writerow([key, val])
        w = csv.writer(open("wgs_dosages.csv", "w"))
        for key, val in wgs_dosages.items():
            w.writerow([key, val])     
       
    #acc_dict['accuracy_per_var']
    #acc_dict['correct_pred_per_sample']

    #f-score per var
    #f1_dict['macro_f1']
    #true positives per sample
    #f1_dict['TP']
    #false positives per sample
    #f1_dict['FP']
    #true negatives per sample
    #f1_dict['FN']
    #sample size per chunk
    #print(snp_ids)

    #return list(imputed_dosages.keys()), f1_dict['macro_f1'], acc_dict['accuracy_per_var'], acc_dict['correct_pred_per_sample'], acc_dict['N'], f1_dict['TP'], f1_dict['FN'], f1_dict['N']

    return list(imputed_dosages.keys()), f1_dict, acc_dict, r2_dict, snp_ids, maf_dict, imp_maf_dict, wgs_maf_dict, rmse_dict

def load_file_chunks(ncores):

    imputed_sample_ids=extract_sample_ids_from_vcf(IMPUTED_file)
    WGS_sample_ids=extract_sample_ids_from_vcf(WGS_file)

    #print(len(imputed_sample_ids), len(WGS_sample_ids))

    #print('Processing ',len(imputed_sample_ids),'imputed samples')

    read_total=0
    chunk_total=0
    calc_total=0
    data = []
    results = []
    ci=0
    cj=0
    if IMPUTED_file.endswith(".gz"):
        file = gzip.open(IMPUTED_file,'rb')
    else:
        file = open(IMPUTED_file)

    while True:
    #for line in file:
        read_start = timeit.default_timer()

        line=file.readline()
        if IMPUTED_file.endswith(".gz"):
            line=line.decode('utf-8')
            #line=line.split('\n')

        if(line!=""):
            if(line[0]!='#'):
                data.append(line)

        read_stop = timeit.default_timer()
        read_total += read_stop-read_start

        if(len(data)>=max_total_rows or line == ""):
            cj+=1
            #print("Processing chunk:",cj,"Max rows per chunk:",max_total_rows)
            chunk_start = timeit.default_timer()
            data=chunk(data)
            chunk_stop = timeit.default_timer()
            chunk_total += chunk_stop-chunk_start

            calc_start = timeit.default_timer()
            
            pool = mp.Pool(ncores)
            results_tmp = pool.map(process_lines,data)
            pool.close()
            pool.join()

            data=[]
            results_tmp = [item for sublist in results_tmp for item in sublist if len(sublist)>1]
            if(len(results_tmp)>0):
                results.append(results_tmp)
                ci+=1
            if(DEBUG!=False):
                print("RESULT length: ", len(results_tmp), "ci:", ci)
                print("RESULT 0: ", len(results_tmp[0]), "ci:", ci)
                print("results_tmp", results_tmp)
            calc_stop = timeit.default_timer()
            calc_total += calc_stop-calc_start

        if(line==""):
            break

    file.close()

    #list(imputed_dosages.keys()), f1_dict, acc_dict

   #imputed_dosages.keys(), 
   #f1_dict['macro_f1'], acc_dict['accuracy_per_var'], acc_dict['correct_pred_per_sample'], acc_dict['N'], 
   #f1_dict['TP'], f1_dict['FN'], f1_dict['N']

    merge_start = timeit.default_timer()
    if(ci==0):
        print("After removing genotype array variants from imputed data, and intersecting WGS, no variants remained")
        exit(0)

    #maf_start = timeit.default_timer()
    #if(REF_file==""):
    #    maf_dict = calculate_MAF(WGS_file)
    #else:
    #    maf_dict = calculate_MAF(REF_file)
    #
    #maf_stop = timeit.default_timer()

    chunk_i=0

    pos=[]
    accuracy_per_var=[]
    macro_f1=[]
    weights=[]
    top=[]
    precision=[]
    recall=[]
    ref_mafs=[]
    wgs_mafs=[]
    imp_mafs=[]
    snp_ids=[]

    P0_per_var=[]
    IQS=[]

    vTP = []
    vTN = []
    vFP = []
    vFN = []
    vTPr = []
    vTNr = []
    vFPr = []
    vFNr = []
    var_precision = []
    var_recall = []
    vRMSE = []

    #positions: index 0
    #f1 result dictionary: index 1
    #accuracy dictionary: 2
    #r2 dictionary: 3
    #snp ids: 4
    #ref/wgs maf: 5
    #imputed maf: 6
    #wgs maf: 7
    pi=0
    fi=1
    ai=2
    ri=3
    vi=4
    mi=5
    mi2=6
    mi3=7
    msi=8 #rsme index
    #total number of items per chunk
    ni=9


    #    results['accuracy_per_var'] = accuracy
    #    results['correct_pred_per_sample'] = correct_pred_per_sample
    #    results['xy_per_sample'] = xy_per_sample
    #    results['P0_per_var'] = P0
    #    results['IQS'] = IQS
    #    results['N'] = len(y)


    for chunk_i in range(ci):
        #print('results shape:', len(results), len(results[chunk_i]),len(results[chunk_i][1]))
        for si in list(range(len(results[chunk_i])))[0::ni]:
            pos_tmp=list(results[chunk_i][si+pi])
            pos.extend(pos_tmp)
            snp_ids_tmp=list(results[chunk_i][si+vi].values())
            mafs_tmp=list(results[chunk_i][si+mi].values())
            if(REF_file!=""):
                ref_mafs.extend(mafs_tmp)
                mafs_tmp=list(results[chunk_i][si+mi3].values())
            wgs_mafs.extend(mafs_tmp)
            mafs_tmp=list(results[chunk_i][si+mi2].values())
            imp_mafs.extend(mafs_tmp)

            snp_ids.extend(snp_ids_tmp)
            results[chunk_i][si+pi] = 0
            results[chunk_i][si+vi] = 0
            results[chunk_i][si+mi] = 0
            accuracy_per_var.extend(list(results[chunk_i][si+ai]['accuracy_per_var']))
            results[chunk_i][si+ai]['accuracy_per_var'] = 0
            P0_per_var.extend(list(results[chunk_i][si+ai]['P0_per_var']))
            results[chunk_i][si+ai]['P0_per_var'] = 0
            IQS.extend(list(results[chunk_i][si+ai]['IQS']))
            results[chunk_i][si+ai]['IQS'] = 0
            macro_f1.extend(list(results[chunk_i][si+fi]['macro_f1']))
            results[chunk_i][si+fi]['macro_f1'] = 0

            vTP.extend(list(results[chunk_i][si+fi]['var_TP']))
            results[chunk_i][si+fi]['var_TP'] = 0
            vTN.extend(list(results[chunk_i][si+fi]['var_TN']))
            results[chunk_i][si+fi]['var_TN'] = 0
            vFP.extend(list(results[chunk_i][si+fi]['var_FP']))
            results[chunk_i][si+fi]['var_FP'] = 0
            vFN.extend(list(results[chunk_i][si+fi]['var_FN']))
            results[chunk_i][si+fi]['var_FN'] = 0

            vTPr.extend(list(results[chunk_i][si+fi]['var_TPr']))
            results[chunk_i][si+fi]['var_TPr'] = 0
            vTNr.extend(list(results[chunk_i][si+fi]['var_TNr']))
            results[chunk_i][si+fi]['var_TNr'] = 0
            vFPr.extend(list(results[chunk_i][si+fi]['var_FPr']))
            results[chunk_i][si+fi]['var_FPr'] = 0
            vFNr.extend(list(results[chunk_i][si+fi]['var_FNr']))
            results[chunk_i][si+fi]['var_FNr'] = 0

            var_precision.extend(list(results[chunk_i][si+fi]['var_precision']))
            results[chunk_i][si+fi]['var_precision'] = 0
            var_recall.extend(list(results[chunk_i][si+fi]['var_recall']))
            results[chunk_i][si+fi]['var_recall'] = 0
            
            vRMSE.extend(list(results[chunk_i][si+msi]['var_rmse']))

            if(X_mode!="False"):
                weights.extend(list(results[chunk_i][si+fi]['weights']))
                results[chunk_i][si+fi]['weights'] = 0
                top.extend(list(results[chunk_i][si+fi]['top']))
                results[chunk_i][si+fi]['top'] = 0
                precision.extend(list(results[chunk_i][si+fi]['precision']))
                results[chunk_i][si+fi]['precision'] = 0
                recall.extend(list(results[chunk_i][si+fi]['recall']))
                results[chunk_i][si+fi]['recall'] = 0

    #print('pos',pos[0:11], len(pos))
    #print('accuracy_per_var',accuracy_per_var[0:11], len(accuracy_per_var))
    #print('macro_f1',macro_f1[0:11], len(macro_f1))
    vRMSE=np.round(vRMSE, decimals=3)

    if(X_mode!="False"):
        weights = np.divide(weights,np.sum(weights))
        weights = np.interp(weights, (weights.min(), weights.max()), (0.0, 1.0))
        top = np.multiply(np.add(weights,1.0), top)
        bottom = np.add(np.multiply(weights,precision), recall)
        w_macro_f1 = np.divide(top,bottom)
        w_macro_f1 = np.round(w_macro_f1, decimals=3)
    macro_f1 = np.round(macro_f1, decimals=3)



    CPT=np.zeros(len(results[chunk_i][ai]['correct_pred_per_sample']))
    #variant count
    CNT=0

    xy_per_sample=np.zeros(len(results[chunk_i][ai]['xy_per_sample']))

    for chunk_i in range(ci):
        for si in list(range(len(results[chunk_i])))[0::ni]:
            CPT=np.add(CPT, results[chunk_i][si+ai]['correct_pred_per_sample'])
            xy_per_sample=np.add(xy_per_sample, results[chunk_i][si+ai]['xy_per_sample'])
            results[chunk_i][si+ai]['correct_pred_per_sample'] = 0
            CNT+=results[chunk_i][si+ai]['N']

    #print('results', len(results[chunk_i][1]['correct_pred_per_sample']) )
    #print('CPT', CPT, len(CPT))
    #print('CNT', CNT)
    #print('chunks', ci)

    accuracy_per_sample=np.divide(CPT,CNT)
    accuracy_per_sample=np.round(accuracy_per_sample, decimals=3)
    P0_per_sample=np.divide(xy_per_sample,CNT)
    P0_per_sample=np.round(P0_per_sample, decimals=3)


    digits = len(str(len(P0_per_sample)))+1
    imp_mafs=np.round(np.asfarray(imp_mafs,float), decimals=digits)
    if(REF_file!=""):
        ref_mafs=np.round(np.asfarray(ref_mafs,float), decimals=digits)
    wgs_mafs=np.round(np.asfarray(wgs_mafs,float), decimals=digits)

    CPT=0
    xy_per_sample=0
    #print('accuracy_per_sample', accuracy_per_sample, len(accuracy_per_sample))

    eps = 1e-8
    chunk_i=0
    si=0

    if(DEBUG!=False):
        print("chunk_i", chunk_i, "si", si, "fi" , fi)
        print("len(results[chunk_i])", len(results[chunk_i]))
        print("len(results[chunk_i][si+fi])", len(results[chunk_i][si+fi]) )        

    N = 0
    TP = np.zeros(len(results[chunk_i][si+fi]['TP']))
    FP = np.zeros(len(results[chunk_i][si+fi]['FP']))
    FN = np.zeros(len(results[chunk_i][si+fi]['FN']))
    TN = np.zeros(len(results[chunk_i][si+fi]['TN']))

    ses = np.zeros(len(results[chunk_i][si+msi]['sample_ses']))
    #results['var_rmse'] = var_rmse
    #results['sample_ses'] = sample_ses
    #results['N'] = sample_ses
    Nses=0

    for chunk_i in range(ci):
        for si in list(range(len(results[chunk_i])))[0::ni]:
            TP=np.add(TP,results[chunk_i][si+fi]['TP'])
            results[chunk_i][si+fi]['TP'] = 0
            FP=np.add(FP,results[chunk_i][si+fi]['FP'])
            results[chunk_i][si+fi]['FP'] = 0
            FN=np.add(FN,results[chunk_i][si+fi]['FN'])
            results[chunk_i][si+fi]['FN'] = 0
            TN=np.add(TN,results[chunk_i][si+fi]['TN'])
            results[chunk_i][si+fi]['TN'] = 0
            N+=results[chunk_i][si+fi]['N']

            ses=np.add(ses,results[chunk_i][si+msi]['sample_ses'])
            results[chunk_i][si+msi]['sample_ses'] = 0
            Nses+=results[chunk_i][si+msi]['N']


    TP0=TP
    TP = np.add(TP, eps)
    precision = np.divide(TP, np.add(TP, FP))
    recall = np.divide(TP, np.add(TP, FN))
    top = np.multiply(precision, recall)
    bottom = np.add(precision, recall)

    per_sample_f1 = np.multiply(2.0, np.divide(top,bottom))
    per_sample_f1 = np.round(per_sample_f1, decimals=3)
    
    sTPr = TP/(TP+FP)
    sFPr = FP/(TP+FP)
    sFNr = FN/(TN+FN)
    sTNr = TN/(TN+FN)
    sTP = TP0
    sFP = FP
    sFN = FN
    sTN = TN

    sTPr = np.round(sTPr, decimals=3)
    sTNr = np.round(sTNr, decimals=3)
    sFPr = np.round(sFPr, decimals=3)
    sFNr = np.round(sFNr, decimals=3)

    precision = np.round(precision, decimals=3)
    recall = np.round(recall, decimals=3)

    TP = 0
    FP = 0
    FN = 0
    TN = 0

    sRMSE=np.sqrt(ses/Nses)
    sRMSE=np.round(sRMSE, decimals=3)

    #print('per_sample_f1', per_sample_f1, len(per_sample_f1))

    #results['x_sum']=x_sum
    #results['y_sum']=y_sum
    #results['xy_sum']=xy_sum
    #results['x_squared_sum']=x_squared_sum
    #results['y_squared_sum']=y_squared_sum

    #results['r2_per_variant']=r2_results
    #results['p_per_variant']=p_results

    chunk_i=0
    si=0

    x_sum=np.zeros(len(results[chunk_i][si+ri]['x_sum']))
    y_sum=np.zeros(len(results[chunk_i][si+ri]['y_sum']))
    xy_sum=np.zeros(len(results[chunk_i][si+ri]['xy_sum']))
    x_squared_sum=np.zeros(len(results[chunk_i][si+ri]['x_squared_sum']))
    y_squared_sum=np.zeros(len(results[chunk_i][si+ri]['y_squared_sum']))
    N=0

    r2_per_variant=[]
    r2_per_variant_m=[]
    p_per_variant=[]
    for chunk_i in range(ci):
        for si in list(range(len(results[chunk_i])))[0::ni]:
            x_sum=np.add(x_sum,results[chunk_i][si+ri]['x_sum'])
            results[chunk_i][si+ri]['x_sum']=0
            y_sum=np.add(y_sum,results[chunk_i][si+ri]['y_sum'])
            results[chunk_i][si+ri]['y_sum']=0
            if(DEBUG==True):
                #print(np.savetxt(sys.stdout, results[chunk_i][si+ri]['xy']))
                print(np.savetxt("xy", results[chunk_i][si+ri]['xy']))

                print("DEBUG r2: cumulative xy_sum", xy_sum, "original", results[chunk_i][si+ri]['xy_sum'])
                print("DEBUG r2: xy_sum", results[chunk_i][si+ri]['xy_sum'], "len(xy_sum)", len(results[chunk_i][si+ri]['xy_sum']))
                print("DEBUG r2: xy", results[chunk_i][si+ri]['xy'], "len(xy)", len(results[chunk_i][si+ri]['xy']))
            xy_sum=np.add(xy_sum,results[chunk_i][si+ri]['xy_sum'])
            results[chunk_i][si+ri]['xy_sum']=0
            x_squared_sum=np.add(x_squared_sum,results[chunk_i][si+ri]['x_squared_sum'])
            results[chunk_i][si+ri]['x_squared_sum']=0
            y_squared_sum=np.add(y_squared_sum,results[chunk_i][si+ri]['y_squared_sum'])
            results[chunk_i][si+ri]['u_squared_sum']=0
            N+=results[chunk_i][si+ri]['N']
            #r2_per_variant.extend(results[chunk_i][si+ri]['r2_per_variant_linregress'])
            #results[chunk_i][si+ri]['r2_per_variant_linregress']=0
            r2_per_variant_m.extend(results[chunk_i][si+ri]['r2_per_variant_manual'])
            results[chunk_i][si+ri]['r2_per_variant_manual']=0
            #p_per_variant.extend(results[chunk_i][si+ri]['p_per_variant'])
            #results[chunk_i][si+ri]['p_per_variant']=0

    #np.subtract(np.multiply(xy_sum, N), np.multiply(x_sum, y_sum) )
    num=np.subtract(np.multiply(xy_sum, N), np.multiply(x_sum, y_sum) )
    den=np.multiply(x_squared_sum, N)
    den=np.subtract(den, np.power(x_sum,2))
    den2=np.multiply(y_squared_sum, N)
    den2=np.subtract(den2, np.power(y_sum,2))
    den=np.sqrt(np.multiply(den, den2))
    r2_per_sample=np.divide(num,den)
    r2_per_sample=np.power(r2_per_sample,2)
    if(DEBUG==True):
        print("DEBUG r2: num", num,"xy_sum", xy_sum, "N", N, "x_squared_sum", x_squared_sum,"x_sum",x_sum, "den", den,"den2", den2, "r2_per_sample", r2_per_sample)
        

    r2_per_sample=np.round(r2_per_sample, decimals=3)
    #accuracy_per_var=np.round(accuracy_per_var, decimals=3)
    r2_per_variant=np.round(r2_per_variant, decimals=3)

    IQS = np.round(IQS, decimals=3)
    IQS = np.nan_to_num(IQS, 0)

    P0_per_var = np.round(P0_per_var, decimals=3)

    #merged_results_per_sample=np.column_stack((imputed_sample_ids,WGS_sample_ids,per_sample_f1,P0_per_sample,r2_per_sample))

    labels_s=['imputed_ids','WGS_ids','F-score','concordance_P0','r2', 'precision', 'recall', 'TP', 'TN', 'FP', 'FN', 'TP_ratio', 'TN_ratio', 'FP_ratio', 'FN_ratio', 'RMSE']
    merged_results_per_sample=np.column_stack((imputed_sample_ids,WGS_sample_ids,per_sample_f1,P0_per_sample,r2_per_sample, precision, recall, sTP, sTN, sFP, sFN, sTPr, sTNr, sFPr, sFNr, sRMSE))

    labels_v=[]
    if(REF_file!=""):
        merged_results_per_variant=np.column_stack((pos,snp_ids,ref_mafs,imp_mafs,wgs_mafs,macro_f1,P0_per_var, IQS, r2_per_variant_m, var_precision, var_recall, vTP, vTN, vFP, vFN, vTPr, vTNr, vFPr, vFNr, vRMSE))
        labels_v=['position','SNP','REF_MAF','IMPUTED_MAF','WGS_MAF', 'F-score', 'concordance_P0','IQS', 'r2', 'precision', 'recall', 'TP', 'TN', 'FP', 'FN', 'TP_ratio', 'TN_ratio', 'FP_ratio', 'FN_ratio', 'RMSE']
    else:
        merged_results_per_variant=np.column_stack((pos,snp_ids,imp_mafs,wgs_mafs,macro_f1,P0_per_var, IQS, r2_per_variant_m, var_precision, var_recall, vTP, vTN, vFP, vFN, vTPr, vTNr, vFPr, vFNr, vRMSE))
        labels_v=['position','SNP','IMPUTED_MAF','WGS_MAF', 'F-score', 'concordance_P0','IQS', 'r2', 'precision', 'recall', 'TP', 'TN', 'FP', 'FN', 'TP_ratio', 'TN_ratio', 'FP_ratio', 'FN_ratio', 'RMSE']

    if(X_mode!="False"):
         merged_results_per_variant=np.column_stack((merged_results_per_variant,accuracy_per_var,w_macro_f1))
         labels_v.extend(['W_F-score', 'accuracy_ratio'])
         merged_results_per_sample=np.column_stack((merged_results_per_sample,accuracy_per_sample))
         labels_s.extend(['accuracy_ratio'])


    merged_results_per_variant=np.vstack((np.asarray(labels_v), merged_results_per_variant))

    merged_results_per_sample=np.vstack((np.asarray(labels_s),merged_results_per_sample))


    outname=IMPUTED_file.split('.gz')[0]
    global vout
    global sout

    if(vout==""):
        vout=outname+'_per_variant_results.txt'
    if(sout==""):
        sout=outname+'_per_sample_results.txt'
    #np.savetxt(outname+'_per_sample_results.txt', merged_results_per_sample, delimiter="\t")
    #merged_results_per_sample.tofile(outname+'_per_sample_results.txt', sep='\t', format='%10.5f')
    #merged_results_per_variant.tofile(outname+'_per_variant_results.txt', sep='\t', format='%10.5f')
    #np.savetxt(outname+'_per_variant_results.txt', merged_results_per_variant, delimiter="\t")
    np.savetxt(sout, merged_results_per_sample, fmt="%s", delimiter="\t")
    np.savetxt(vout, merged_results_per_variant, fmt="%s", delimiter="\t")

    #print('r2_per_sample', r2_per_sample[0:11], len(r2_per_sample))
    #print('r2_per_variant', r2_per_variant[0:11], len(r2_per_variant))
    #print('p_per_variant', p_per_variant[0:11], len(p_per_variant))

    merge_stop = timeit.default_timer()

    merge_total = merge_stop-merge_start

    #maf_total = maf_stop-maf_start

    #print('1 Read imputed file time: ', read_total)

    #print('2 Chunking time: ', chunk_total)

    #print('3 Calculation time: ', calc_total)

    #print('4 Merging calculations per sample time: ', merge_total)

    #print('5 MAF calculation time: ', maf_total)

    #results = [item for sublist in results for item in sublist]

    #print('Results per sample at:', sout)
    #print('Results per variant at:', vout)

    return results


'''
my_ncores=mp.cpu_count()
start = timeit.default_timer()
load_file(my_ncores)
stop = timeit.default_timer()
print('Time to load the data by line (sec), numpy chunking, using ', my_ncores, 'cores: ', stop - start)
'''

def main():

    global WGS_file
    global IMPUTED_file
    global REF_file
    global disable_DS
    global max_total_rows
    global n_max
    global n_min
    global sout
    global vout
    global X_mode

    parser = argparse.ArgumentParser(usage='%(prog)s --ga <input_genotype_array.vcf.gz> --imputed <imputed_file.vcf.gz> --wgs <whole_genome_file.vcf.gz>\nUse -h or --help to display help.')
    parser.add_argument('--ga', dest='ga', type=str, required=False, default="", nargs=1, help='(optional for low pass) path to genotype array file in vcf.gz format, with tbi')
    parser.add_argument('--wgs', dest='wgs', type=str, required=True, nargs=1, help='path to whole genome file in vcf.gz format, with tbi')
    parser.add_argument('--bwgs', dest='bwgs', type=str, required=False, default="", nargs=1, help='path to whole genome file in bcf.gz format, with cbi index')
    parser.add_argument('--disable_DS', dest='disableDS', type=str, required=False, default="False", nargs=1, help='when DS is not present in the FORMAT')
    parser.add_argument('--bimputed', dest='bimputed', type=str, required=False, default="", nargs=1, help='path to imputed file in bcf.gz format, with cbi index')
    parser.add_argument('--imputed', dest='imputed', type=str, required=True, nargs=1, help='path to imputed file in vcf.gz format, with tbi')
    parser.add_argument('--ref', dest='ref', default="", type=str, required=False, nargs=1, help='optional, path to reference panel file in vcf.gz format, with tbi. Used for MAF calculation. WGS file will be used if no reference file is provided.')
    parser.add_argument('--max_total_rows', dest='max_total_rows', default=max_total_rows, type=int, required=False, nargs=1, help='maximun number of rows or variants to be loaded simultaneously, summing all chunks loaded by all cores')
    parser.add_argument('--max_per_core', dest='max_per_core', default=n_max, type=int, required=False, nargs=1, help='maximun number of variants per chunk per core, lower it to avoid RAM overload')
    parser.add_argument('--min_per_core', dest='min_per_core', default=n_min, type=int, required=False, nargs=1, help='minimun number of variants per chunk per core, increase to avoid interprocess communication overload')
    parser.add_argument('--sout', dest='sout', default="", type=str, required=False, nargs=1, help='optional output file path/name per sample, default is the same as the imputed file with _per_sample_results.txt suffix')
    parser.add_argument('--vout', dest='vout', default="", type=str, required=False, nargs=1, help='optional output file path/name per variant, default is the same as the imputed file with _per_variant_results.txt suffix')
    parser.add_argument('--xmode', dest='xmode', default="False", type=str, required=False, nargs=1, help='Option for developers, print additional scores.')

    try:
        args = parser.parse_args()
    except:
        #parser.print_help()
        sys.exit(0)

    #print(args)

    if(args.disableDS!="False"):
        disable_DS = True

    if(args.ga!=""):
        global GA_file
        GA_file = args.ga[0]
    else:
        global multimask_mode
        multimask_mode = True

    WGS_file = args.wgs[0]
    IMPUTED_file = args.imputed[0]
    
    if(args.bwgs!=""):
        global WGS_file_bcf
        WGS_file_bcf = args.bwgs[0]
    if(args.bimputed!=""):
        global IMPUTED_file_bcf
        IMPUTED_file_bcf = args.bimputed[0]
    if(args.max_total_rows!=max_total_rows):
        max_total_rows = args.max_total_rows[0]
    if(args.max_per_core!=n_max):
        n_max = args.max_per_core[0]
    if(args.min_per_core!=n_min):
        n_min = args.min_per_core[0]
    if(args.sout!=""):
        sout = args.sout[0]
    if(args.vout!=""):
        vout = args.vout[0]
    if(args.ref!=""):
        REF_file = args.ref[0]
    if(args.xmode!='False'):
        X_mode = args.xmode[0]

    #print(max_total_rows, n_max, n_min)

    start = timeit.default_timer()
    my_ncores=mp.cpu_count()
    results = load_file_chunks(my_ncores)
    stop = timeit.default_timer()
    #print('Total run time (sec):', stop - start)

#    for result in results:
#        print(result)


main()
