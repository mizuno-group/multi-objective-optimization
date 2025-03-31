import os
import pandas as pd
from tqdm import tqdm
from pathlib import Path
import pubchempy as pcp
import time
import random
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
from scipy.stats import norm
from sklearn.preprocessing import robust_scale

from util import set_logger, getLogger, file_checker, robust_z

def prep_for_09(test, lig="ago"):
    if test != "07":
        output_path = f"../../data/processed/09{test}_{lig}/cas_sev.tsv"
    else:
        output_path = f"../../data/processed/09{test}/cas_sev.tsv"
    overwrite = False
    if file_checker(output_path, overwrite):
        df_AR_ago = pd.read_csv("../../data/processed/tox21/tox21_AR_ago.tsv", sep="\t", header=None)
        df_AR_anta = pd.read_csv("../../data/processed/tox21/tox21_AR_anta.tsv", sep="\t", header=None)
        df_ERa_ago = pd.read_csv("../../data/processed/tox21/tox21_ERa_ago.tsv", sep="\t", header=None)
        df_ERa_anta = pd.read_csv("../../data/processed/tox21/tox21_ERa_anta.tsv", sep="\t", header=None)
        df_ERb_ago = pd.read_csv("../../data/processed/tox21/tox21_ERb_ago.tsv", sep="\t", header=None)
        df_ERb_anta = pd.read_csv("../../data/processed/tox21/tox21_ERb_anta.tsv", sep="\t", header=None)

        if test == "01" and lig == "ago":
            data = [df_ERa_ago, df_ERb_ago]
        elif test == "01" and lig == "anta":
            data = [df_ERa_anta, df_ERb_anta]
        elif test == "02" and lig == "ago":
            data = [df_ERa_ago]
        elif test == "02" and lig == "anta":
            data = [df_ERa_anta]
        elif test == "04" and lig == "ago":
            data = [df_AR_ago]
        elif test == "04" and lig == "anta":
            data = [df_AR_anta]
        elif test == "05" and lig == "ago":
            data = [df_AR_ago]
        elif test == "05" and lig == "anta":
            data = [df_AR_anta]
        elif test == "07":
            data = [df_ERa_ago, df_ERa_anta]

        cas_sev = dict()
        for x in data:
            df = x
            for i in range(len(df)):
                if df.iloc[i,0] not in cas_sev.keys():
                    cas_sev[df.iloc[i,0]] = set()
                    cas_sev[df.iloc[i,0]].add(df.iloc[i,1])
                else:
                    cas_sev[df.iloc[i,0]].add(df.iloc[i,1])  
        
        cas_sev_use = dict()
        for i in cas_sev.keys():
            if list(cas_sev[i]) == [0, 1]:
                cas_sev_use[i] = 1
            else:
                cas_sev_use[i] = list(cas_sev[i])[0]

        cas = dict()
        tsv = []
        for x in data:
            df = x
            for i in range(len(df)):
                if df.iloc[i,0] not in cas.keys():
                    cas[df.iloc[i,0]] = 0
                    col = []
                    for n in range(len(df.iloc[i])):
                        if n == 1:
                            col.append(cas_sev_use[df.iloc[i,0]])
                            continue
                        col.append(df.iloc[i,n])
                    tsv.append(col)
                else:
                    continue
        if test != "07":
            os.makedirs(f"../../data/processed/09{test}_{lig}", exist_ok=True)
        else:
            os.makedirs(f"../../data/processed/09{test}", exist_ok=True)
        if test != "07":
            pd.DataFrame(tsv).to_csv(f"../../data/processed/09{test}_{lig}/cas_sev.tsv", sep="\t", header=None, index=False)
        else:
            pd.DataFrame(tsv).to_csv(f"../../data/processed/09{test}/cas_sev.tsv", sep="\t", header=None, index=False)

def prep_for_ga(test_num, lig="ago"):
    if "09" in test_num and "07" not in test_num:
        output_path = f"../../data/processed/{test_num}_{lig}/for_GA.tsv"
    else:
        output_path = f"../../data/processed/{test_num}/for_GA.tsv"

    overwrite = False
    if file_checker(output_path, overwrite):
        if "09" in test_num and "07" not in test_num:
            whole_df = pd.read_csv(f"../../data/processed/{test_num}_{lig}/cas_sev.tsv", sep="\t", header=None)
            val_df = pd.read_csv(f"../../data/raw/validation_test_cas/{test_num}_{lig}.csv", header=None)

        else:
            whole_df = pd.read_csv(f"../../data/processed/{test_num}/cas_sev.tsv", sep="\t", header=None)
            val_df = pd.read_csv(f"../../data/raw/validation_test_cas/{test_num}.csv", header=None)
        
        cas_val = set(val_df.iloc[:,0])

        tsv = []
        for i in range(len(whole_df)):
            if whole_df.iloc[i,0] in cas_val:
                continue
            col = []
            for n in range(len(whole_df.iloc[i])):
                col.append(whole_df.iloc[i,n])
            tsv.append(col)

        if "09" in test_num and "07" not in test_num:
            pd.DataFrame(tsv).to_csv(f"../../data/processed/{test_num}_{lig}/for_GA.tsv", sep="\t", header=None, index=False)
        else:
            pd.DataFrame(tsv).to_csv(f"../../data/processed/{test_num}/for_GA.tsv", sep="\t", header=None, index=False)

def prep_validation(test_num, lig="ago"):
    if "09" in test_num and "07" not in test_num:
        output_path = f"../../data/processed/{test_num}_{lig}/validation.tsv"
    else:
        output_path = f"../../data/processed/{test_num}/validation.tsv" 

    overwrite = False
    if file_checker(output_path, overwrite):
        if "09" in test_num and "07" not in test_num:
            whole_df = pd.read_csv(f"../../data/processed/{test_num}_{lig}/cas_sev.tsv", sep="\t", header=None)
            val_df = pd.read_csv(f"../../data/raw/validation_test_cas/{test_num}_{lig}.csv", header=None)

        else:
            whole_df = pd.read_csv(f"../../data/processed/{test_num}/cas_sev.tsv", sep="\t", header=None)
            val_df = pd.read_csv(f"../../data/raw/validation_test_cas/{test_num}.csv", header=None)

        cas_val = set(val_df.iloc[:,0])

        tsv = []
        for i in range(len(whole_df)):
            if whole_df.iloc[i,0] not in cas_val:
                continue
            col = []
            for n in range(len(whole_df.iloc[i])):
                col.append(whole_df.iloc[i,n])
            tsv.append(col)
        
        df_tsv = pd.DataFrame(tsv)
        if "09" in test_num and "07" not in test_num:
            pd.DataFrame(tsv).to_csv(f"../../data/processed/{test_num}_{lig}/validation.tsv", sep="\t", header=None, index=False)
        else:
            pd.DataFrame(tsv).to_csv(f"../../data/processed/{test_num}/validation.tsv", sep="\t", header=None, index=False)
        
        cas_in_whole = set(df_tsv.iloc[:,0])
        cas_not_in_whole = list(cas_val - cas_in_whole)
        print(cas_not_in_whole)

        tsv_2 = []
        for cas in cas_not_in_whole:
            tsv_2.append([cas, ""])
            
        tsv_2 = sorted(tsv_2, key=lambda x: x[0])

        if "09" in test_num and "07" not in test_num:
            pd.DataFrame(tsv_2).to_csv(f"../../data/processed/{test_num}_{lig}/cas_not_in_whole.tsv", sep="\t", header=None, index=False)
        else:
            pd.DataFrame(tsv_2).to_csv(f"../../data/processed/{test_num}/cas_not_in_whole.tsv", sep="\t", header=None, index=False)

def prep_pubchem_bycas(cas_list, property):
    set_logger()
    logger = getLogger(__name__)
    
    result_df = pd.DataFrame()
    error = []

    for cas in tqdm(cas_list):
        try:
            results = pcp.get_properties(property, cas, 'name', as_dataframe=True)
            results["CAS"] = cas
            result_df = pd.concat([results, result_df],axis=0,join='outer',sort=True)
        except:
            cas_1 = "CAS-" + cas
            logger.info("cas = %s", cas_1)
            try:
                results = pcp.get_properties(property, cas, 'name', as_dataframe=True)
                results["CAS"] = cas
                result_df = pd.concat([result, result_df],axis=0,join='outer',sort=True)
            except Exception as e:
                error.append(cas)
                logger.error("An error occurred: %s", e)
                continue

        time.sleep(random.uniform(1,3))
    
    result_df = result_df.reset_index()
    return result_df, error

def for_lookup(data):
    # get ecfp
    smiles = data.iloc[:,2]
    ecfp = []
    radius = 2
    nBits = 2048
    for i in range(len(smiles)):
        mol = Chem.MolFromSmiles(smiles[i])
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
        ecfp.append(fp)
    
    # scaling phisical property
    logp = list(data.iloc[:,3])
    tpsa = list(data.iloc[:,4])
    logp_z = robust_z(logp) 
    tpsa_z = robust_z(tpsa)
    
    phisi = []
    for i in tqdm(range(len(logp_z))):
        for n in range(len(logp_z)):
            if i >= n:
                continue
            v_i = np.array([logp_z[i], tpsa_z[i]])
            v_n = np.array([logp_z[n], tpsa_z[n]])
            dist =  np.linalg.norm(v_i - v_n)
            phisi.append(dist)
    
    max_phisi = max(phisi)

    phisi_m = []
    for i in range(len(phisi)):
        phisi_m.append(phisi[i] / max_phisi)

    cas = data.iloc[:,1]
    
    num = 0
    output = dict()
    for i in tqdm(range(len(data))):
        for n in range(len(data)):
            if i >= n:
                continue
            output[f"{cas[i]}, {cas[n]}"] = []
            # structure
            output[f"{cas[i]}, {cas[n]}"].append(DataStructs.TanimotoSimilarity(ecfp[i],ecfp[n]))
            # physical property
            dist = phisi_m[num]
            output[f"{cas[i]}, {cas[n]}"].append(dist)
            num += 1
    
    return output

def for_lookup_val(data):
    # get ecfp
    smiles = data.iloc[:,2]
    ecfp = []
    radius = 2
    nBits = 2048
    for i in range(len(smiles)):
        mol = Chem.MolFromSmiles(smiles[i])
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
        ecfp.append(fp)
    
    # scaling phisical property
    logp = list(data.iloc[:,4])
    tpsa = list(data.iloc[:,4])
    logp_z = robust_z(logp) 
    tpsa_z = robust_z(tpsa)
    
    phisi = []
    for i in tqdm(range(len(logp_z))):
        for n in range(len(logp_z)):
            if i >= n:
                continue
            v_i = np.array([logp_z[i], tpsa_z[i]])
            v_n = np.array([logp_z[n], tpsa_z[n]])
            dist =  np.linalg.norm(v_i - v_n)
            phisi.append(dist)
    
    max_phisi = max(phisi)

    phisi_m = []
    for i in range(len(phisi)):
        phisi_m.append(phisi[i] / max_phisi)

    cas = data.iloc[:,0]
    
    num = 0
    output = dict()
    for i in tqdm(range(len(data))):
        for n in range(len(data)):
            if i >= n:
                continue
            output[f"{cas[i]}, {cas[n]}"] = []
            # structure
            output[f"{cas[i]}, {cas[n]}"].append(DataStructs.TanimotoSimilarity(ecfp[i],ecfp[n]))
            # physical property
            dist = phisi_m[num]
            output[f"{cas[i]}, {cas[n]}"].append(dist)
            num += 1
    
    return output