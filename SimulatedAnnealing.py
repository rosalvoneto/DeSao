import numpy as np
from math import exp
from LibDeepFMPO.global_parameters import MAX_SWAP, TIMES
from rewards2 import evaluate_mol, modify_fragment

from LibDeepFMPO.rewards import decode
from rdkit import Chem
import json
import time
from save_results import save_results

def coolTemperature(To, k, alpha):
    return To/(1 + alpha*k)

def getAcceptanceProbability(delta, t):
    return exp(-delta/t)

def search(X, decodings, MaxIterations, temperatura_inicial, alpha, S, s_temp, SampleSize):

    start_time = time.time()

    mols       = {}
    dic_temp   = {}
    dic_delta  = {}
    dic_prob   = {}
    dic_accept = {}
    dic_news   = {}
    dic_news_ok = {}
    s_count = 0
    
    print('*'*45)
    print(f'Max. iterarion: {MaxIterations}')
    print(f'S: {S}')
    print(f's_temp: {s_temp}')
    print(f'temperatura_inicial: {temperatura_inicial}')
    print(f'alpha: {alpha}')
    print(f'SampleSize: {SampleSize}')
    print('*'*45)

    # Build population
    rand_n = np.random.randint(0,X.shape[0],SampleSize)
    LL = X[rand_n].copy()   

    temp = temperatura_inicial
    
    with open('Results/New_Mols.txt', 'w') as arquivo:
    # For every iteration
        for e in range(MaxIterations):

            dic_temp[e] = temp            

            # For all modification steps
            for t in range(TIMES):

                # Select actions
                for i in range(SampleSize):

                    a = np.random.randint(0, 62)

                    a = int(a // MAX_SWAP)

                    if a == 12:
                        a = 11

                    s = a % MAX_SWAP                                        
                    mol_orriginal = LL[i,a].copy()
                    mol_orriginal_av = LL[i].copy()
                    LL[i,a] = modify_fragment(LL[i,a], s)                    
                    fr, score, isNew = evaluate_mol(LL[i], e, decodings)                    
                    if not(isNew):
                        s_count +=1
                        continue
                    s_count = 0
                    dic_news[e] = dic_news.get(e, 0) + 1                  
                    if all(fr):
                        mol_new = decode(LL[i], decodings)
                        smiles_code = Chem.MolToSmiles(mol_new)
                        new = mols.get(smiles_code, 0)
                        if new == 0:
                            mols[smiles_code] =  1
                            arquivo.write(f'{smiles_code}\n')
                            dic_news_ok[e] = dic_news_ok.get(e, 0) + 1
                        else:
                            mols[smiles_code] =  new + 1                            
                    else:
                        _, score_old,_ = evaluate_mol(mol_orriginal_av, e, decodings)                    
                        delta = (score_old - score)
                        dic_delta[e] = delta
                        if (delta > 0):
                            acceptanceProbability = getAcceptanceProbability(delta, temp)                            
                            dic_prob[e] = acceptanceProbability
                            if (np.random.rand()<=acceptanceProbability):
                                LL[i,a] = mol_orriginal
                                dic_accept[e] = 1                                    
                    if (temp < s_temp):
                        break
            print (f"Iteration {e+1}")
            temp = coolTemperature(temperatura_inicial, e+1, alpha)
            if (s_count >= S):
                if (temp < s_temp):
                    print(f'Stop, because {S} iteration without new molecules!')
                    break
    end_time = time.time()
    duration = (end_time - start_time) /60

    print(f"Execution time: {duration:.2f} minutes")
    print(f"Total number of new molecules that meet the criteria: {sum(dic_news_ok.values())}")

    save_results(
        dic_temp=dic_temp,
        dic_delta=dic_delta,
        dic_prob=dic_prob,
        dic_accept=dic_accept,
        dic_news_ok=dic_news_ok,
        dic_news=dic_news,
        results_dir="Results"
    )

    return True
