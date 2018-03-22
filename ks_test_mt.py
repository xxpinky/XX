# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 15:25:23 2017

@author: whtop
"""

import os
import sys
import glob

import numpy as np
import pandas as pd
from scipy import stats


def ks_2samp(data1, data2, alternative='less'):
    
    data1, data2 = map(np.asarray, (data1, data2))
    n1 = data1.shape[0]
    n2 = data2.shape[0]
    en = np.sqrt(n1 * n2 / float(n1 + n2))
    data1 = np.sort(data1)
    data2 = np.sort(data2)
    data_all = np.concatenate([data1,data2])
    cdf1 = np.searchsorted(data1,data_all,side='right')/(1.0*n1)
    cdf2 = (np.searchsorted(data2,data_all,side='right'))/(1.0*n2)
    
    if alternative == 'two-sided':
        d = np.max(np.absolute(cdf1 - cdf2))
        # Note: d absolute not signed distance
        try:
            prob = stats.distributions.kstwobign.sf((en + 0.12 + 0.11 / en) * d)
        except:
            prob = 1.0
    elif alternative == 'greater':
        d = np.max(cdf1 - cdf2)
    elif alternative == 'less':
        d = np.max(cdf2 - cdf1)
    if alternative in ['less', 'greater']:
        lambd = (en + 0.12 + 0.11 / en) * d
        prob = np.exp(-2 * lambd * lambd)
    
    return d, prob


def ks_test(dt_file, rank_file):
    
    drug_target = pd.read_csv(dt_file ,header=None, index_col=0)[1]
    gene_rank = pd.read_csv(rank_file, header=None)
    
    ks_p = []
    for d in drug_target.index.unique():
        d_t = drug_target[d]
        ps = []
        
        for i in gene_rank.columns:
            rank = gene_rank[i]
            if isinstance(d_t, str):
                sample_rank = rank[rank==d_t]
            elif isinstance(d_t, float):
                ps.append(0)
                continue         
            else:
                sample_rank = rank[rank.isin(d_t)]
            try:
                #p = stats.mstats.ks_2samp(rank.index, sample_rank.index, 'less')[1]
                p = ks_2samp(rank.index, sample_rank.index, 'less')[1]
                ps.append(p)
            except:
                ps.append(0)
        if sample_rank.shape[0]==0: print(d)
        ks_p.append(ps)
        
    w = open(os.path.basename(dt_file)[:-4] +'_out.csv', "w")
    for i in ks_p:
        line = []
        for j in i:
            line.append("{:.5g}".format(j))
        w.write(",".join(line) + "\n")
    w.close()
    
if __name__ == "__main__":
    
    
    name, directory, rank_file = sys.argv[1:]
    app, cli = glob.glob(os.path.join(directory, name+"*"))
    ks_test(app, rank_file)
    ks_test(cli, rank_file)
    

    
