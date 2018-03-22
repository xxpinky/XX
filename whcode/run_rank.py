# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 17:49:32 2017

@author: whtop
"""


import time
import sys

import numpy as np
import pandas as pd
import scipy as sp

from collections import OrderedDict


def generank(net, fold, d):
    """ 
    input data: 
    net: connectivity  matrix (zero/one, symmetric with zero diag)
    fold: vector of expression levels (non-negative)                
    d: parameter in algorithm
    
    output is   
    r: vector of rankings 
    """
    fold = abs(fold)
    norm_ex = fold/max(fold)
    degrees = net.sum(axis=1)
    degrees[degrees==0] = 1
    degrees = np.array(degrees)[:, 0]
    D1 = sp.sparse.csr_matrix(np.diag(1./degrees))
    A = np.eye(net.shape[0]) - d*np.dot(net.T, D1)
    b = (1-d)*norm_ex
    r = np.linalg.solve(A, b)
    
    return r

if __name__ == "__main__":
    
    expr_file, ppi_file,rank_file = sys.argv[1:]
    start = time.time()
    expr = pd.read_csv(expr_file, index_col=0)
    genes = OrderedDict({j.strip().upper():i \
                        for i, j in enumerate(expr.index)})
    conn = np.load(ppi_file).item()
    assert len(genes) == conn.shape[0]
    
    print("running gene rank...")
    rank = pd.DataFrame()
    for i in expr.columns:
        print("processing sample: {}".format(i))
        fold = expr[i].values
        rank_new = generank(conn, fold, 0.5)
        pro_rank = sorted(zip(rank_new, genes), reverse=True)       
        rank[i] = [i for _,i in pro_rank]
        
    rank.to_csv(rank_file, header=False, index=False)
    end = time.time()
    print("{}h elased".format((end - start)/3600))
        

        
