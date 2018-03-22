# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 12:36:00 2017
Modified generank to use sparse matrix.
Save memory
@author: whtop
"""


import time

import numpy as np
import pandas as pd
import scipy as sp
import networkx as nx

from itertools import combinations
from config import Config
from collections import OrderedDict


def read_ppi():
    
    g = nx.Graph()
    ppi = open(Config.PPI)
    for line in ppi:
        g.add_edges_from([[i.upper() for i in line.split()[:2]]])
    ppi.close()
    
    return g


def to_conn_matrix(ppi, genes):
    
    n = len(genes)
    c=0
    conn = sp.sparse.lil_matrix((n, n))
    for i,j in combinations(genes, 2):
        c+=1
        if ppi.has_edge(i, j):
            conn[genes[i], genes[j]]=1
    
    return conn
    
        
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
    D1 = sp.sparse.diags(1./degrees, format='csr')
    A = sp.sparse.identity(net.shape[0]) - d*np.dot(net.T, D1)
    b = (1-d)*norm_ex
    r = sp.sparse.linalg.spsolve(A, b)
    
    return r

if __name__ == "__main__":
    
    start = time.time()
    
    print("processing data...")
    expr = pd.read_csv("PSO-com-fd.csv", index_col=0, header=None)[1]
    genes = OrderedDict([(j.strip().upper(),i) \
                        for i, j in enumerate(expr.index)]) 
    ppi = read_ppi()
    conn = to_conn_matrix(ppi, genes)
    conn = conn + conn.T
    
    print("running gene rank...")
    rank_new = generank(conn, expr.values, 0.5)
    pro_rank = [i[1] for i in sorted(zip(rank_new, genes), reverse=True)]       

    end = time.time()
    print("{}h elased".format((end - start)/3600))
