import os
import sys
import time

import numpy as np
import pandas as pd
import scipy as sp
import networkx as nx

from itertools import repeat
from collections import OrderedDict
from itertools import combinations
from multiprocessing import Pool


def read_ppi(f='/home/xxu/generank/BC/string_ppi_400.txt'):
    
    g = nx.Graph()
    for line in open(f):
        g.add_edges_from([[i.upper() for i in line.split()]])
    
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
    D1 = sp.sparse.csr_matrix(np.diag(1./degrees))
    A = np.eye(net.shape[0]) - d*np.dot(net.T, D1)
    b = (1-d)*norm_ex
    r = np.linalg.solve(A, b)
    
    return r


def run_rank(expr, ppi):
    
    start = time.time()
    name = os.path.basename(expr)[:-6]
    expr = pd.read_csv(expr, index_col=0,header =0)
    genes = OrderedDict({j:i \
                        for i, j in enumerate(expr.index)}) 
    try:
        conn = to_conn_matrix(ppi, genes)
        assert len(genes) == conn.shape[0]
    except Exception as exp:
        print("ERROR: {}: {}".format(name, exp))
        sys.stdout.flush()
        return
    
    conn = conn + conn.T
    conn = conn.tocsr()
    np.save(name, conn)
    print("{} ppi matrix:{}".format(name, repr(conn)))
    sys.stdout.flush()
    
    print('GeneRank process data:',name)
    try:
        rank = pd.DataFrame()
        for i in expr.columns:
            fold = expr[i].values
            rank_new = generank(conn, fold, 0.5)
            pro_rank = sorted(zip(rank_new, genes), reverse=True)       
            rank[i] = [i for _,i in pro_rank]
    except Exception as exp:
        print("ERROR: {}: {}".format(name, exp))
        sys.stdout.flush()
        return
        
    rank.to_csv("{}-rank.txt".format(name), header=False, index=False)
    end = time.time()
    print("{}: {}h elased".format(name, (end - start)/3600))
    sys.stdout.flush()
    
if __name__ == "__main__":
    
    ppi = read_ppi()
    dir_ = '/home/xxu/generank/54combin-data/'
    inputs = [os.path.join(dir_, i) for i in os.listdir(dir_)]
    

    p = Pool(5)
    p.starmap(run_rank, zip(inputs, repeat(ppi)))
    p.close()
    p.join()
        

        

