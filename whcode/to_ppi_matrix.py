# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 17:46:42 2017

@author: whtop
"""


import sys

import numpy as np
import scipy as sp
import networkx as nx

from collections import OrderedDict
from itertools import combinations


def read_ppi(f='/home/xxu/generank/BC/string_ppi_400.txt'):
    
    g = nx.Graph()
    for line in open(f):
        g.add_edges_from([[i.upper() for i in line.split()]])
    
    return g


def read_gene(f='expr.csv'):
    
    f = open(f)
    next(f)
    return OrderedDict({j.split(',')[0].upper():i \
                        for i, j in enumerate(f)})


def to_conn_matrix(ppi, genes):
    
    n = len(genes)
    c=0
    conn = sp.sparse.lil_matrix((n, n))
    for i,j in combinations(genes, 2):
        c+=1
        if ppi.has_edge(i, j):
            conn[genes[i], genes[j]]=1
    
    return conn

if __name__ == "__main__":
    
    expr_file, out_name = sys.argv[1:]
    genes = read_gene(expr_file)
    ppi = read_ppi()
    conn = to_conn_matrix(ppi, genes)
    conn = conn + conn.T
    conn = conn.tocsr()
    print(repr(conn))
    np.save(out_name, conn)
