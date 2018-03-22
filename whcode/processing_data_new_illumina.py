# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 20:21:12 2017

@author: whtop
"""

import os
import argparse

import pandas as pd
import numpy as np

from bs4 import BeautifulSoup

def data_integration(folder):
    
    # dict that store the corresponding genes for each probe
    probe_gene = {}
    # stores expression data for each sample
    expr = {}
    # unique genes in the expression data
    genes = set()
    pos_header = []
    neg_header = []
    # remove genes that do not overlap in different series
    remove_list = set()
    sample_files = []
    
    for file in  os.listdir(folder):
    
        if file.endswith('xml'):                      
            # first get the sample imformation
            sample_info_file = file 
                
        elif file.startswith('GPL'):          
            # next get the probes and genes
            probe_info_file = file
                
        else:
            # final the gene expression for each sample
            sample_files.append(file)
        
    # process sample_info_file
    soup = BeautifulSoup(open(os.path.join(folder,sample_info_file),'rb'), 'xml')
    samples = soup.find_all('Sample')
    
    for sample in samples:
        # characteristics = sample.find_all("Characteristics") 
        if " ".join(args.positive) in sample.Title.text:
            pos_header.append(sample['iid'])

        if " ".join(args.negtive) in sample.Title.text:
            neg_header.append(sample['iid'])

    # process probe_info_file
    for i in open(os.path.join(folder,probe_info_file)):
        items = i.split(sep='\t')
        p = items[0]
        if items[10]:
            gene_list = [i.strip() for i in items[10].split(sep='///')]
            probe_gene[p]=gene_list
                         
            genes.update(gene_list)
    
    genes = list(genes)
        
    # process sample files
    for file in sample_files:                 
        sample = file.split(sep='-')[0]
        if sample in neg_header or sample in pos_header:
            print("processing {}......".format(sample))
            expr[sample] = []
            single_sample_expr = {}
            
            for i in open(os.path.join(folder,file)): 
                probe, expression = i.strip().split()[:2]
                #if float(expression) < 0: expression=0
                
                if probe_gene.get(probe,0):
                    for g in probe_gene[probe]:
                        if g in single_sample_expr:
                            single_sample_expr[g].append(float(expression))
                        else:
                            single_sample_expr[g] = [float(expression)]
    
    
            c = -1
            for g in genes:
                c += 1
                if g not in single_sample_expr:
                    #print ("gene {} not in sample {}".format(g, sample))
                    remove_list.add(c)
                    expr[sample].append(0)
                    continue
                expr[sample].append(np.mean(single_sample_expr[g]))
                
    return genes, expr, pos_header, neg_header, remove_list

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser('preprocessing the expression profile and output the result')
    
    parser.add_argument('folder',
                        help='data folder for expression data')
    parser.add_argument('-o', '--output', default='output',
                        help='output prefix for processed case and control data')
    parser.add_argument('-i', '--interaction', default='ppi.txt',
                        help='output name of protein ppi')
    # Note that both key words for positive samples and negative samples are needed
    parser.add_argument('-p', '--positive', 
                        help='keywords to select positive samples', nargs='+')
    parser.add_argument('-n', '--negtive', default=None, 
                        help='keywords to select negative samples', nargs='+')
    args = parser.parse_args()
    
    genes, expr, pos_header, neg_header, remove_list = data_integration(args.folder)
    
    # first remove genes that do not overlap
    for sample in expr:
        for index in sorted(remove_list, reverse=True):
            expr[sample].pop(index)
            
    for index in sorted(remove_list, reverse=True):
        genes.pop(index)
                  
    # reorganize files into the input format required by ngf   
    
    pos = pd.DataFrame(index=genes, columns=pos_header)
    neg = pd.DataFrame(index=genes, columns=neg_header)
    
    
    for sample in neg_header:
        neg[sample] = expr[sample]

    for sample in pos_header:
        pos[sample] = expr[sample]
    
    unique_genes = set()
    wh = open(args.interaction, 'w')
    with open('string_ppi_gene_id.txt') as f:
        for i in f:
            j, k, s= i.split()
            if int(s)<770: continue
            if j in genes:
                if k in genes:
                    wh.write(j+'\t'+k+'\n')
                    unique_genes.update([j,k])
    wh.close()
    
    for g in genes:
        if g not in unique_genes:
            neg.drop(g, inplace=True)
            pos.drop(g, inplace=True)
    
    pos.to_csv(args.output + '.case', index_label='#gene', float_format='%.7f', sep='\t')
    neg.to_csv(args.output + '.control', index_label='#gene', float_format='%.7f', sep='\t')
