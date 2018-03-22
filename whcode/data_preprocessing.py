# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:17:44 2017

@author: whtop
"""

import os
import argparse


import mygene
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
    
    for file in  os.listdir(folder):
    
        if file.endswith('xml'):
            # first get the sample imformation
            soup = BeautifulSoup(open(os.path.join(folder,file),'rb'), 'xml')
            samples = soup.find_all('Sample')
            
            for sample in samples:
                # characteristics = sample.find_all("Characteristics") 
                if " ".join(args.positive) in sample.Title.text:
                    pos_header.append(sample['iid'])
                if " ".join(args.negtive) in sample.Title.text:
                    neg_header.append(sample['iid'])
                
        elif file.startswith('GPL'):
            # next get the probes and genes
            for i in open(os.path.join(folder,file)):
                items = i.split(sep='\t')
                p = items[0]
                probe_gene[p] = set()
                if items[11]:
                    gene_list = [i.strip() for i in items[11].split(sep='///')]
                    probe_gene[p].update(gene_list)
                                 
                    genes.update(gene_list)
            
            genes = list(genes)
                
        else:
            # final the gene expression for each sample
            sample = file.split(sep='-')[0]
            expr[sample] = []
            single_sample_expr = {}
            
            for i in open(os.path.join(folder,file)): 
                probe, expression = i.strip().split()[:2]
                
                if probe_gene.get(probe,0):
                    for g in probe_gene[probe]:
                        if g in single_sample_expr:
                            single_sample_expr[g].append(float(expression))
                        else:
                            single_sample_expr[g] = [float(expression)]
    
            for g in genes:    
                expr[sample].append(np.mean(single_sample_expr[g]))

    return probe_gene, genes, expr, pos_header, neg_header


def report(ensemble_id, ensemble={}):
    if isinstance(ensembl_id, list):
        for p_id in ensembl_id:
            if p_id not in ensembl:
                ensembl[p_id] = entrez
            elif ensembl[p_id] != entrez:
                print("{} has two different ids: {} {}".format(p_id, ensembl[p_id], entrez))
    else:
        if ensembl_id not in ensembl:
            ensembl[ensembl_id] = entrez
        elif ensembl[ensembl_id] != entrez:
            print("{} has two different ids: {} {}".format(ensembl_id, ensembl[ensembl_id], entrez)) 
            


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
    
    probe_gene, genes, expr, pos_header, neg_header = data_integration(args.folder)
    query = mygene.MyGeneInfo()
    
    # reorganize files into the input format required by ngf   
    
    pos = pd.DataFrame(index=genes, columns=pos_header)
    neg = pd.DataFrame(index=genes, columns=neg_header)
    
    
    df = query.getgenes(genes, fields='ensembl.protein', as_dataframe=True)
    ensembl = {}

    for entrez, proteins in df.ensembl.dropna().items():
        if isinstance(proteins, list):
            for protein in proteins:
                ensembl_id = protein['protein']
                report(ensembl_id, ensembl)
        
        else:
            ensembl_id = proteins['protein']
            report(ensembl_id, ensembl)
    
    for sample in neg_header:
        neg[sample] = expr[sample]
    
    for sample in pos_header:
        pos[sample] = expr[sample]
    
    unique_genes = set()
    wh = open(args.interaction, 'w')
    with open('human_string.txt') as f:
        next(f)
        for i in f:
            j, k, s= i.replace("9606.", "").strip().split()
            if int(s)<770: continue
            if j in ensembl:
                if k in ensembl:
                    wh.write(j+'\t'+k+'\n')
                    unique_genes.update([j,k])
    wh.close()
    
    pos_e = pd.DataFrame(columns=pos_header)
    neg_e = pd.DataFrame(columns=neg_header)
    
    for ensembl_id in ensembl:
        if ensembl_id in unique_genes:
            n = neg.loc[ensembl[ensembl_id]]
            p = pos.loc[ensembl[ensembl_id]]
            n.name = ensembl_id
            p.name = ensembl_id
            neg_e = neg_e.append(n)
            pos_e = pos_e.append(p)    
                
    neg_e.to_csv(args.output + '.case', index_label='#gene', float_format='%.7f', sep='\t')
    pos_e.to_csv(args.output + '.control', index_label='#gene', float_format='%.7f', sep='\t')
