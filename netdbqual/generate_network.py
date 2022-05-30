# Imports
from Bio import SeqIO
import gzip
import progressbar
import glob
import re 

from os import path
import sys
sys.path.append('/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval')

import netdbqual.classify_acc as ca
import netdbqual.parser as parser
from importlib import reload  # Python 3.4+
data_dir = '/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/'


f= data_dir+'uniprot_sprot.dat'
#f=data_dir+'tmp_EC_3_4_11_4_partial_sequence.gpk'
f=data_dir+'tmp_EC_3_4_11_4_partial_sequence_append.gp'

#gbs_iter=SeqIO.parse(gzip.open(f, "rt"), "gb")
#gbs_iter=SeqIO.parse(f, "swiss")
gbs_iter=SeqIO.parse(f, "gb")
gb_file_no_ext=path.basename(f).split(".")[0]
#print(i)

nodes=[]
edges=[]
with open(data_dir + gb_file_no_ext+'_node.csv', "w") as node_file:
    with open(data_dir + gb_file_no_ext+'_edge.csv', "w") as edge_file:            
        for j,r in enumerate(gbs_iter):
            # look at all CDS 
            if j % 10000 == 0:
                print(j)

            (db,seq_type)=ca.classify_acc(r.id)
            nodes,edges,key=parser.parseRecord(r, db, seq_type) 

            if edges:
                if(j==0):
                    edge_file.write("trg\tsrc\n")
                edge_file.write("\n".join(["\t".join(e) for e in edges])+"\n")
            if nodes:
                if(j==0):
                    node_file.write("\t".join(key)+"\n")
                node_file.write("\n".join(["\t".join([str(x) for x in n]) for n in nodes])+"\n")
            #1
            #if j > 10000:
            #    break