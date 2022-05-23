# Imports
from Bio import SeqIO
import gzip
import progressbar
import glob
import re 

from os import path
import sys
sys.path.append('./')

import netdbqual.classify_acc as ca
import netdbqual.parser as parser
from importlib import reload  # Python 3.4+
data_dir = './data/'


f= data_dir+'uniprot_sprot.dat'
    
#gbs_iter=SeqIO.parse(gzip.open(f, "rt"), "gb")
gbs_iter=SeqIO.parse(f, "swiss")
gb_file_no_ext=path.basename(f).split(".")[0]
#print(i)

nodes=[]
edges=[]
with open(data_dir + gb_file_no_ext+'_node.csv', "w") as node_file:
    with open(data_dir + gb_file_no_ext+'_edge.csv', "w") as edge_file:            
        for j,r in enumerate(gbs_iter):
            # look at all CDS 
            print(j)
            (db,seq_type)=ca.classify_acc(r.id)
            nodes,edges=parser.parseRecord(r, db, seq_type) 
            print(j)

            if edges:
                edge_file.write('\n'.join(edges)+"\n")
            if nodes:
                    node_file.write('\n'.join(nodes)+"\n")
            #1
            