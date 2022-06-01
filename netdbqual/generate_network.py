# Imports
from Bio import SeqIO
import gzip

from os import path
#sys.path.append('/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval')

import classify_acc as ca
import rec_parser as rec_parser
from importlib import reload  # Python 3.4+
from docopt import docopt

def main():
    doc =  """generate_network

    Usage: generate_network.py  [-i INPUT] [-f FORMAT]

    Options:
        -i INPUT   --input INPUT     input file [default: /Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/bacteria.187.genomic.gbff.gz]
        -f FORMAT   --format FORMAT    [default: gb]
    """

    args = docopt(doc, argv=None, help=True, version=None, options_first=False)
    
    print(args)


    f=args['-i']
    fmt=args['-f']
    if fmt not in ['gb', 'swiss']:
        raise RuntimeError("File format not understood: '{}'".format(fmt))

    data_dir = path.dirname(f)##'/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/'
    #f= data_dir+'uniprot_sprot.dat'
    # f=data_dir+'tmp_EC_3_4_11_4_partial_sequence.gpk'
    #f=data_dir+'tmp_EC_3_4_11_4_partial_sequence_append.gp'

    #gbs_iter=SeqIO.parse(gzip.open(f, "rt"), "gb")
    #gbs_iter=SeqIO.parse(f, "swiss")
    
    if f[-2:]=="gz":
        gbs_iter=SeqIO.parse(gzip.open(f, "rt", encoding='utf8', errors='ignore'), fmt)
    else:
        gbs_iter=SeqIO.parse(f, fmt)


    gb_file_no_ext=path.splitext(f)[0]#ath.basename(f).split(".")[0]
    #print(i)

    nodes=[]
    edges=[]
    node_filename=gb_file_no_ext+'_node.csv'
    pc_edge_filename=gb_file_no_ext+'_pc_edge.csv'
    xref_edge_filename=gb_file_no_ext+'_xref_edge.csv'

    with open(node_filename, "w") as node_file:
        with open(pc_edge_filename, "w") as pc_edge_file:       
            with open(xref_edge_filename, "w") as xref_edge_file:     
                for j,r in enumerate(gbs_iter):
                    # look at all CDS 
                    if j % 10000 == 0:
                        print(j)

                    (db,seq_type)=ca.classify_acc(r.id)
                    if seq_type == "unknown":
                        a=1
                    nodes,xref_strs, parent_child_edges,key=rec_parser.parseRecord(r, db, seq_type) 
                    node_ids=[n['id'] for n in nodes if n['id'][0:3]!='WP_']
                    if len(set(node_ids))!=len(list(node_ids)):
                        raise RuntimeError("Duplicate identifiers for multiple entities")

                    if xref_strs:
                        if(j==0):
                            print("Writing edges to {}".format(xref_edge_filename))
                            xref_edge_file.write("\t".join(['trg', 'trg_ver' 'src', 'src_ver', 'trg_type', 'src_type','edge_type'])+"\n")#"trg\tsrc\t\n")
                        xref_edge_file.write("\n".join(["\t".join(e) for e in xref_strs])+"\n")
                    if parent_child_edges:
                        if(j==0):
                            print("Writing edges to {}".format(pc_edge_filename))
                            pc_edge_file.write("\t".join(['trg', 'trg_ver' 'src', 'src_ver', 'trg_type', 'src_type','edge_type'])+"\n")#"trg\tsrc\t\n")
                        pc_edge_file.write("\n".join(["\t".join(e) for e in parent_child_edges])+"\n")
                    if nodes:
                        if(j==0):
                            print("Writing nodes to {}".format(node_filename))
                            node_file.write("\t".join(key)+"\n")
                        node_file.write("\n".join(["\t".join([str(x) for x in n.values()]) for n in nodes])+"\n")
                    #1
                    #if j > 10000:
                    #    break

   
if __name__ == "__main__":
    main()           