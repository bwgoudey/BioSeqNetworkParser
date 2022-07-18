# Imports
from Bio import SeqIO
import gzip

from os import path
#sys.path.append('/Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval')

import classify_acc as ca
import rec_parser as rec_parser
from importlib import reload  # Python 3.4+
from docopt import docopt


#-i INPUT   --input INPUT     input file [default: /Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/uniprot_sprot.dat.bgz]
def main():
    doc =  """generate_network

    Usage: generate_network.py  [-i INPUT] [-f FORMAT]

    Options:
        -i INPUT   --input INPUT     input file [default: /Users/bgoudey/Research/BioDbNetQual/BioDbPropagationEval/data/gbbct21.seq.gz]
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

    #nodes=[]
    #edges=[]
    node_filename=gb_file_no_ext+'_node.csv'
    pc_edge_filename=gb_file_no_ext+'_pc_edge.csv'
    xref_edge_filename=gb_file_no_ext+'_xref_edge.csv'
    func_edge_filename=gb_file_no_ext+'_func_edge.csv'
    annot_edge_filename=gb_file_no_ext+'_annot_edge.csv'   
    with open(node_filename, "w") as node_file, \
         open(pc_edge_filename, "w") as pc_edge_file, \
         open(xref_edge_filename, "w") as xref_edge_file, \
         open(annot_edge_filename, "w") as annot_edge_file, \
         open(func_edge_filename, "w") as func_edge_file:     
            for j,r in enumerate(gbs_iter):
                # look at all CDS 
                if j % 10000 == 0:
                    print(j)

                (db,seq_type)=ca.classify_acc(r.id)
                if seq_type == "unknown":
                    a=1
                #nodes,xref_strs, parent_child_edges,key,func_edges=rec_parser.parseRecord(r, db, seq_type) 
                g=rec_parser.parseRecord(r, db, seq_type) 

                node_ids=[n['id'] for n in g['nodes'] if n['id'][0:3]!='WP_' and n['id']]                
                if len(set(node_ids))!=len(list(node_ids)):
                    raise RuntimeError("Duplicate identifiers for multiple entities")

                if g['xref_strs']:
                    if(j==0):
                        print("Writing edges to {}".format(xref_edge_filename))
                        xref_edge_file.write("\t".join(['trg', 'trg_ver', 'src', 'src_ver', 'trg_type', 'src_type','edge_type'])+"\n")#"trg\tsrc\t\n")
                    xref_edge_file.write("\n".join(["\t".join(e) for e in g['xref_strs']])+"\n")

                if g['func_edges']:
                    if(j==0):
                        print("Writing edges to {}".format(func_edge_filename))
                        func_edge_file.write("\t".join(['id', 'annot', 'eco', 'source'])+"\n")#"trg\tsrc\t\n")
                    func_edge_file.write("\n".join(["\t".join(e.values()) for e in g['func_edges']])+"\n")

                if g['parent_child_edges']:
                    if(j==0):
                        print("Writing edges to {}".format(pc_edge_filename))
                        cols=['trg', 'trg_ver', 'src', 'src_ver', 'trg_type', 'src_type','edge_type']
                        pc_edge_file.write("\t".join(cols)+"\n")#"trg\tsrc\t\n")
                    pc_edge_file.write("\n".join(["\t".join(e) for e in g['parent_child_edges']])+"\n")

                if g['nodes']:
                    if(j==0):
                        print("Writing nodes to {}".format(node_filename))
                        node_file.write("\t".join(g['key'])+"\n")
                    node_file.write("\n".join(["\t".join([str(x) for x in n.values()]) for n in g['nodes']])+"\n")

                #HACK
                g['annot_edges']=list(filter(lambda x: x is not None, g['annot_edges']))
                if g['annot_edges']:
                    cols=['trg_id','trg_seq_ver', 'acc', 'seq_version', 'identity', 'db', 'type', 'model']
                    col_names=['trg_id','trg_seq_ver', 'src_acc', 'src_seq_ver', 'identity', 'db', 'type', 'model']
                    if(j==0):
                        print("Writing edges to {}".format(annot_edge_filename))
                        annot_edge_file.write("\t".join(col_names)+"\n")#"trg\tsrc\t\n")
                    annot_edge_file.write("\n".join(["\t".join([str(e[k]) for k in cols]) for e in g['annot_edges'] if e is not None])+"\n")
   
if __name__ == "__main__":
    main()           