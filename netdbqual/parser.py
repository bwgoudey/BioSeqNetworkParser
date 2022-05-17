from __future__ import annotations
import netdbqual.utils as u
import re 



class RecNode: 
    def __init__(self, r, db, seq_type):
        db=db
        seq_type=seq_type
        edges=[]
            



def extract_nodes(r):
    db,seq_type=u.classify_acc(p.id)[0]
    rec_type=["rec", "subrec"][hasattr(r, "qualifiers")]
    
    #Extract the edges and node from the given record
    nodes=[]
    edges=[]
    node,edge=RecNode(r, db, seq_type).getNodeEdgesStr()
    nodes += node
    edges += edge

    # If we have a nucleotide record, extract nodes/edges from all CDS as well
    if rec_type == "rec" and seq_type=="nuc":
        for p in r.features:
                node,edge=parse_gb_protein(p)
                nodes += node
                edges += edge
    
    return(nodes, edges)




def parse_gb_protein_gff(p):
    db=u.classify_acc(p.id)[0]
    if 'molecule_type' not in p.annotations:
        return ([],[])

    prod=p.description.split(' [')[0] 
    if db == "uniprot":
        prod=p.description.split(";")[0].split("=")[-1]
        
    if 'organism' not in p.annotations:
        return  ([],[])
    
    edge=[]
    
    #if db!='genbank':
    if 'db_source' in p.annotations:
        #if p.id[0:3]!="WP_":        
        if 'xrefs' in p.annotations['db_source']:
            m=re.search("xrefs: (.*?)[^,] .+[a-z]", p.annotations['db_source'])
            if m:
                edge = m.groups()[0].strip(".").split(", ")        
        elif 'accession ' in p.annotations['db_source']:
            edge=[p.annotations['db_source'].split("accession ")[-1]]

        else:
            print("no xrefs")
        


    id,seq_ver=p.id.split('.')

    node="\t".join([       
                p.id,  
                seq_ver,           
                prod,
                p.annotations['organism'],
                db,
                str(p.seq)
            ])
    
    if not edge and p.dbxrefs:
        edge=p.dbxrefs
            
    return ([node], [p.id+"\t"+e for e in edge])


def parse_gb_protein(p):

    if p.type!="CDS" or 'pseudo'  in p.qualifiers or 'pseudogene'  in p.qualifiers:
        return ([],[])

    if 'protein_id' not in p.qualifiers:
        return  ([],[])

    prod=""
    if 'product' in p.qualifiers:
        prod=p.qualifiers['product'][0] 

    if 'organism' not in r.annotations:
        return  ([],[])

    db=u.classify_acc(p.qualifiers['protein_id'][0])[0]
    if db!='genbank':
        a=1
        print("XXX")

    id,seq_ver=p.qualifiers['protein_id'][0].split('.')
    node="\t".join([       
                p.qualifiers['protein_id'][0],  
                seq_ver,            
                prod,
                p.annotations['organism'],
                u.classify_acc(p.qualifiers['protein_id'][0])[0],
                p.qualifiers['translation'][0]     
            ])
    return (node, [])


def parse_gb_nuc(r):
    # look at all CDS 
    #print(j)
    nodes=[]
    edges=[]

    for p in r.features:
        node,edge=parse_gb_protein(p)
        nodes=nodes+node
        edges=edges+edge
    return(nodes, edges)