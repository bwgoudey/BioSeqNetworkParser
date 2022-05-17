from __future__ import annotations
import netdbqual.utils as u
import re 
from typing import List, Tuple


def extractProduct(r, rec_type:str, seq_type):
    # For a genbank record
    if rec_type == "gb":
        1
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        1
    else:
        raise

def extractOrganism(r, rec_type:str, seq_type):
    # For a genbank record
    if rec_type == "gb":
        1
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        1
    else:
        raise


def extractEdges(r, rec_type:str, seq_type):
    edge=[]
    # For a genbank record
    if rec_type == "gb":
        if 'db_source' in r.annotations:
            #if r.id[0:3]!="WP_":        
            if 'xrefs' in r.annotations['db_source']:
                m=re.search("xrefs: (.*?)[^,] .+[a-z]", r.annotations['db_source'])
                if m:
                    edge = m.groups()[0].strip(".").split(", ")        
            elif 'accession ' in r.annotations['db_source']:
                edge=[r.annotations['db_source'].split("accession ")[-1]]
            
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        edge=[]
    else:
        raise
    
    return edge



def extractSeqVersion(r, rec_type:str, seq_type):
    # For a genbank record
    if rec_type == "gb":
        id,seq_ver=r.id.split('.')
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        id,seq_ver=r.qualifiers['protein_id'][0].split('.')
    else:
        raise

    return seq_ver

def extractSeq(r, rec_type:str, seq_type:str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        db (str): String describing database {genbank, refseq, uniprot}
        seq_type (str): _description_

    Returns:
        str: nucleotide or amino acid sequence
    """ 
    # For a genbank record
    if rec_type == "gb":
        seq=str(r.seq)
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        seq=str(r.qualifiers['translation'][0])
    else:
        raise
    return(seq)


def isRelevant(r, rec_type:str, seq_type:str) -> bool:
    """Determine whether a record is one that we want
     to extract nodes/edges from. Irrelevant records are ignored by our analysis. 

    Args:
        r (_type_): Record from Bio.Seq object
        db (str): String describing database {genbank, refseq, uniprot}
        seq_type (str): _description_

    Returns:
        bool: True/False depending on record's relevance
    """    
    
    # For a genbank record
    if rec_type == "gb":
        if 'molecule_type' not in r.annotations:
            return False
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub" and seq_type == "protein":
        if (
            r.type != "CDS" or
            'pseudo' in r.qualifiers or
            'pseudogene' in r.qualifiers or
            'protein_id' not in r.qualifiers
        ):
            return False
    # throw error if case not captured
    else:
        raise

    return True


def parseRecord(r, rec_type:str, seq_type:str)-> Tuple(List[str], List[str]):
    if hasattr(r, 'annotations'):
        rec_type='gb'
    elif hasattr(r, 'qualifiers'):
        rec_type="gb_sub"
    else:
        raise

    #If its not a CDS or if its a pseudoprotein
    if not isRelevant(r, rec_type, seq_type):
        return
    
    id=r.id
    #prod=extractProduct(r, rec_type, seq_type)
    #organism=extractOrganism(r, rec_type, seq_type)
    edge=extractEdges(r, rec_type, seq_type)
    #seq_version=extractSeqVersion(r, rec_type, seq_type)
    seq=extractSeq(r, rec_type, seq_type)
    node  = "\t".join([       
                id,  
                #seq_version,           
                #prod,
                #organism,
                db,
                seq
            ])
    return ([node], [id+"\t"+e for e in edge])    

            




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
                node,edge=RecNode(p, db, seq_type).getNodeEdgesStr()
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