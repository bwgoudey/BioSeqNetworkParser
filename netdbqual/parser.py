from __future__ import annotations
import netdbqual.classify_acc as ca
import re
from typing import List, Tuple


def extractOrganism(r):
    org = ""
    # For a protein that is a genbank record subset
    if 'organism' in r.annotations:       
        org=r.annotations['organism']
    
    return(org)


def extractProduct(r, rec_type: str, db:str):
    # For a genbank record
    prod=""
    if rec_type == "gb":
        if r.description[0:7] == "RecName":
            prod = r.description.split(";")[0].split("=")[-1]
        else:
            prod = r.description.split(' [')[0]        
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        if 'product' in r.qualifiers:
            prod = r.qualifiers['product'][0]
    else:
        raise
    
    prod=prod.split("MULTISPECIES: ")[-1]
    return(prod)
    
def get_dbxrefs(dbx):
    return [id.split(":")[1] for id in dbx if id[0:6]=="RefSeq" or id[0:4]=="EMBL"]

def extractEdges(r, rec_type: str, seq_type):
    edge = []
    # For a genbank record
    if rec_type == "gb":
        if 'db_source' in r.annotations:
            # if r.id[0:3]!="WP_":
            if 'xrefs' in r.annotations['db_source']:
                m = re.search("xrefs: (.*?)[^,] .+[a-z]",
                              r.annotations['db_source'])
                if m:
                    edge = m.groups()[0].strip(".").split(", ")
            elif 'accession ' in r.annotations['db_source']:
                edge = [r.annotations['db_source'].split("accession ")[-1]]
        elif hasattr(r, 'dbxrefs'):
            edge=get_dbxrefs(r.dbxrefs)

    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        edge = []
    else:
        raise

    return edge


def extractSeqVersion(r, rec_type: str, seq_type) -> int:
    """Extract the sequence version from a given record. 
    If no version exists, return 

    Args:
        r (_type_): Record from Bio.Seq object
        db (str): String describing database {genbank, refseq, uniprot}
        seq_type (str): _description_

    Returns:
        int: version of sequence
    """
    seq_ver=""
    # For a genbank record
    if rec_type == "gb" :
        if  "." in r.id :
            id, seq_ver = r.id.split('.')
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        id, seq_ver = r.qualifiers['protein_id'][0].split('.')
    else:
        raise

    #If a version exists, it should be numeric. 
    #This is a hack to confirm this. Could use an assert
    if(seq_ver):
        seq_ver=int(seq_ver)

    return seq_ver


def extractSeq(r, rec_type: str, seq_type: str) -> str:
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
        seq = str(r.seq)
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        seq = str(r.qualifiers['translation'][0])
    else:
        raise
    return(seq)


def isRelevant(r, rec_type: str, seq_type: str) -> bool:
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


def extractSeq(r,rec_type:str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    raise 
    # For a genbank record
    if rec_type == "gb":
        1
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        1
    else:
        raise
    return(1)


def extractDateUpload(r,rec_type:str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    raise 
    # For a genbank record
    if rec_type == "gb":
        1
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        1
    else:
        raise
    return(1)


def extractDateModified(r,rec_type:str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    raise 
    # For a genbank record
    if rec_type == "gb":
        1
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        1
    else:
        raise
    return(1)



def extractNumProducts(r,rec_type:str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    raise 
    # For a genbank record
    if rec_type == "gb":
        1
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        1
    else:
        raise
    return(1)



def extractGO(r,rec_type:str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    raise 
    # For a genbank record
    if rec_type == "gb":
        1
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        1
    else:
        raise
    return(1)





def extractEC(r,rec_type:str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    raise 
    # For a genbank record
    if rec_type == "gb":
        1
    # For a protein that is a genbank record subset
    elif rec_type == "gb_sub":
        1
    else:
        raise
    return(1)











def parseRecord(r, db: str, seq_type: str) -> Tuple(List[str], List[str]):
    id = r.id
    if "UOI52910.1" in id:
        a=1

    if hasattr(r, 'annotations'):
        rec_type = 'gb'
    elif hasattr(r, 'qualifiers'):
        rec_type = "gb_sub"
    else:
        raise

    # If its not a CDS or if its a pseudoprotein
    if not isRelevant(r, rec_type, seq_type):
        return

    id = r.id
    prod=extractProduct(r, rec_type, db)
    organism=extractOrganism(r, rec_type)
    date_upload=extractDateUpload(r, rec_type)
    date_modified=extractDateModified(r, rec_type)
    n_products=extractNumProducts(r, rec_type)
    taxa_id=extractNumProducts(r, rec_type)
    go=extractGO(r, rec_type)
    ec=extractEC(r, rec_type)

    edge = extractEdges(r, rec_type, seq_type)
    seq_version=extractSeqVersion(r, rec_type, seq_type)
    seq = extractSeq(r, rec_type, seq_type)
    node = "\t".join([
        id,
        db,
        seq_type, 
        n_products, 
        date_upload, 
        date_modified, 
        go,
        ec,
        organism,
        taxa_id,
        prod,
        str(seq_version),
        seq
    ])
    if db=="uniprot":
        a=1
    return ([node], ["\t".join([id, e, seq_type,ca.classify_acc(e)[1]]) for e in edge])












