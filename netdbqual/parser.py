from __future__ import annotations
#from types import NoneType
import netdbqual.classify_acc as ca
import re
from typing import List, Tuple
from Bio import SeqRecord
from datetime import datetime
from collections import defaultdict

#
# Taxa
#



def extractOrganism(r_annots):
    org = ""
    # For a protein that is a genbank record subset
    if 'organism' in r_annots:       
        org=r_annots['organism']
    
    return(org)


def extractTaxonomy(r_annots, n_level=6):
    tax=[]
    if 'taxonomy' in r_annots:
        tax=r_annots['taxonomy'][0:min(len(r_annots['taxonomy']), n_level)]
    
    tax=",".join(tax+[""]*(n_level-len(tax)))
    return tax


def extractTaxaID(source_feat):

    #taxid={f.type:f for f in r.features}['source'].f.
    if source_feat.type!="source":
        raise
    return(source_feat.qualifiers['db_xref'][0])


#
# Dates
#

def gpkDateStrToNum(date_str):
    d=datetime.strptime(date_str, '%d-%b-%Y')
    return(d.year*10000+d.month*100+d.day)


def extractDateFromJournalStr(date_str):
    date_regex=re.compile(r"Submitted \((\d{2}-[A-Z]{3}-[0-9]{4})\)")
    upload_date=date_regex.match(date_str)
    if(not upload_date):
        return -1
    
    return(gpkDateStrToNum(upload_date.groups()[0]))
    

def extractDateModified(r,rec_type:str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
     
    # For a genbank record
    modify_dates=[]
    if rec_type == "top":
        #We make the assumption that the first indication of 
        if 'references' in r.annotations and r.annotations['references']:
            possible_dates=[extractDateFromJournalStr(s.journal) for s in r.annotations['references']]
            modify_dates=list(filter(lambda x: x>0, possible_dates))
        else:
            raise
    # For a protein that is a genbank record subset
    elif rec_type == "product":
        1
    else:
        raise
    return([modify_dates[-1], modify_dates[0], len(modify_dates)])


def extractDateLastModified(date_str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    
    # For a genbank record
    return(gpkDateStrToNum(date_str))







def extractProduct(feature_dict, rec_type):
    # # For a genbank record
    # prod=""
    # if rec_type == "top":
    #     if r.description[0:7] == "RecName":
    #         prod = r.description.split(";")[0].split("=")[-1]
    #     else:
    #         prod = r.description.split(' [')[0]        
    # # For a protein that is a genbank record subset
    # elif rec_type == "product":
    #     if 'product' in r.qualifiers:
    #         prod = r.qualifiers['product'][0]
    # else:
    #     raise
    
    # prod=prod.split("MULTISPECIES: ")[-1]
    # return(prod)
    if rec_type=="nucleotide":
        feat=feature_dict['CDS']
    else:
        feat=feature_dict['Protein']

    return(feat.qualifiers['product'])
    

def get_dbxrefs(dbx):
    return [id.split(":")[1] for id in dbx if id[0:6]=="RefSeq" or id[0:4]=="EMBL"]


def extractParentEdges(r_annot):
    if 'db_source' in r_annot:
            # if r.id[0:3]!="WP_":
            if 'xrefs' in r_annot['db_source']:
                m = re.search("xrefs: (.*?)[^,] .+[a-z]",
                              r_annot['db_source'])
                if m:
                    edge = m.groups()[0].strip(".").split(", ")
            elif 'accession ' in r_annot['db_source']:
                edge = [r_annot['db_source'].split("accession ")[-1]]
    elif hasattr(r, 'dbxrefs'):
            edge=get_dbxrefs(r.dbxrefs)


def extractRelatedEdges(cds_features):

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
    if rec_type == "top" :
        if  "." in r.id :
            id, seq_ver = r.id.split('.')
    # For a protein that is a genbank record subset
    elif rec_type == "product":
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
    if rec_type == "top":
        seq = str(r.seq)
    # For a protein that is a genbank record subset
    elif rec_type == "product":
        seq = str(r.qualifiers['translation'][0])
    else:
        raise
    return(seq)


def isPseudo(r: SeqRecord, rec_type: str, seq_type: str) -> bool:
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
    if rec_type == "top":
        if 'molecule_type' not in r.annotations:
            return False
    # For a protein that is a genbank record subset
    elif rec_type == "product" and seq_type == "protein":
        if (
            r.type != "CDS" or
            'pseudo' in r.qualifiers or
            'pseudogene' in r.qualifiers or
            'protein_id' not in r.qualifiers
        ):
            return True
    # throw error if case not captured
    
    else:
        if rec_type not in ["top","product"]:
            print("Bad rec_type specified not (top, product)")
        raise

    return False


def identifyProteins(r):
    feats=r.features[1:]
    fd=defaultdict(dict)
    for f in feats:
        key=(f.location.start, f.location.end)
        fd[key][f.type]=f
    return({k:v for k,v in fd.items() if 'CDS' in v})


def extractNumProducts(r,rec_type:str) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    
    # For a genbank record
    if rec_type == "top":
        return(len([f for f in r.features if f.type=="Protein" or f.type=="CDS"]))    # For a protein that is a genbank record subset
    elif rec_type == "product":
        1
    else:
        raise
    return(1)



def extractNcbiGO(p) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    gos=[]
    relevant_feats=set(p.keys()).intersection(set(['CDS', 'Protein']))

    # For a genbank record
    gos=[re.findall(r'(GO:\d{7})', f.qualifiers['note'][0]) for f in [p['CDS'], ] if 'note' in f.qualifiers]
    gos=sorted(list(set([x for go in gos for x in go])))

    return(gos)

def extractUniprotGo(dbsource) -> str:
    """Extract the sequence from a given record

    Args:
        dbsource (str): string 
    Returns:
        str: 
    """

    # For a genbank record
    return(sorted(re.findall(r'(GO:\d{7})',dbsource)))
    

def extractGO(r, p, db):
    if(db=="uniprot"):
        return extractUniprotGo(r.annotations['dbsource'])
    return extractNcbiGO(p)



#def extractEC(r,rec_type:str) -> str:
def extractEC(product) -> str:
    """Extract the sequence from a given record
    """
    if 'EC_number' in product.qualifiers:
        ec=product.qualifiers['EC_number']
        assert(len(ec==1))
        return (ec[0])
    return("")


def determineRecordType(r):
    if hasattr(r, 'annotations'):
        rec_type = 'top'
    elif hasattr(r, 'qualifiers'):
        rec_type = "product"
    else:
        raise

    return(rec_type)







def parseRecord(r, db: str, seq_type: str) -> Tuple(List[str], List[str]):
    id = r.id
    if "UOI52910.1" in id:
        a=1

    rec_type=determineRecordType(r)

    # If its not a CDS or if its a pseudoprotein
    if isPseudo(r, rec_type, seq_type):
        return

    id = r.id
    organism=extractOrganism(r.annotations)
    taxonomy=extractTaxonomy(r.annotations)
    taxa_id=extractTaxaID(r.features[0])

    date_upload=extractDateModified(r, rec_type)
    date_modified=extractDateLastModified(r, rec_type)
    seq_version=extractSeqVersion(r, rec_type, seq_type)
    #n_products=extractNumProducts(r, rec_type)

    proteins=identifyProteins(r.annotations)
    n_products=len(proteins)

    core_fields="\t".join([
        id, 
        seq_version,
        date_upload, 
        date_modified, 
        seq_version,
        organism, 
        taxa_id,
        taxonomy, 
        n_products
    ])

    nodes=[]
    if rec_type=="nucleotide":
        nodes.append([
            core_fields,
            "","","",
            str(r.seq)
        ])

    for p in identifyProteins(r.annotations):
        prod=extractProduct(p, rec_type)
        go=extractGO(r, p, db)
        ec=extractEC(r, rec_type)
        seq = extractSeq(r, rec_type, seq_type)
        nodes.append("\t".join([prod, go, ec, seq]))

  

    parent_edge=extractParentEdges(nodes, r.annotations, db)
    edge = extractEdges(r, rec_type)


    return ([node], ["\t".join([id, e, seq_type,ca.classify_acc(e)[1]]) for e in edge])












