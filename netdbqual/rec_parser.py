from __future__ import annotations
from html import entities
from os import get_blocking
from ssl import HAS_ALPN
from xml.dom.minicompat import NodeList
#from types import NoneType

import classify_acc as ca
import re
from typing import List, Tuple
from Bio import SeqRecord
from datetime import datetime
from collections import defaultdict
import copy

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


def extractTaxaID(source_feat, r_annot={}):
    if 'ncbi_taxid' in r_annot: 
        return r_annot['ncbi_taxid'][0]
    #taxid={f.type:f for f in r.features}['source'].f.
    if source_feat.type!="source":
        raise
    if 'db_xref' in source_feat.qualifiers:
        return source_feat.qualifiers['db_xref'][0].split(":")[1]
        return ""


#
# Dates
#

def gpkDateStrToNum(date_str):
    d=datetime.strptime(date_str, '%d-%b-%Y')
    return(d.year*10000+d.month*100+d.day)

def uniprotDateStrToNum(date_str):
    d=datetime.strptime(date_str, '%b %d, %Y')
    return(d.year*10000+d.month*100+d.day)



def extractDateFromJournalStr(date_str):
    date_regex=re.compile(r"Submitted \((\d{2}-[A-Z]{3}-[0-9]{4})\)")
    upload_date=date_regex.match(date_str)
    if(not upload_date):
        return -1
    
    return(gpkDateStrToNum(upload_date.groups()[0]))
    

def extractDateModified(r,uniprot_dbsource="") -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
     
    # For a genbank record
    modify_dates=[]
    #We make the assumption that the first indication of 
    if uniprot_dbsource:
        upload=uniprotDateStrToNum(uniprot_dbsource['created'])
        modify=uniprotDateStrToNum(uniprot_dbsource['annotation updated'])
        return({
                "first_upload":upload,
                "last_modify":modify,
                "num_modify":-1 })    
    elif 'date_last_annotation_update' in r.annotations:
        upload=gpkDateStrToNum(r.annotations['date_last_annotation_update'])
        modify=gpkDateStrToNum(r.annotations['date_last_annotation_update'])
        return({
                "first_upload":upload,
                "last_modify":modify,
                "num_modify":-1 })    
    elif 'references' in r.annotations and r.annotations['references']:
        possible_dates=[extractDateFromJournalStr(s.journal) for s in r.annotations['references']]
        modify_dates=list(filter(lambda x: x>0, possible_dates))
    
    if r.id[0:3]=='WP_' or not modify_dates:
        return({
                "first_upload":"",
                "last_modify":"",
                "num_modify":-1 }) 
    #else:
    #    raise

    return({"first_upload":modify_dates[-1],
    "last_modify":modify_dates[0],
    "num_modify":len(modify_dates) })
    


def extractDateLastModified(r_annot) -> str:
    """Extract the sequence from a given record

    Args:
        r (_type_): Record from Bio.Seq object
        rec_type (tr): describes the type of record we are considering
    Returns:
        str: 
    """
    if 'date_last_annotation_update' in r_annot:
        return gpkDateStrToNum(r_annot['date_last_annotation_update'])
    # For a genbank record
    return gpkDateStrToNum(r_annot['date'])


def extractDescription(r_desc):
    # For a genbank record
    desc=""
    if r_desc[0:7] == "RecName":
        desc = r_desc.split(";")[0].split("=")[-1]
    else:
        desc = r_desc.split(' [')[0]       
    
    desc=desc.split("MULTISPECIES: ")[-1]
    return desc




def extractProduct(feature_dict, seq_type):
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
    if seq_type=="nucleotide" or seq_type=="contig":
        feat=feature_dict['CDS']
    else:
        feat=feature_dict['Protein']

    id=""
    product_name=""
    if 'protein_id' in feat.qualifiers:
        id=feat.qualifiers['protein_id'][0]
    if 'product' in feat.qualifiers:
        product_name=feat.qualifiers['product'][0]

    return (id, product_name)
    

def get_dbxrefs(dbx):
    return [id.split(":")[1] for id in dbx if id[0:6]=="RefSeq" or id[0:4]=="EMBL"]


def extractRefSeqParentEdge(comment_str):
    comment_str=comment_str.replace('\n', ' ').replace(' and ', ', ')
    m=re.search(
        "The reference sequence was derived from ([_A-Z0-9\.,and ]+)\.",  comment_str)
    if m:
        return(m.groups()[0].split(", "))
    m=re.search(
        "The reference sequence is identical to ([_A-Z0-9\.,and ]+)\.", comment_str)
    if m:
        return(m.groups()[0].split(", "))
    m=re.search(
        "This record is derived from a genomic sequence \(([_A-Z0-9\.,and ]+)\)", comment_str)
    if m:
        return(m.groups()[0].split(", "))        
    else:
        raise

    


def extractParentEdges(r, db, uniprot_dbsource=""):
    r_annot=r.annotations
    if db=="refseq":
        return(extractRefSeqParentEdge(r_annot['comment']))
    
    edge=[]
    if uniprot_dbsource and 'xrefs' in uniprot_dbsource:
        edge=uniprot_dbsource['xrefs'].strip().split(', ')
    elif 'db_source' in r_annot:
            # # if r.id[0:3]!="WP_":
            # if 'xrefs' in r_annot['db_source']:
            #     m = re.search("xrefs: (.*?)[^,] .+[a-z]",
            #                   r_annot['db_source'])
            #     if m:
            #         edge = m.groups()[0].strip(".").split(", ")
            if 'accession ' in r_annot['db_source']:
                edge = [r_annot['db_source'].split("accession ")[-1]]
    elif hasattr(r, 'dbxrefs'):
            edge=get_dbxrefs(r.dbxrefs)

    return edge

def extractEdges(r, rec_type, node, db, ):
    edge=extractParentEdges(r, db)
    for n in node[1:]:
        [r.id, ]
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
    return seq


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
    return {k:v for k,v in fd.items() if 'CDS' in v or 'Protein' in v}






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
    gos=[re.findall(r'(GO:\d{7})', p[f].qualifiers['note'][0]) for f in relevant_feats if 'note' in p[f].qualifiers]
    gos=sorted(list(set([x for go in gos for x in go])))

    return gos

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
    if db=="uniprot":
        if hasattr(r, 'ncbi_taxid'):
            return [x[3:] for x in r.dbxrefs if x[0:3]=='GO:']
        elif 'db_source' in r.annotations:
            return extractUniprotGo(r.annotations['db_source'])
        else:
            return []   
    return extractNcbiGO(p)



#def extractEC(r,rec_type:str) -> str:
def extractEC(r, product) -> str:
    """Extract the sequence from a given record
    """
    if product and 'EC_number' in product.qualifiers:
        ec=product.qualifiers['EC_number']
        #assert(len(ec)==1)
        return (','.join(ec))

    #TODO: This fails a lot. Need a better test
    uniprot_format = hasattr(r, 'dbxrefs') and len(r.dbxrefs)
    if uniprot_format: 
        desc=r.description
        tokens= desc.split("; ")
        ECs=[l[3:].rstrip(";") for l in tokens if l[0:3]=="EC="]
        if not ECs:
            return ""
        return ECs


    return ""


def determineRecordType(r):
    if hasattr(r, 'annotations'):
        rec_type = 'top'
    elif hasattr(r, 'qualifiers'):
        rec_type = "product"
    else:
        raise

    return rec_type

def extractChildren(r, parent, seq_type, db):
    ps=identifyProteins(r)
    
    nodes=[]
    edges=[]

    if seq_type=="nucleotide" or seq_type == "contig":
        #product_field={'nucleotide':'CDS', 'protein':'Protein'}
        for p in ps.values():
            if ('gene' in p and 'pseudo' in p['gene'].qualifiers) or \
                ('CDS' in p and 'pseudo' in p['CDS'].qualifiers or \
                'pseudogene' in p['CDS'].qualifiers ):
                continue

            child=copy.deepcopy(parent)

            child['id'], child['name']=extractProduct(p, seq_type)
            child['go']=extractGO(r, p, db)
            child['ec']=extractEC(r, p['CDS'])
            child['n_products']=0
            if 'translation' in p['CDS'].qualifiers:
                child['seq']=str(p['CDS'].qualifiers['translation'][0])
            #child['parent']=parent['id']
            #extractEC(p[product_field[rec_type]])
            # if len(ps)>1:
            #     child['id']=p['CDS'].qualifier['protein_id']
            #     child['seq'] = extractSeq(r, rec_type, seq_type)
            #     echild_p['parent']=[r.id]
            #     entities=[e_p]
            #     break
            nodes.append(child.values())
            edges.append((child['id'], parent['id']))
    elif len(ps)>1:
        raise
    return {"n":nodes, "e":edges}

def createTopLevelNode(r, rec_type, seq_type, db,uniprot_dbsource=""):
    e={}
    e['id'] = r.id
    e['seq_version']=extractSeqVersion(r, rec_type, seq_type)

    modified_info=extractDateModified(r,uniprot_dbsource)
    e['date_first_upload']=modified_info['first_upload']
    e['num_modified']=modified_info['num_modify']
    e['date_last_modified']=extractDateLastModified(r.annotations)

    e['organism']=extractOrganism(r.annotations)
    e['taxa_id']=extractTaxaID(r.features[0], r.annotations)
    e['taxonomy']=extractTaxonomy(r.annotations)
    
    #n_products=extractNumProducts(r, rec_type)

    proteins=identifyProteins(r)
    e['n_products']=len(proteins)

    e['name']=extractDescription(r.description)
    
    if seq_type=="protein":
        p=list(identifyProteins(r).values())
        if len(p)>1:
            raise
        
        source_annot=""
        if len(p)==1:
            e['go']=extractGO(r, p[0], db)
            e['ec']=extractEC(r, p[0]['Protein'])
        elif len(p)==0:
            e['go']=extractGO(r, "", db)
            e['ec']=extractEC(r,  "")
    else:
        e['go']=""
        e['ec']=""

    
    try:
        e['seq']=str(r.seq)
    except:
        e['seq']=""
    edges=extractParentEdges(r,db,uniprot_dbsource)
    edges=[(e['id'], edge) for edge in edges]
    return {'n':e, 'e':edges}



def processUniProtDBsource(dbsource_str):
    preprocess_str=dbsource_str.replace("; ", ". ").replace("xrefs (n", ". xrefs (n")
    preprocess_str=preprocess_str.replace(' created', '. created').replace('extra accessions:', 'extra accessions: ').split(". ")
    return dict([x.split(": ") for x in preprocess_str if x])



def parseRecord(r, db: str, seq_type: str) -> Tuple(List[str], List[str]):
    id = r.id
    if "UOI52910.1" in id:
        a=1

    rec_type=determineRecordType(r)

    # If its not a CDS or if its a pseudoprotein
    if isPseudo(r, rec_type, seq_type):
        return
    uniprot_dbsource=""
    if db=="uniprot" and 'db_source' in r.annotations:
        uniprot_dbsource=processUniProtDBsource(r.annotations['db_source'])
   
    parent=createTopLevelNode(r, rec_type, seq_type, db,uniprot_dbsource)
    children=extractChildren(r, parent['n'], seq_type, db)
    #edge = extractEdges(r, rec_type, nodes)
    node_strs=[parent['n'].values()]+children['n']
    edge_strs=parent['e']+children['e']

    #return ([node], ["\t".join([id, e, seq_type,ca.#classify_acc(e)[1]]) for e in edge])
    return node_strs, edge_strs, list(parent['n'].keys())












