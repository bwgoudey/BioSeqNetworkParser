
import re
try:
    from netdbqual.classify_acc import classify_acc 
except:
    from classify_acc import classify_acc





try:
    from netdbqual.parser_utils import identifyProteins 
except:
    from parser_utils import identifyProteins



def parseInference(inf, note=""):
    no_inf="non-experimental evidence, no additional details recorded"
    if inf[0:20]=="COORDINATES: similar" or inf[0:22]=="similar to AA sequence":
        match=inf.split("sequence:")
        if len(match)!=2:
            return
            #raise RuntimeError("Unexpected pattern")

        if "," in match[1]:
            #TODO: we just trim this for now. 
            match[1] = match[1].split(",")[0]
        
        if ":" in match[1]:
            db, acc_ver = match[1].split(":")
        else:
            acc_ver=match[1]
            db=classify_acc(acc_ver)

        if "." in acc_ver:
            acc, seq_version=acc_ver.split(".")
            seq_version=int(seq_version)
        else:
            acc=acc_ver
            seq_version=""
        return({'acc':acc,
                'seq_version':seq_version,
                'identity': "",
                'db':db, 
                'type':'Hom',
                 'model':""})

    elif inf[0:20]=="COORDINATES: ab init":
        match=inf.split("prediction:")
        if len(match)!=2:
            raise RuntimeError("Unexpected pattern")
        model = match[1]    
        return({'acc':"",
                'seq_version':"",
                'identity': "",
                'db':'', 
                'type':'AbInit', 
                'model':model})
                
    elif inf[0:21]=="ab initio prediction:":
        match=inf.split("prediction:")
        if len(match)!=2:
            raise RuntimeError("Unexpected pattern")
        model = match[1]    
        return({'acc':"",
                'seq_version':"",
                'identity': "",
                'db':'', 
                'type':'AbInit', 
                'model':model})

    elif inf[0:14]=="protein motif:" or inf[0:27]=="COORDINATES: protein motif:":
        match=inf.split("motif:")
        if len(match)!=2:
            raise RuntimeError("Unexpected pattern")
        if "," in match[1]:
            #TODO: we just trim this for now. 
            match[1] = match[1].split(",")[0]
        db_acc = match[1].split(":")
        if len(db_acc) == 2:
            db, acc = db_acc
        else:
            raise RuntimeError("Unexpected pattern")
        

        return({'acc':acc,
                'seq_version':"",
                'identity': "",
                'db':db, 
                'type':'Motif', 
                'model':""})

    elif inf==no_inf and (not note or len(note)<6):
        return        
    else:
        return#raise RuntimeError("Not implemented")
    return 


def extractInferenceEdgesP(p, db, p_id, p_seq_ver):
    if db=="uniprot":
        raise 

    inferences=[]

    relevant_feats=set(p.keys()).intersection(set(['CDS', 'Protein']))
    # For a genbank record
    
    for f in relevant_feats:
        if 'inference' in p[f].qualifiers:
            note = p[f].qualifiers['note'][0] if 'note' in p[f].qualifiers else ""
            inf=p[f].qualifiers['inference']
            if note and len(p[f].qualifiers['note'])>1:
                raise RuntimeError("Inference is returning some list")
            for i in inf:
                inferences.append(parseInference(i, note)) 
    for i in inferences:
        if i is None:
            continue
        i['trg_id']=p_id
        i['trg_seq_ver']=p_seq_ver
        

    return inferences

def extractInferenceEdges(r, db):
    if db=="uniprot":
       raise 

    inferences=[]
    ps=list(identifyProteins(r).values())
    for p in ps: 
        relevant_feats=set(p.keys()).intersection(set(['CDS', 'Protein']))
        # For a genbank record
        
        for f in relevant_feats:
            if 'inference' in p[f].qualifiers:
                note=""
                if 'inference' in p[f].qualifiers:
                    note = p[f].qualifiers['note'][0] if 'note' in p[f].qualifiers else ""
                    if note and len( p[f].qualifiers['note']) > 1:
                        raise RuntimeError("Inference is returning some list")
                    inf=p[f].qualifiers['inference']
                    #if(len(inf)>1):
                    #    raise RuntimeError("Inference is returning some list")
                    for i in inf:
                        inferences.append(parseInference(i, note)) 
    return inferences
