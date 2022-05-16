
def parse_gb_protein_gff(p):

    if 'molecule_type' not in p.annotations:
        return ([],[])


    prod=p.description.split(' [')[0] 

    if 'organism' not in p.annotations:
        return  ([],[])

    db=classify_acc(p.id)[0]
    if db!='genbank':
        a=1
        print("XXX")

    id,seq_ver=p.id.split('.')

    node=",".join([       
                id,  
                seq_ver,           
                prod,
                p.annotations['organism'],
                db,
                str(p.seq)
            ])

            
    return ([node], [])


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

    db=classify_acc(p.qualifiers['protein_id'][0])[0]
    if db!='genbank':
        a=1
        print("XXX")

    id,seq_ver=p.qualifiers['protein_id'][0].split('.')
    node=",".join([       
                id,  
                seq_ver,            
                prod,
                r.annotations['organism'],
                classify_acc(p.qualifiers['protein_id'][0])[0],
                p.qualifiers['translation'][0]     
            ])
    return (node, [])


def parse_gb_nuc(r):
    # look at all CDS 
    print(j)
    nodes=[]
    edges=[]

    for p in r.features:
        node,edge=parse_gb_protein(p)
        nodes=nodes+node
        edges=edges+edges
    return(nodes, edges)