from collections import defaultdict

def identifyProteins(r):
    feats=r.features[1:]
    fd=defaultdict(dict)
    for f in feats:
        key=(f.location.start, f.location.end)
        fd[key][f.type]=f
    return {k:v for k,v in fd.items() if 'CDS' in v or 'Protein' in v}
