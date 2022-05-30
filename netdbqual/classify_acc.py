import re 

def map_acc(acc, nuc, prot, nuc_str, prot_str):
    nuc_vals = [nuc_str]*len(nuc)
    prot_vals = [prot_str]*len(prot)
    type_map = dict(zip(nuc+prot,nuc_vals+prot_vals))
    #print(type_map)
    #print(acc)
    if acc not in type_map:
        return ("uniprot", "protein")
    return type_map[acc]


# Using the first 3 characters of the accession, 
# determine which type of refseq record this is
def classify_refseq_acc(acc):
    nuc=["AC_", "NC_", "NG_", "NT_", "NW_", "NZ_"]
    prot=["AP_", "NP_", "YP_", "XP_", "WP_"]
    return map_acc(acc[0:3], nuc, prot, ("refseq", "nucleotide"), ("refseq", "protein"))


def count_start_letters(acc):
    m=re.match("^[A-Z]+", acc)
    if m==None:
        return None
    return m.span()[1]

def count_end_num(acc):
    #Reverse to make counting easier
    m=re.match("^[0-9]+", acc[::-1])
    if m==None:
        return None
    return m.span()[1]


# Use patterns of letters/numners to classify
# See following for details
# https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
def classify_gb_acc(acc):
    acc=acc.split(".")[0]
    gb=[(1,5), (2,6), (2,8), (4,8)]
    prot=[(3,5), (3,7)]
    #
    n_alpha=count_start_letters(acc)
    #print(acc)
    n_num=count_end_num(acc)
    #print (n_alpha, n_num)
    return map_acc((n_alpha, n_num), gb, prot, ("genbank", "nucleotide"), ("genbank", "protein"))

def classify_acc(acc):

    acc_regexs={
       ("uniprot", "protein"):re.compile("^(([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9]))(.\\d+)?$"),
       ("genbank", "protein"):re.compile("^[A-Z]{3}(\\d{5}|\\d{7})(\\.\\d+)?$"),
       ("genbank", "nucleotide"):re.compile("^(([U]\\d{5})|([A-Z]{2}\\d{6}))(\\.\\d+)?$|([A-Z]\\d{5}\\.\\d+)"),
       ("genbank", "contig"):re.compile("^([A-Z]{4}\\d{8}(\\d+)?)|([A-Z]{6}\\d{9}(\\d+)?)$"),
       ("refseq", "protein"):re.compile("^((AP|NP|XP|YP|ZP)_\\d+)(\\.\\d+)?$"),
       ("refseq.anon", "protein"):re.compile("^((WP)_\\d+)(.\\d+)?$"),
       ("refseq", "nucleotide"):re.compile("^(((AC|NC|NM|NR|XM|XR)_\\d+)|(NZ\\_[A-Z]{1}\\d+))(.\\d+)?$"),
       ("refseq", "RNA"):re.compile("^(NM|NR|XM|XR)_\\d+"),
       ("refseq", "contig"):re.compile("^(((NG|NT|NW|NZ)_\\d+)|(NZ\\_[A-Z]{2,6}\\d+))(\\.\\d+)?$"),
       ("pdb", "protein"):re.compile("^[0-9][A-Za-z0-9]{3}(_([ABCLX]?))?$"),
       ("pir", "protein"):re.compile("^(([A-Z]{2}\\d{4})|([A-OR-TV-Z]{1}\\d{5}))$")  
    }
 
    for k,r in acc_regexs.items():
        if r.match(acc):
            return k
    return("unknown", "unknown")
