# Read in python file

# Loop over records and record information
    # Extract Field
    # 
    # 


from Bio import SeqIO
import re
from collections import defaultdict
from progressbar import progressbar




input_file='../'

for p_rec in progressbar(SeqIO.parse(input_file, "genbank")):
    if sp_rec.description[0:12]!="MULTISPECIES":
        continue

    sp_rec.description
    id.append(sp_rec.id)

