from docopt import docopt
from Bio import SeqIO
from io import StringIO
import sys

def main():
    doc =  """genbank2fasta

    Usage: genbank2fasta.py [-o FASTA]  [-i GBK]

    Options:
        -o FASTA --output FASTA  FASTA file to write output to
        -i GBK   --input GBK     input Genbank file
        

    """

    args = docopt(doc, argv=None, help=True, version=None, options_first=False)
    handle = StringIO("")
    if 'FASTA' in args: 
        handle=args['FASTA']

    
    if 'GBK' in args: 
        input=args['GBK']
    else:
        input = sys.stdin


    SeqIO.convert(input, "genbank", handle, "fasta")   
    print(handle.getvalue())

if __name__ == "__main__":
    main()


