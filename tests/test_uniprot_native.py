from sqlite3 import dbapi2
import unittest   # The test framework
from importlib import reload  # Python 3.4+
from netdbqual import parser
from Bio import SeqIO
from io import StringIO

class Test_Parser(unittest.TestCase):
    
    def setUp(self):
        self.rec_str="""ID   001R_FRG3G              Reviewed;         256 AA.
AC   Q6GZX4;
DT   28-JUN-2011, integrated into UniProtKB/Swiss-Prot.
DT   19-JUL-2004, sequence version 1.
DT   02-JUN-2021, entry version 42.
DE   RecName: Full=Putative transcription factor 001R;
DE            EC=2.3.1.122;
DE            EC=2.3.1.20;
GN   ORFNames=FV3-001R;
OS   Frog virus 3 (isolate Goorha) (FV-3).
OC   Viruses; Varidnaviria; Bamfordvirae; Nucleocytoviricota; Megaviricetes;
OC   Pimascovirales; Iridoviridae; Alphairidovirinae; Ranavirus.
OX   NCBI_TaxID=654924;
OH   NCBI_TaxID=30343; Dryophytes versicolor (chameleon treefrog).
OH   NCBI_TaxID=8404; Lithobates pipiens (Northern leopard frog) (Rana pipiens).
OH   NCBI_TaxID=45438; Lithobates sylvaticus (Wood frog) (Rana sylvatica).
OH   NCBI_TaxID=8316; Notophthalmus viridescens (Eastern newt) (Triturus viridescens).
RN   [1]
RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].
RX   PubMed=15165820; DOI=10.1016/j.virol.2004.02.019;
RA   Tan W.G., Barkman T.J., Gregory Chinchar V., Essani K.;
RT   "Comparative genomic analyses of frog virus 3, type species of the genus
RT   Ranavirus (family Iridoviridae).";
RL   Virology 323:70-84(2004).
CC   -!- FUNCTION: Transcription activation. {ECO:0000305}.
CC   ---------------------------------------------------------------------------
CC   Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms
CC   Distributed under the Creative Commons Attribution (CC BY 4.0) License
CC   ---------------------------------------------------------------------------
DR   EMBL; AY548484; AAT09660.1; -; Genomic_DNA.
DR   RefSeq; YP_031579.1; NC_005946.1.
DR   SwissPalm; Q6GZX4; -.
DR   GeneID; 2947773; -.
DR   KEGG; vg:2947773; -.
DR   Proteomes; UP000008770; Genome.
DR   GO; GO:0046782; P:regulation of viral transcription; IEA:InterPro.
DR   InterPro; IPR007031; Poxvirus_VLTF3.
DR   Pfam; PF04947; Pox_VLTF3; 1.
PE   4: Predicted;
KW   Activator; Reference proteome; Transcription; Transcription regulation.
FT   CHAIN           1..256
FT                   /note="Putative transcription factor 001R"
FT                   /id="PRO_0000410512"
SQ   SEQUENCE   256 AA;  29735 MW;  B4840739BF7D4121 CRC64;
     MAFSAEDVLK EYDRRRRMEA LLLSLYYPND RKLLDYKEWS PPRVQVECPK APVEWNNPPS
     EKGLIVGHFS GIKYKGEKAQ ASEVDVNKMC CWVSKFKDAM RRYQGIQTCK IPGKVLSDLD
     AKIKAYNLTV EGVEGFVRYS RVTKQHVAAF LKELRHSKQY ENVNLIHYIL TDKRVDIQHL
     EKDLVKDFKA LVESAHRMRQ GHMINVKYIL YQLLKKHGHG PDGPDILTVK TGSKGVLYDD
     SFRKIYTDLG WKFTPL
//"""
        self.rec=SeqIO.parse(StringIO(self.rec_str), "swiss").__iter__().__next__()
        self.db='uniprot'
        self.seq_type="protein"
        self.rec_type="top"
        self.dbsource_str=""

        ##########
         # What sort of record? Its a protein, described at top-level 
    # (i.e whole record is about this protein))
    def test_determineRecordType(self):
        obs_rec_type=parser.determineRecordType(self.rec)
        exp_rec_type="top"
        self.assertEqual(obs_rec_type, exp_rec_type)        


    # Is it a real protein? Discard Psuedo proteins
    def test_is_pseudo(self):
        obs_pseudo=parser.isPseudo(self.rec, self.rec_type, self.seq_type)
        exp_pseudo=False
        self.assertEqual(obs_pseudo, exp_pseudo)  

    def test_extractName(self):
        obs_name=parser.extractDescription(self.rec.description)
        exp_name='Putative transcription factor 001R'
        self.assertEqual(obs_name, exp_name)

    # What organism are we looking at
    def test_extractOrganism(self):
        obs_pseudo=parser.extractOrganism(self.rec.annotations)
        exp_pseudo="Frog virus 3 (isolate Goorha) (FV-3)"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    def test_extractTaxaID(self):
        obs_taxid=parser.extractTaxaID(self.rec.features[0], self.rec.annotations)
        exp_taxid="654924"
        self.assertEqual(obs_taxid, exp_taxid)

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_extractTaxonomy(self):
        obs_pseudo=parser.extractTaxonomy(self.rec.annotations)
        exp_pseudo='Viruses,Varidnaviria,Bamfordvirae,Nucleocytoviricota,Megaviricetes,Pimascovirales' 
        self.assertEqual(obs_pseudo, exp_pseudo)  

    #Figure out when thi was first submitted
    def test_extractDateModified(self):
        obs=parser.extractDateModified(self.rec, "")
        exp_upload=20210602
        exp_last_mod=20210602
        exp_nmodify=-1
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)   


    def test_extractDateLastModified(self):
        obs_date=parser.extractDateLastModified(self.rec.annotations)
        exp_date=20210602
        self.assertEqual(obs_date, exp_date)  

  
    def test_extractGo(self):
        prot1=""#{f.type:f for f in (r.features[1],r.features[4])}
        obs_go=parser.extractGO(self.rec, prot1, self.db)
        #obs_go=parser.extractGO(self.rec, self.rec_type)
        exp_go=["GO:0046782"]
        self.assertEqual(obs_go, exp_go)  

    def test_extractEC(self):
        obs_ec=parser.extractEC(self.rec, "")
        exp_ec=['2.3.1.122', '2.3.1.20']
        self.assertEqual(obs_ec, exp_ec)  

    #There are some limitations to the current UNiProt parser. In particular, it seems to only capture a subset of accesssions. 
    def test_extractParentEdges(self):
        obs_edges=parser.extractParentEdges(self.rec,self.db, "")
        #exp_edges=["AY548484.1", "AAT09660.1", "YP_031579.1", "NC_005946.1"]
        exp_edges=sorted(["AY548484",  "YP_031579.1"])
        self.assertEqual(sorted(obs_edges), exp_edges)        


    def test_createTopLevelNode(self):
            
            obs=parser.createTopLevelNode(self.rec, self.rec_type, self.seq_type, self.db,self.dbsource_str)
            obs['n']['seq']=obs['n']['seq'][0:10]
            
            exp_node={'id': 'Q6GZX4', 
                'seq_version': '',
                'date_first_upload': 20210602,
                'num_modified': -1, 
                'date_last_modified': 20210602,
                'organism': 'Frog virus 3 (isolate Goorha) (FV-3)', 
                'taxa_id': '654924', 
                'taxonomy': 'Viruses,Varidnaviria,Bamfordvirae,Nucleocytoviricota,Megaviricetes,Pimascovirales',
                'n_products': 0, 
                'name': 'Putative transcription factor 001R',
                'go': ['GO:0046782'], 
                'ec': ['2.3.1.122', '2.3.1.20'], 
                'seq': 'MAFSAEDVLK'}

            self.assertEqual(obs['n'], exp_node)

            exp_edges=[('Q6GZX4', 'AY548484'), ('Q6GZX4', 'YP_031579.1')]
            self.assertEqual(obs['e'], exp_edges)


    def test_uniprot_parseRecord(self):
        obs_node, obs_edge=parser.parseRecord(self.rec, self.db, self.seq_type)
        exp_edge=[('Q6GZX4', 'AY548484'), ('Q6GZX4', 'YP_031579.1')]
        exp_nodes=1
        self.assertEqual(len(obs_node), exp_nodes)
        self.assertEqual(obs_edge, exp_edge)        