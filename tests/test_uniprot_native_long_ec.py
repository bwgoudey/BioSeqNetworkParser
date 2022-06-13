from sqlite3 import dbapi2
import unittest   # The test framework
from importlib import reload  # Python 3.4+
from netdbqual import rec_parser
from Bio import SeqIO
from io import StringIO

class Test_Parser(unittest.TestCase):
    
    def setUp(self):
        self.rec_str="""ID   12OLP_LISIN             Reviewed;        1086 AA.
AC   Q92AT0;
DT   01-APR-2015, integrated into UniProtKB/Swiss-Prot.
DT   01-DEC-2001, sequence version 1.
DT   02-JUN-2021, entry version 90.
DE   RecName: Full=1,2-beta-oligoglucan phosphorylase {ECO:0000305};
DE            EC=2.4.1.333 {ECO:0000269|PubMed:24647662};
GN   OrderedLocusNames=lin1839 {ECO:0000312|EMBL:CAC97070.1};
OS   Listeria innocua serovar 6a (strain ATCC BAA-680 / CLIP 11262).
OC   Bacteria; Firmicutes; Bacilli; Bacillales; Listeriaceae; Listeria.
OX   NCBI_TaxID=272626;
RN   [1]
RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].
RC   STRAIN=ATCC BAA-680 / CLIP 11262 {ECO:0000312|Proteomes:UP000002513};
RX   PubMed=11679669; DOI=10.1126/science.1063447;
RA   Glaser P., Frangeul L., Buchrieser C., Rusniok C., Amend A., Baquero F.,
RA   Berche P., Bloecker H., Brandt P., Chakraborty T., Charbit A.,
RA   Chetouani F., Couve E., de Daruvar A., Dehoux P., Domann E.,
RA   Dominguez-Bernal G., Duchaud E., Durant L., Dussurget O., Entian K.-D.,
RA   Fsihi H., Garcia-del Portillo F., Garrido P., Gautier L., Goebel W.,
RA   Gomez-Lopez N., Hain T., Hauf J., Jackson D., Jones L.-M., Kaerst U.,
RA   Kreft J., Kuhn M., Kunst F., Kurapkat G., Madueno E., Maitournam A.,
RA   Mata Vicente J., Ng E., Nedjari H., Nordsiek G., Novella S., de Pablos B.,
RA   Perez-Diaz J.-C., Purcell R., Remmel B., Rose M., Schlueter T., Simoes N.,
RA   Tierrez A., Vazquez-Boland J.-A., Voss H., Wehland J., Cossart P.;
RT   "Comparative genomics of Listeria species.";
RL   Science 294:849-852(2001).
RN   [2]
RP   FUNCTION, CATALYTIC ACTIVITY, SUBUNIT, AND BIOPHYSICOCHEMICAL PROPERTIES.
RX   PubMed=24647662; DOI=10.1371/journal.pone.0092353;
RA   Nakajima M., Toyoizumi H., Abe K., Nakai H., Taguchi H., Kitaoka M.;
RT   "1,2-beta-Oligoglucan phosphorylase from Listeria innocua.";
RL   PLoS ONE 9:E92353-E92353(2014).
CC   -!- FUNCTION: Catalyzes the reversible phosphorolysis of beta-(1->2)-D-
CC       glucans. The minimum length of the substrate for the phosphorolytic
CC       reaction is 3 D-glucose units. {ECO:0000269|PubMed:24647662}.
CC   -!- CATALYTIC ACTIVITY:
CC       Reaction=[(1->2)-beta-D-glucosyl](n) + phosphate = [(1->2)-beta-D-
CC         glucosyl](n-1) + alpha-D-glucose 1-phosphate; Xref=Rhea:RHEA:47684,
CC         Rhea:RHEA-COMP:11881, Rhea:RHEA-COMP:11882, ChEBI:CHEBI:27517,
CC         ChEBI:CHEBI:43474, ChEBI:CHEBI:58601; EC=2.4.1.333;
CC         Evidence={ECO:0000269|PubMed:24647662};
CC   -!- BIOPHYSICOCHEMICAL PROPERTIES:
CC       Kinetic parameters:
CC         KM=8.5 mM for (beta-(1->2)-D-glucans) x 2 units in the synthetic
CC         reaction {ECO:0000269|PubMed:24647662};
CC         KM=6.0 mM for (beta-(1->2)-D-glucans) x 3 units in the synthetic
CC         reaction {ECO:0000269|PubMed:24647662};
CC         KM=6.8 mM for (beta-(1->2)-D-glucans) x 4 units in the synthetic
CC         reaction {ECO:0000269|PubMed:24647662};
CC         KM=1.2 mM for alpha-D-glucose 1-phosphate in the synthetic reaction
CC         {ECO:0000269|PubMed:24647662};
CC         KM=1.0 mM for (beta-(1->2)-D-glucans) x 3 units in the phosphorolytic
CC         reaction {ECO:0000269|PubMed:24647662};
CC         KM=1.2 mM for (beta-(1->2)-D-glucans) x 4 units in the phosphorolytic
CC         reaction {ECO:0000269|PubMed:24647662};
CC         KM=1.8 mM for (beta-(1->2)-D-glucans) x 5 units in the phosphorolytic
CC         reaction {ECO:0000269|PubMed:24647662};
CC         Note=kcat is 97 sec(-1) for (beta-(1->2)-D-glucans) x 2 units in the
CC         synthetic reaction. kcat is 110 sec(-1) for (beta-(1->2)-D-glucans) x
CC         3 units in the synthetic reaction. kcat is 90 sec(-1) for (beta-
CC         (1->2)-D-glucans) x 4 units in the synthetic reaction. kcat is 43
CC         sec(-1) for alpha-D-glucose 1-phosphate in the synthetic reaction.
CC         kcat is 19 sec(-1) for (beta-(1->2)-D-glucans) x 3 units in the
CC         phosphorolytic reaction. kcat is 30 sec(-1) for (beta-(1->2)-D-
CC         glucans) x 4 units in the phosphorolytic reaction. kcat is 31 sec(-1)
CC         for (beta-(1->2)-D-glucans) x 5 units in the phosphorolytic reaction.
CC         {ECO:0000269|PubMed:24647662};
CC       pH dependence:
CC         Optimum pH is 7.5-8.0. {ECO:0000269|PubMed:24647662};
CC       Temperature dependence:
CC         Optimum temperature is 37-45 degrees Celsius.
CC         {ECO:0000269|PubMed:24647662};
CC   -!- SUBUNIT: Monomer. {ECO:0000305|PubMed:24647662}.
CC   -!- SIMILARITY: Belongs to the glycosyl hydrolase 94 family. {ECO:0000305}.
CC   ---------------------------------------------------------------------------
CC   Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms
CC   Distributed under the Creative Commons Attribution (CC BY 4.0) License
CC   ---------------------------------------------------------------------------
DR   EMBL; AL596170; CAC97070.1; -; Genomic_DNA.
DR   PIR; AF1662; AF1662.
DR   RefSeq; WP_010990982.1; NC_003212.1.
DR   SMR; Q92AT0; -.
DR   STRING; 272626.lin1839; -.
DR   CAZy; GH94; Glycoside Hydrolase Family 94.
DR   EnsemblBacteria; CAC97070; CAC97070; CAC97070.
DR   KEGG; lin:lin1839; -.
DR   eggNOG; COG3459; Bacteria.
DR   HOGENOM; CLU_009079_0_0_9; -.
DR   OMA; VPHGLEQ; -.
DR   OrthoDB; 279804at2; -.
DR   BRENDA; 2.4.1.333; 3044.
DR   Proteomes; UP000002513; Chromosome.
DR   GO; GO:0016740; F:transferase activity; IEA:UniProtKB-KW.
DR   GO; GO:0005975; P:carbohydrate metabolic process; IEA:InterPro.
DR   Gene3D; 1.50.10.10; -; 1.
DR   InterPro; IPR008928; 6-hairpin_glycosidase_sf.
DR   InterPro; IPR012341; 6hp_glycosidase-like_sf.
DR   InterPro; IPR033432; GH36_catalytic.
DR   Pfam; PF17167; Glyco_hydro_36; 1.
DR   SUPFAM; SSF48208; SSF48208; 1.
PE   1: Evidence at protein level;
KW   Transferase.
FT   CHAIN           1..1086
FT                   /note="1,2-beta-oligoglucan phosphorylase"
FT                   /id="PRO_0000432711"
FT   ACT_SITE        738
FT                   /note="Nucleophile"
FT                   /evidence="ECO:0000250"
SQ   SEQUENCE   1086 AA;  122895 MW;  5D87750286582095 CRC64;
     MTMLKEIKKA DLSAAFYPSG ELAWLKLKDI MLNQVIQNPL ENRLSQIYVR AHVGDKIEIY
     PLLSRDAEVG FNENGVEYRG VVGPFRYSVQ MHFHTRGWFY DVTVDGDLEF DLVYLQDLGL
     AEQAAVRTNE AYMSQYIDYH VTEGATGFTV QARQNQPQNE RFPAVQIGAL TKIVGYATDG
     FDIYGTNYKL TSELANLKEK SLPSRVYQYE FAQISLQTEL FTNHGETIFY GYATENQPKA
     SGAPFENLAE LKSNISEQPY QPSTKAILNK HIGTPITGET ISDSWLQENF PDRIQEEQQN
     GALLSFFTPN YAHVVMREKE AELERPHGSI LLDKVDVLNP EATLSATTYM YGAFLSQLVA
     GNTNMNKWNS HARNPLNILQ TSGLRIYIEL DSELRLLGVP SVWETSTNYS TWYYQWNGDL
     ITVQTTLTAD SKEAFVTVHS EKGHSYKLVL TNQVTMGTNE YDTTVKKEIK DGIVTYFPAE
     DSPILETYPA LQFRVDGTYN ELTDERYFAK DYVGTAGLDV FVFEPSDKAT FHVQAKLSDE
     FSKPTEDLEA NNKVIRASYD ELTAQFHLNH QSTTAEKLNL TVYWYAHQML VHYASPHGLE
     QYSGAAWGTR DVSQGPFEFF LATGNKAVLR KLVLTIFSHQ YQDTGDWPQW FMFDKYTSIQ
     QEESHGDVIV WPLKIIGDYL EMSGDAGILE EAIPFVDRES KTFTKEQGTL LEHIELAVKT
     IEARFMKGTA LSNYGDGDWD DTLQPANAQL KKNMVSSWTV ALTYQTFKRL AAFLPVGEKY
     ETLAKNVQAD FAKYMTNDTD VIPGFLYLEE GKAPVWMIHP EDKETNIKYR LIPLTRSVIS
     ELVDKKQASR NFEIIGEHLL HPDGVRLMSE PAHYAGGVST HFKRAEQAAN FGREVGLQYV
     HAHIRYIEAL AKIGDKSAWH MLDVINPINI KEVVPNAALR QSNTYFSSSD AAFLDRYQAQ
     NEFGRVKEGS IPVKGGWRIY SSGPGIYLHQ LISSVLGIRQ TEDALIFDPI LPEELDGLEC
     HIELDNYPLD LTFESADEGS IVVNGEKQPV ENGANLYRTG ALILPKKNLT TKCSQITIKF
     QKNNRL
//
"""
        self.rec=SeqIO.parse(StringIO(self.rec_str), "swiss").__iter__().__next__()
        self.db='uniprot'
        self.seq_type="protein"
        self.rec_type="top"
        self.dbsource_str=""

        ##########
         # What sort of record? Its a protein, described at top-level 
    # (i.e whole record is about this protein))
    def test_determineRecordType(self):
        obs_rec_type=rec_parser.determineRecordType(self.rec)
        exp_rec_type="top"
        self.assertEqual(obs_rec_type, exp_rec_type)        


    # Is it a real protein? Discard Psuedo proteins
    def test_is_pseudo(self):
        obs_pseudo=rec_parser.isPseudo(self.rec, self.rec_type, self.seq_type)
        exp_pseudo=False
        self.assertEqual(obs_pseudo, exp_pseudo)  

    def test_extractName(self):
        obs_name=rec_parser.extractDescription(self.rec.description)
        exp_name='1,2-beta-oligoglucan phosphorylase'
        self.assertEqual(obs_name, exp_name)

    # What organism are we looking at
    def test_extractOrganism(self):
        obs_pseudo=rec_parser.extractOrganism(self.rec.annotations)
        exp_pseudo="Listeria innocua serovar 6a (strain ATCC BAA-680 / CLIP 11262)"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    def test_extractTaxaID(self):
        obs_taxid=rec_parser.extractTaxaID(self.rec.features[0], self.rec.annotations)
        exp_taxid="272626"
        self.assertEqual(obs_taxid, exp_taxid)

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_extractTaxonomy(self):
        obs_pseudo=rec_parser.extractTaxonomy(self.rec.annotations)
        exp_pseudo='Bacteria,Firmicutes,Bacilli,Bacillales,Listeriaceae,Listeria' 
        self.assertEqual(obs_pseudo, exp_pseudo)  

    #Figure out when thi was first submitted
    def test_extractDateModified(self):
        obs=rec_parser.extractDateModified(self.rec, "")
        exp_upload=20210602
        exp_last_mod=20210602
        exp_nmodify=-1
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)   


    def test_extractDateLastModified(self):
        obs_date=rec_parser.extractDateLastModified(self.rec.annotations)
        exp_date=20210602
        self.assertEqual(obs_date, exp_date)  

  
    def test_extractGo(self):
        prot1=""#{f.type:f for f in (r.features[1],r.features[4])}
        obs_go=rec_parser.extractGO(self.rec, prot1, self.db)
        #obs_go=rec_parser.extractGO(self.rec, self.rec_type)
        exp_go=['GO:0016740', 'GO:0005975']
        self.assertEqual(obs_go, exp_go)  

    def test_extractEC(self):
        obs_ec=rec_parser.extractEC(self.rec, "")
        exp_ec='2.4.1.333'
        self.assertEqual(obs_ec, exp_ec)  

    #There are some limitations to the current UNiProt rec_parser. In particular, it seems to only capture a subset of accesssions. 
    def test_extractParentEdges(self):
        obs_edges=rec_parser.extractParentEdges(self.rec,self.db, "")
        #exp_edges=["AY548484.1", "AAT09660.1", "YP_031579.1", "NC_005946.1"]
        exp_edges=[['AL596170', ''], ['WP_010990982', '1']]
        self.assertEqual(sorted(obs_edges), exp_edges)        


    def test_createTopLevelNode(self):
            
            obs=rec_parser.createTopLevelNode(self.rec, self.rec_type, self.seq_type, self.db,self.dbsource_str)
            obs['n']['seq']=obs['n']['seq'][0:10]
            
            exp_node={'id': 'Q92AT0', 
                'seq_version': '',
                'date_first_upload': 20210602,
                'num_modified': -1, 
                'date_last_modified': 20210602,
                'organism': 'Listeria innocua serovar 6a (strain ATCC BAA-680 / CLIP 11262)', 
                'taxa_id': '272626', 
                'taxonomy': 'Bacteria,Firmicutes,Bacilli,Bacillales,Listeriaceae,Listeria',
                'n_products': 0, 
                'name': '1,2-beta-oligoglucan phosphorylase',
                'go': ['GO:0016740', 'GO:0005975'],
                'ec': '2.4.1.333',
                'seq': 'MTMLKEIKKA',
                'db': 'u',
                'seq_type': 'p'

                }

            self.assertEqual(obs['n'], exp_node)

            exp_edges=[('Q92AT0', '', 'AL596170', '', 'p', 'g', 'xref'),
                       ('Q92AT0', '', 'WP_010990982', '1', 'p', 'r', 'xref')]
            self.assertEqual(obs['e'], exp_edges)


    def test_uniprot_parseRecord(self):
        obs_nodes,obs_xref_strs, obs_parent_child_edges,obs_key=rec_parser.parseRecord(self.rec, self.db, self.seq_type)
        exp_xref_strs=[('Q92AT0', '', 'AL596170', '', 'p', 'g', 'xref'),
                  ('Q92AT0', '', 'WP_010990982', '1', 'p', 'r', 'xref')]
        exp_xref_nodes=2
        self.assertEqual(len(obs_xref_strs), exp_xref_nodes)
        self.assertEqual(obs_xref_strs, exp_xref_strs)        