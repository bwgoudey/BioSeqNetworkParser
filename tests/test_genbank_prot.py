
from sqlite3 import dbapi2
import unittest   # The test framework
from importlib import reload  # Python 3.4+
from netdbqual import rec_parser
from Bio import SeqIO
from io import StringIO

class Test_Parser(unittest.TestCase):
    
    def setUp(self):
        self.rec_str="""
LOCUS       UOH70005                 412 aa            linear   BCT 07-APR-2022
DEFINITION  peptidase T [Staphylococcus epidermidis].
ACCESSION   UOH70005
VERSION     UOH70005.1
DBLINK      BioProject,PRJNA821968
            BioSample: SAMN27161489
DBSOURCE    accession CP094865.1
KEYWORDS    .
SOURCE      Staphylococcus epidermidis
  ORGANISM  Staphylococcus epidermidis
            Bacteria; Firmicutes; Bacilli; Bacillales; Staphylococcaceae;
            Staphylococcus.
REFERENCE   1  (residues 1 to 412)
  AUTHORS   Kopczyk,M., Mulroney,K., Chakera,A., Ramsay,J.P., Colombi,E.,
            Coombs,G.W. and Mowlaboccus,S.
  TITLE     Sequencing of Staphylococcus epidermidis isolated from peritoneal
            cavity
  JOURNAL   Unpublished
REFERENCE   2  (residues 1 to 412)
  AUTHORS   Kopczyk,M., Mulroney,K., Chakera,A., Ramsay,J.P., Colombi,E.,
            Coombs,G.W. and Mowlaboccus,S.
  TITLE     Direct Submission
  JOURNAL   Submitted (31-MAR-2022) School of Pharmacy & Biomedical Sciences,
            Curtin University, Kent street, Bentley, WA 6102, Australia
COMMENT     The annotation was added by the NCBI Prokaryotic Genome Annotation
            Pipeline (PGAP). Information about PGAP can be found here:
            https://www.ncbi.nlm.nih.gov/genome/annotation_prok/
            
            ##Genome-Annotation-Data-START##
            Annotation Provider               :: NCBI
            Annotation Date                   :: 04/01/2022 12:38:54
            Annotation Pipeline               :: NCBI Prokaryotic Genome
                                                 Annotation Pipeline (PGAP)
            Annotation Method                 :: Best-placed reference protein
                                                 set; GeneMarkS-2+
            Annotation Software revision      :: 6.1
            Features Annotated                :: Gene; CDS; rRNA; tRNA; ncRNA;
                                                 repeat_region
            Genes (total)                     :: 2,411
            CDSs (total)                      :: 2,327
            Genes (coding)                    :: 2,267
            CDSs (with protein)               :: 2,267
            Genes (RNA)                       :: 84
            rRNAs                             :: 7, 6, 6 (5S, 16S, 23S)
            complete rRNAs                    :: 7, 6, 6 (5S, 16S, 23S)
            tRNAs                             :: 61
            ncRNAs                            :: 4
            Pseudo Genes (total)              :: 60
            CDSs (without protein)            :: 60
            Pseudo Genes (ambiguous residues) :: 0 of 60
            Pseudo Genes (frameshifted)       :: 34 of 60
            Pseudo Genes (incomplete)         :: 33 of 60
            Pseudo Genes (internal stop)      :: 20 of 60
            Pseudo Genes (multiple problems)  :: 21 of 60
            ##Genome-Annotation-Data-END##
FEATURES             Location/Qualifiers
     source          1..412
                     /organism="Staphylococcus epidermidis"
                     /strain="C100"
                     /isolation_source="pertoneal cavity"
                     /host="Homo sapiens"
                     /db_xref="taxon:1282"
                     /country="Australia: Perth"
                     /lat_lon="31.92851574 S 115.86399442 E"
                     /collection_date="2017-04-19"
                     /collected_by="Harry Perkins Institute of Medical
                     Research"
     Protein         1..412
                     /product="peptidase T"
                     /EC_number="3.4.11.4"
                     /note="GO_function: GO:0045148 - tripeptide aminopeptidase
                     activity [Evidence IEA];
                     GO_process: GO:0006508 - proteolysis [Evidence IEA]"
     CDS             1..412
                     /gene="pepT"
                     /locus_tag="MUG64_09775"
                     /coded_by="CP094865.1:2064704..2065942"
                     /inference="COORDINATES: similar to AA
                     sequence:RefSeq:WP_002468678.1"
                     /note="Derived by automated computational analysis using
                     gene prediction method: Protein Homology.
                     GO_function: GO:0045148 - tripeptide aminopeptidase
                     activity [Evidence IEA];
                     GO_process: GO:0006508 - proteolysis [Evidence IEA]"
                     /transl_table=11
ORIGIN      
        1 mkkqiierlt ryvkidtqsn pdskttpstn kqwdllnlle eelqslglkt dmdehgylfa
       61 tlesninynv ptvgflahvd tspdfnashv npqiieayng qpiklgesqr ildpdvfpel
      121 nkvvghtlmv tdgtsllgad dkagvveime gikylidhpd ikhgtirvgf tpdeeigrgp
      181 hqfdvsrfna dfaytmdgsq lgelqfesfn aaevtvtchg vnvhpgsakn amvnaislgq
      241 qfnsllpshe vpertegyeg fyhlmnftgn vekatlqyii rdhdkeqfel rkkrmmeird
      301 dinvhynhfp ikvdvhdqyf nmaekieplk hiidipkrvf ealdivpnte pirggtdgsq
      361 lsfmglptpn iftgcgnfhg pfeyasidvm ekavhvvvgi aqevanshqs yk
//
"""
        self.rec=SeqIO.parse(StringIO(self.rec_str), "gb").__iter__().__next__()
        self.db='genbank'
        self.seq_type="protein"
        self.rec_type="top"


    # Is it a real protein? Discard Psuedo proteins
    def test_is_pseudo(self):
        obs_pseudo=rec_parser.isPseudo(self.rec, self.rec_type, self.seq_type)
        exp_pseudo=False
        self.assertEqual(obs_pseudo, exp_pseudo)  


    # What organism are we looking at
    def test_extractOrganism(self):
        obs_pseudo=rec_parser.extractOrganism(self.rec.annotations)
        exp_pseudo="Staphylococcus epidermidis"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_extractTaxonomy(self):
        obs_pseudo=rec_parser.extractTaxonomy(self.rec.annotations)
        exp_pseudo='Bacteria,Firmicutes,Bacilli,Bacillales,Staphylococcaceae,Staphylococcus'
        
        obs_pseudo_short=rec_parser.extractTaxonomy({'taxonomy':['Bacteria', 'Firmicutes', 'Bacilli']})
        exp_pseudo_short='Bacteria,Firmicutes,Bacilli,,,'
        
        self.assertEqual(obs_pseudo, exp_pseudo)  
        self.assertEqual(obs_pseudo_short, exp_pseudo_short)  

    #Figure out when thi was first submitted
    def test_extractDateModified(self):
        obs=rec_parser.extractDateModified(self.rec, "")
        exp_upload=20220331
        exp_last_mod=20220331
        exp_nmodify=1
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)   


    def test_extractDateLastModified(self):
        obs_date=rec_parser.extractDateLastModified(self.rec.annotations)
        exp_date=20220407
        self.assertEqual(obs_date, exp_date)  

    def test_extractGo(self):
       
        r=self.rec
        prot1={f.type:f for f in (r.features[1],r.features[2])}
        obs_go=rec_parser.extractGO(self.rec, prot1, self.db)
        exp_go=sorted(['GO:0006508', 'GO:0045148'])
        self.assertEqual(obs_go, exp_go)  

    def test_extractEC(self):
        r=self.rec
        prot1={f.type:f for f in (r.features[1],r.features[2])}
        obs_ec=rec_parser.extractEC(r, prot1['Protein'])
        exp_ec="3.4.11.4"
        self.assertEqual(obs_ec, exp_ec)  

    def test_extractTaxaID(self):
        obs_taxid=rec_parser.extractTaxaID(self.rec.features[0])
        exp_taxid="1282"
        self.assertEqual(obs_taxid, exp_taxid)

    def test_extractParentEdges(self):
        obs_edges=rec_parser.extractParentEdges(self.rec, self.db)
        exp_edges=[['CP094865', '1']]
        self.assertEqual(obs_edges, exp_edges)

    def test_extractInferenceEdges(self):
        obs_edges=rec_parser.extractInferenceEdges(self.rec, self.db)
        exp_edges=[['CP094865', 1, "g", 'WP_002468678', 1, "r", "Protein Homology"]]
        self.assertEqual(obs_edges, exp_edges)

    # def test_extractEdges(self):
    #     obs_edges=rec_parser.extractEdges(self.rec, self.db)
    #     exp_edges=[(self.rec.id, 'CP094865.1')]
    #     self.assertEqual(obs_edges, exp_edges)

    def test_top_createTopLevelNode(self):
        r=self.rec  
        rec_type=self.rec_type
        seq_type=self.seq_type
        db=self.db
        obs=rec_parser.createTopLevelNode(r, rec_type, seq_type, db)
        obs['n']['seq']=obs['n']['seq'][0:10]
        
        exp_node={'id': 'UOH70005', 
        'seq_version': 1, 
        'date_first_upload': 20220331, 
        'num_modified': 1, 
        'date_last_modified': 20220407,
        'organism': 'Staphylococcus epidermidis', 
        'taxa_id': '1282', 
        'taxonomy': 'Bacteria,Firmicutes,Bacilli,Bacillales,Staphylococcaceae,Staphylococcus', 
        'n_products': 1, 
        'name': 'peptidase T', 
        'go': ['GO:0006508', 'GO:0045148'], 
        'ec': '3.4.11.4', 
        'seq': 'MKKQIIERLT',
        #'parent': ['CP094865.1']
        'db': 'g',
        'seq_type': 'p'
        }

        self.assertEqual(obs['n'], exp_node)
        self.assertEqual(obs['e'], [('UOH70005', '1', 'CP094865', '1', 'p', 'g', 'xref')])

