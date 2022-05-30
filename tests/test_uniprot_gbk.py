from sqlite3 import dbapi2
import unittest   # The test framework
from importlib import reload  # Python 3.4+
from netdbqual import parser
from Bio import SeqIO
from io import StringIO

class Test_Parser(unittest.TestCase):
    
    def setUp(self):
      
        self.rec_str="""
LOCUS       PEPT_LISW6               410 aa            linear   BCT 23-FEB-2022
DEFINITION  RecName: Full=Peptidase T; AltName: Full=Aminotripeptidase;
            Short=Tripeptidase; AltName: Full=Tripeptide aminopeptidase.
ACCESSION   A0AJN4
VERSION     A0AJN4.1
DBSOURCE    UniProtKB: locus PEPT_LISW6, accession A0AJN4;
            class: standard.
            created: Jan 15, 2008.
            sequence updated: Nov 28, 2006.
            annotation updated: Feb 23, 2022.
            xrefs: AM263198.1, CAK21216.1, WP_011702575.1
            xrefs (non-sequence databases): SMR:A0AJN4, STRING:386043.lwe1798,
            MEROPS:M20.003, PRIDE:A0AJN4, EnsemblBacteria:CAK21216,
            EnsemblBacteria:CAK21216, EnsemblBacteria:lwe1798, GeneID:61189699,
            KEGG:lwe:lwe1798, eggNOG:COG2195, HOGENOM:CLU_053676_0_0_9,
            OMA:GHNFHGK, OrthoDB:1015417at2, Proteomes:UP000000779, GO:0005737,
            GO:0008237, GO:0045148, GO:0008270, GO:0043171, HAMAP:MF_00550,
            InterPro:IPR001261, InterPro:IPR036264, InterPro:IPR002933,
            InterPro:IPR011650, InterPro:IPR010161, PANTHER:PTHR42994:SF1,
            Pfam:PF07687, Pfam:PF01546, PIRSF:PIRSF037215, SUPFAM:SSF55031,
            TIGRFAMs:TIGR01882, PROSITE:PS00758, PROSITE:PS00759
KEYWORDS    Aminopeptidase; Cytoplasm; Hydrolase; Metal-binding;
            Metalloprotease; Protease; Zinc.
SOURCE      Listeria welshimeri serovar 6b str. SLCC5334
  ORGANISM  Listeria welshimeri serovar 6b str. SLCC5334
            Bacteria; Firmicutes; Bacilli; Bacillales; Listeriaceae; Listeria.
REFERENCE   1  (residues 1 to 410)
  AUTHORS   Hain,T., Steinweg,C., Kuenne,C.T., Billion,A., Ghai,R.,
            Chatterjee,S.S., Domann,E., Karst,U., Goesmann,A., Bekel,T.,
            Bartels,D., Kaiser,O., Meyer,F., Puhler,A., Weisshaar,B.,
            Wehland,J., Liang,C., Dandekar,T., Lampidis,R., Kreft,J., Goebel,W.
            and Chakraborty,T.
  TITLE     Whole-genome sequence of Listeria welshimeri reveals common steps
            in genome reduction with Listeria innocua as compared to Listeria
            monocytogenes
  JOURNAL   J Bacteriol 188 (21), 7405-7415 (2006)
   PUBMED   16936040
  REMARK    NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].;
            STRAIN=ATCC 35897 / DSM 20650 / CIP 8149 / NCTC 11857 / SLCC 5334 /
            V8
COMMENT     [FUNCTION] Cleaves the N-terminal amino acid of tripeptides.
            {ECO:0000255|HAMAP-Rule:MF_00550}.
            [CATALYTIC ACTIVITY] Reaction=Release of the N-terminal residue
            from a tripeptide.; EC=3.4.11.4;
            Evidence={ECO:0000255|HAMAP-Rule:MF_00550}.
            [COFACTOR] Name=Zn(2+); Xref=ChEBI:CHEBI:29105;
            Evidence={ECO:0000255|HAMAP-Rule:MF_00550}; Note=Binds 2 Zn(2+)
            ions per subunit. {ECO:0000255|HAMAP-Rule:MF_00550}.
            [SUBCELLULAR LOCATION] Cytoplasm {ECO:0000255|HAMAP-Rule:MF_00550}.
            [SIMILARITY] Belongs to the peptidase M20B family.
            {ECO:0000255|HAMAP-Rule:MF_00550}.
FEATURES             Location/Qualifiers
     source          1..410
                     /organism="Listeria welshimeri serovar 6b str. SLCC5334"
                     /strain="SLCC5334"
                     /serovar="6b"
                     /type_material="type strain of Listeria welshimeri"
                     /db_xref="taxon:386043"
     gene            1..410
                     /gene="pepT"
                     /locus_tag="lwe1798"
     Protein         1..410
                     /product="Peptidase T"
                     /EC_number="3.4.11.4"
                     /note="Aminotripeptidase; Tripeptide aminopeptidase;
                     Tripeptidase"
                     /UniProtKB_evidence="Inferred from homology"
     Region          1..410
                     /region_name="Mature chain"
                     /note="Peptidase T. /id=PRO_1000017849."
     Region          2..409
                     /region_name="PRK05469"
                     /note="tripeptide aminopeptidase PepT"
                     /db_xref="CDD:235484"
     Site            79
                     /site_type="metal-binding"
                     /note="Zinc 1. /evidence=ECO:0000255|HAMAP-Rule:MF_00550."
     Site            81
                     /site_type="active"
                     /note="/evidence=ECO:0000255|HAMAP-Rule:MF_00550."
     Site            142
                     /site_type="metal-binding"
                     /note="Zinc 1. /evidence=ECO:0000255|HAMAP-Rule:MF_00550."
     Site            142
                     /site_type="metal-binding"
                     /note="Zinc 2. /evidence=ECO:0000255|HAMAP-Rule:MF_00550."
     Site            176
                     /site_type="active"
                     /note="Proton acceptor.
                     /evidence=ECO:0000255|HAMAP-Rule:MF_00550."
     Site            177
                     /site_type="metal-binding"
                     /note="Zinc 2. /evidence=ECO:0000255|HAMAP-Rule:MF_00550."
     Site            199
                     /site_type="metal-binding"
                     /note="Zinc 1. /evidence=ECO:0000255|HAMAP-Rule:MF_00550."
     Site            381
                     /site_type="metal-binding"
                     /note="Zinc 2. /evidence=ECO:0000255|HAMAP-Rule:MF_00550."
ORIGIN
        1 mkeellkrft kyvkvdtqsn eestvcpttk gqmelanilv telkeigmqe vtvdefgyvm
       61 atlpsnttke vpvigflahl dtatdltgkn vqpqvhenyd gkdivlnkel nvvlspkqfp
      121 elaqykgktl ittdgttllg addkagitei mvamnylinh peikhgkiri aftpdeeigr
      181 gperfdveaf gakyaytmdg gplgeleyes fnaaaakitf ngnsvhpgta knkmvnavkm
      241 amefdaripk eeapehtdgy egfyhlisln gdveqaksyy iirdfdhlkf verkthiatv
      301 akeleekygk gsvelklndq yynmrekiep vkeivdivsa amrnldiepk ispirggtdg
      361 aqlsykglpt pnifgggenf hgkfeyvale smvkatevii evarlfeeke
//
"""
        self.rec=SeqIO.parse(StringIO(self.rec_str), "gb").__iter__().__next__()
        self.db='uniprot'
        self.seq_type="protein"
        self.rec_type="top"
        self.dbsource_str=parser.processUniProtDBsource(self.rec.annotations['db_source'])


      
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


    # What organism are we looking at
    def test_extractOrganism(self):
        obs_pseudo=parser.extractOrganism(self.rec.annotations)
        exp_pseudo="Listeria welshimeri serovar 6b str. SLCC5334"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_extractTaxonomy(self):
        obs_pseudo=parser.extractTaxonomy(self.rec.annotations)
        exp_pseudo='Bacteria,Firmicutes,Bacilli,Bacillales,Listeriaceae,Listeria' 
        self.assertEqual(obs_pseudo, exp_pseudo)  

    #Figure out when thi was first submitted
    def test_extractDateModified(self):
        db_str=parser.processUniProtDBsource(self.rec.annotations['db_source'])
        obs=parser.extractDateModified(self.rec, db_str)
        exp_upload=20080115
        exp_last_mod=20220223
        exp_nmodify=-1
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)   


    def test_extractDateLastModified(self):
        obs_date=parser.extractDateLastModified(self.rec.annotations)
        exp_date=20220223
        self.assertEqual(obs_date, exp_date)  

  
    def test_extractGo(self):
        r=self.rec
        prot1={f.type:f for f in (r.features[1],r.features[4])}
        obs_go=parser.extractGO(self.rec, prot1, self.db)
        #obs_go=parser.extractGO(self.rec, self.rec_type)
        exp_go=['GO:0005737', 'GO:0008237', 'GO:0008270', 'GO:0043171', 'GO:0045148']
        self.assertEqual(obs_go, exp_go)  

    def test_extractEC(self):
        r=self.rec
        prot1={f.type:f for f in (r.features[2],r.features[4])}
        
        obs_ec=parser.extractEC(r, prot1['Protein'])
        #obs_ec=parser.extractEC(self.rec, self.rec_type)
        exp_ec="3.4.11.4"
        self.assertEqual(obs_ec, exp_ec)  

    def test_extractTaxaID(self):
        obs_taxid=parser.extractTaxaID(self.rec.features[0])
        exp_taxid="386043"
        self.assertEqual(obs_taxid, exp_taxid)

    def test_extractParentEdges(self):
        db_str=parser.processUniProtDBsource(self.rec.annotations['db_source'])
        r=self.rec
        obs_edges=parser.extractParentEdges(r,self.db, db_str)
        exp_edges=['AM263198.1', 'CAK21216.1', 'WP_011702575.1']
        self.assertEqual(obs_edges, exp_edges)        


    def test_createTopLevelNode(self):
            r=self.rec  
            rec_type=self.rec_type
            seq_type=self.seq_type
            db=self.db
            dbsrc_str=self.dbsource_str
            obs=parser.createTopLevelNode(r, rec_type, seq_type, db,dbsrc_str)
            obs['n']['seq']=obs['n']['seq'][0:10]
            
            exp_node={'id': 'A0AJN4.1',
             'seq_version': 1, 
             'date_first_upload': 20080115, 
             'num_modified': -1, 
             'date_last_modified': 20220223, 
             'organism': 'Listeria welshimeri serovar 6b str. SLCC5334',
             'taxa_id': '386043', 
             'taxonomy': 'Bacteria,Firmicutes,Bacilli,Bacillales,Listeriaceae,Listeria', 
             'n_products': 1, 
             'name': 'Peptidase T', 
             'go': ['GO:0005737', 'GO:0008237', 'GO:0008270', 'GO:0043171', 'GO:0045148'], 
             'ec': '3.4.11.4', 
             'seq': 'MKEELLKRFT'#, 
             #'parent': ['AM263198.1', 
             ##'CAK21216.1', 
             #'WP_011702575.1']
             }

            self.assertEqual(obs['n'], exp_node)

            exp_edges=[('A0AJN4.1','AM263198.1'), 
             ('A0AJN4.1','CAK21216.1'), 
             ('A0AJN4.1','WP_011702575.1')]
            self.assertEqual(obs['e'], exp_edges)


    def test_uniprotDateStrToNum(self):
        date_str='Jan 15, 2008'
        obs=parser.uniprotDateStrToNum(date_str)
        exp=20080115
        self.assertEqual(obs, exp)

    def test_uniprot_parseRecord(self):
        r=self.rec  
        db=self.db
        seq_type=self.seq_type
        obs_node, obs_edge=parser.parseRecord(r, db, seq_type)
        exp_edge=[('A0AJN4.1', 'AM263198.1'), ('A0AJN4.1', 'CAK21216.1'), ('A0AJN4.1', 'WP_011702575.1')]
        exp_nodes=1
        self.assertEqual(len(obs_node), exp_nodes)
        self.assertEqual(obs_edge, exp_edge)



