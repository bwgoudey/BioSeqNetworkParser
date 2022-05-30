from sqlite3 import dbapi2
import unittest   # The test framework
from importlib import reload  # Python 3.4+
from netdbqual import parser
from Bio import SeqIO
from io import StringIO

class Test_Parser(unittest.TestCase):
    
    def setUp(self):
       
        self.rec_str="""
LOCUS       NC_005946             4861 bp    DNA     linear   VRL 20-DEC-2020
DEFINITION  Frog virus 3, complete genome.
ACCESSION   NC_005946 NC_012637
VERSION     NC_005946.1
DBLINK      BioProject: PRJNA485481
KEYWORDS    RefSeq.
SOURCE      Frog virus 3
  ORGANISM  Frog virus 3
            Viruses; Varidnaviria; Bamfordvirae; Nucleocytoviricota;
            Megaviricetes; Pimascovirales; Iridoviridae; Alphairidovirinae;
            Ranavirus.
REFERENCE   1  (bases 1 to 4861)
  AUTHORS   Tan,W.G., Barkman,T.J., Gregory Chinchar,V. and Essani,K.
  TITLE     Comparative genomic analyses of frog virus 3, type species of the
            genus Ranavirus (family Iridoviridae)
  JOURNAL   Virology 323 (1), 70-84 (2004)
   PUBMED   15165820
REFERENCE   2  (bases 1 to 4861)
  CONSRTM   NCBI Genome Project
  TITLE     Direct Submission
  JOURNAL   Submitted (18-AUG-2004) National Center for Biotechnology
            Information, NIH, Bethesda, MD 20894, USA
REFERENCE   3  (bases 1 to 4861)
  AUTHORS   Tan,W.G.H., Barkman,T.J., Chinchar,V.G. and Essani,K.
  TITLE     Direct Submission
  JOURNAL   Submitted (18-FEB-2004) Department of Biological Sciences, Western
            Michigan University, Room 3441, Wood Hall, 1903 West Michigan Ave,
            Kalamazoo, MI 49008, USA
COMMENT     PROVISIONAL REFSEQ: This record has not yet been subject to final
            NCBI review. The reference sequence was derived from AY548484.
            On Jul 12, 2012 this sequence version replaced NC_012637.1.
            COMPLETENESS: full length.
FEATURES             Location/Qualifiers
     source          1..4601
                     /organism="Frog virus 3"
                     /mol_type="genomic DNA"
                     /host="Rana pipiens"
                     /db_xref="taxon:10493"
     gene            272..1042
                     /locus_tag="FV3gorf1R"
                     /db_xref="GeneID:2947773"
     CDS             272..1042
                     /locus_tag="FV3gorf1R"
                     /note="putative DNA packing protein; orf1R"
                     /codon_start=1
                     /product="putative replicating factor"
                     /protein_id="YP_031579.1"
                     /db_xref="GeneID:2947773"
                     /translation="MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRV
                     QVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRR
                     YQGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQY
                     ENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHG
                     HGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL"
     gene            complement(1649..2611)
                     /locus_tag="FV3gorf2L"
                     /db_xref="GeneID:2947774"
     CDS             complement(1649..2611)
                     /locus_tag="FV3gorf2L"
                     /note="orf2L"
                     /codon_start=1
                     /product="putative myristylated membrane protein"
                     /protein_id="YP_031580.1"
                     /db_xref="GeneID:2947774"
                     /translation="MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGE
                     QTCASGFCTSQPLCARIKKTQVCGLRYSSKGKDPLVSAEWDSRGAPYVRCTYDADLID
                     TQAQVDQFVSMFGESPSLAERYCMRGVKNTAGELVSRVSSDADPAGGWCRKWYSAHRG
                     PDQDAALGSFCIKNPGAADCKCINRASDPVYQKVKTLHAYPDQCWYVPCAADVGELKM
                     GTQRDTPTNCPTQVCQIVFNMLDDGSVTMDDVKNTINCDFSKYVPPPPPPKPTPPTPP
                     TPPTPPTPPTPPTPPTPRPVHNRKVMFFVAGAVLVAILISTVRW"
     gene            3418..4734
                     /locus_tag="FV3gorf3R"
                     /db_xref="GeneID:2947775"
     CDS             3418..4734
                     /locus_tag="FV3gorf3R"
                     /note="similar to IIV6 orf229L protein; orf3R"
                     /codon_start=1
                     /product="orf229L protein-like protein"
                     /protein_id="YP_031581.1"
                     /db_xref="GeneID:2947775"
                     /translation="MARPLLGKTSSVRRRLESLSACSIFFFLRKFCQKMASLVFLNSP
                     VYQMSNILLTERRQVDRAMGGSDDDGVMVVALSPSDFKTVLGSALLAVERDMVHVVPK
                     YLQTPGILHDMLVLLTPIFGEALSVDMSGATDVMVQQIATAGFVDVDPLHSSVSWKDN
                     VSCPVALLAVSNAVRTMMGQPCQVTLIIDVGTQNILRDLVNLPVEMSGDLQVMAYTKD
                     PLGKVPAVGVSVFDSGSVQKGDAHSVGAPDGLVSFHTHPVSSAVELNYHAGWPSNVDM
                     SSLLTMKNLMHVVVAEEGLWTMARTLSMQRLTKVLTDAEKDVMRAAAFNLFLPLNELR
                     VMGTKDSNNKSLKTYFEVFETFTIGALMKHSGVTPTAFVDRRWLDNTIYHMGFIPWGR
                     DMRFVVEYDLDGTNPFLNTVPTLMSVKRKAKIQEMFDNMVSRMVTS"
ORIGIN      
        1 aagctttaac agattcatga aattgtattt tattaaggga taactggtaa cccgagaaga
       61 cgtcaaagaa ctttgacagc aatgcgcaaa aaagagagat taagactaat ctctctcaat
      121 actaagaatg tctttggact aaaaatgtgg tacactacat ttctagtcta tctctgtcag
      181 aaaatcttaa gattctatcc ctttcagaag atattatctt gatacaccga taacctcaag
      241 acctagttct gcgacgctac ggtaactata aatggcattc tcggcagaag atgtgctgaa
      301 ggagtacgac aggagacgga ggatggaggc cctcttgctc agcctgtact acccaaacga
      361 ccgcaagctc ctagactaca aagagtggtc tccgcccagg gttcaggtag agtgtcccaa
      421 agcccccgtg gagtggaaca accctccgtc agaaaagggt ctcatcgtgg ggcactttag
      481 cggcataaag tacaaggggg aaaaggctca ggcatccgag gtagacgtca acaagatgtg
      541 ctgctgggtg tccaagttta aagacgccat gaggaggtac cagggcatac agacttgcaa
      601 gatccccggc aaggtcctgt cggacctcga cgccaaaata aaggcttaca acctcaccgt
      661 tgagggcgta gagggtttcg tgaggtactc acgagtgacc aagcagcacg tagcagcttt
      721 cctcaaggag ctcaggcact ctaagcagta cgaaaacgtc aacctcatcc actacatcct
      781 caccgacaag agggtagaca ttcagcacct ggaaaaggat cttgtcaagg attttaaggc
      841 gctggtggaa tctgctcaca ggatgaggca gggccacatg atcaacgtaa agtacatact
      901 ctaccagctc ctcaagaagc acggtcacgg gccagacggt ccagacatcc tgaccgtaaa
      961 gactggaagc aagggagtct tgtacgacga ttcctttcgc aagatttaca cggacctcgg
      1021 gtggaagttt acccccctat gaaggtcggg gattgagagt atatcctaca agagtatctc
      1081 atacaaacac gtatccgaca agcgtctata gctcttgtcg gatctgagat gctaacgcgc
      1141 tgcaactgtt ttaggacaag aatgaggtgt aacatattaa aatactttac tatctcaaga
      1201 acctttatcg taagatactt tcaggagaag aaggagatgt aacatcttaa gatatagctt
      1261 ctctcagaaa gtttctagga caagaaggag atgtaacaac atttatagtc ttaagatata
      1321 gcttctctca gaaagtttat aggatgagaa ggagatgtaa caacatttat agaatatctt
      1381 acgatagtaa tcttaagata ctttactact atcgtaagat actctctcaa gattaagact
      1441 gtcttttggt ctaaaaacgt ggtacaacat ttatagtctt ttacgctctc agaaactttc
      1501 taggacaaga cggagatgta atatcttagg ataaagtctt aagagatctt tactctttca
      1561 gaagactttc tctcagaaga ctctcgacat tatatctcta gatcttgtaa attcttgtcc
      1621 tactcacagg accagaattt ataaccgctt accatctcac tgtagagatc agtatagcca
      1681 ccagcactgc tccggccaca aagaacatta ccttcctgtt gtgcaccggc ctgggagtgg
      1741 ggggagtggg gggagtggga ggagtgggag gagtgggagg agtgggagga gtgggaggag
      1801 tgggtttggg aggcggcgga ggagggacgt acttggaaaa gtcacagttt atagtgttct
      1861 tgacgtcgtc catggtcacg ctcccgtcgt cgagcatgtt aaacactatc tggcacactt
      1921 gggtaggaca gtttgtcggc gtgtccctct gggtccccat ctttagctct ccgacgtcag
      1981 ccgcgcaggg cacgtaccag cactggtcgg ggtaggcgtg gagagtcttt accttctggt
      2041 agacgggatc gctggccctg tttatgcact tgcagtctgc cgctccgggg ttcttgatac
      2101 agaaggagcc cagggcggcg tcctggtcgg ggcccctgtg ggccgagtac cactttctgc
      2161 accagccccc ggcggggtcc gcgtccgagg acaccctgct cacgagctcc ccggccgtgt
      2221 tctttacccc cctcatgcag tacctctcgg ccaggctggg agactcgccg aacatggaca
      2281 caaactggtc cacctgggcc tgggtgtcta tcaggtccgc gtcgtatgtg cacctcacgt
      2341 aaggggcccc cctcgagtcc cactcggcgc tcaccagagg gtccttaccc ttggaagagt
      2401 acctgaggcc gcacacttgg gtctttttta tgcgggcgca caggggctga gaggtgcaga
      2461 agcccgaggc gcacgtctgc tcccccaggt cccagtcctt gccgcaggtg cccctgggtg
      2521 tgaaggcgct gcagcccccg gcgtagcagg ggccggccga gtacgtgtcg gacttgtcgt
      2581 tctgtaatct ggtcgctccg atgatggaca ttttaatagc ctttttctct cacgagggag
      2641 ggtaaatttt aacagtctgg ggtccaggac acaatctctg ctagagcgca accgcgggga
      2701 gggttcctga gctctggaaa gtacctggcg cacacggcgt cgggcctagt ggcccagccg
      2761 gcgtcctgcg aggacaggtc ggcggcctct atatccgagg cgatcctcgc tgcggtcttt
      2821 gacagaggcg agaaggccct ggcgacctcc gcggtgcagt ctgcaaactc cttgacgggg
      2881 aacatcctgg gcatcctcct gtccctcctg gaggccggga cgacgtcggt gtctgagact
      2941 ccagcgtccc ctctccagcc acaggccccc ctcttggtcc tcacagcgac acccacagga
      3001 gtgcccctct cgccgtggtc cgggtcgatc accgcgtcgt acacgacagt gtcgtggaag
      3061 ccgtcgagca tgacgtgcgg ggcgtcgcac ctcaggtgga actctctgtc ccccatcctg
      3121 tagagcctga aagagttgcg gtgcctgggg tcgggagccg ccagcctcac tgccgatggc
      3181 ctgagcctgg agtgcgccat ggccatctcg ccgtgcatga cggctgtggc catgagggcc
      3241 tgctggaggc agttgtccac cgtctcctca gaaggggaaa agtgctccag ccagtcgccg
      3301 agggtgtagt cgtagtgctc gtgaggaccg gcgtagaccc ctacacccgc cccggtgttt
      3361 aggaggacgt ccgaggcgta cctcctgtgg tcctcgtgct tggacctgtg cctgtccatg
      3421 gccaggcccc tcttgggaaa gacgtcgtct gtgcggcggc gtctggaaag tctgtcggcg
      3481 tgcagcatat tctttttttt aagaaagttc tgtcagaaaa tggcatcgct agtgtttcta
      3541 aacagccccg tctaccagat gagtaacatc ctcctcacgg agaggaggca ggtggacagg
      3601 gcgatgggag ggtcggacga cgacggcgtc atggtcgtgg ccctgtcccc ctcggacttt
      3661 aagactgtcc tggggtctgc cctcctcgcc gtggagaggg acatggtaca cgtcgtacca
      3721 aagtacctgc agacgccggg gatcctccac gacatgctcg tgttgcttac ccccatattc
      3781 ggcgaggccc tgtcggttga catgagtggc gccacggacg tcatggttca gcagattgcc
      3841 acggcggggt tcgtagacgt agacccgctt cactcttctg tgtcgtggaa ggacaacgtg
      3901 tcttgcccgg tggcgttgct ggccgtctcc aacgccgtca ggaccatgat gggacagccc
      3961 tgtcaggtca ccctcatcat agacgtcggc acccagaaca tcctgaggga cctggtgaac
      4021 ctcccagtgg agatgtccgg agacctccag gttatggcat acaccaagga ccccctgggg
      4081 aaggttccag ccgtaggagt gtccgtgttt gacagcgggt cggtgcagaa gggggacgct
      4141 cactccgtgg gagcccccga cggcctcgtg tccttccaca cccaccctgt gtcgtccgcg
      4201 gtcgaactaa actaccacgc ggggtggccc tccaacgtag acatgtcctc cctgctcacc
      4261 atgaagaacc tcatgcacgt ggttgtcgcg gaggagggcc tgtggactat ggcgaggacc
      4321 ctgtccatgc agaggctcac caaggtcctc acggatgcgg aaaaggacgt catgagggcc
      4381 gcggccttta acctgtttct ccccctcaac gaactcaggg tgatggggac caaggacagc
      4441 aacaacaagt ccctcaagac ttactttgag gtgtttgaga cgtttaccat aggagccctg
      4501 atgaagcact ctggcgtgac ccccaccgcc tttgtggaca ggaggtggct ggacaacacc
      4561 atttaccaca tgggcttcat cccgtggggc agagacatga ggttcgtggt agagtacgac
      4621 ctggacggca ccaacccatt cctgaacact gtccccaccc tcatgtccgt caagagaaag
      4681 gccaaaatcc aggagatgtt tgacaacatg gtgagcagga tggtcacctc ataataattt
      4741 tttcccacac tctggggaaa aaaaaccagt cgagatgaac gcaaaatacg acacagatca
      4801 gggcgtcggt cgcatgcttt tccttggtac gatcggcctc gccgtagttg tcggaggcct
//
    """
        self.rec=SeqIO.parse(StringIO(self.rec_str), "gb").__iter__().__next__()
        self.db='refseq'
        self.seq_type="nucleotide"
        self.rec_type="top"



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
        exp_pseudo="Frog virus 3"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_extractTaxonomy(self):
        obs_pseudo=parser.extractTaxonomy(self.rec.annotations)
        exp_pseudo='Viruses,Varidnaviria,Bamfordvirae,Nucleocytoviricota,Megaviricetes,Pimascovirales' 
        self.assertEqual(obs_pseudo, exp_pseudo)  

    #Figure out when thi was first submitted
    def test_extractDateModified(self):
        obs=parser.extractDateModified(self.rec, "")
        exp_upload=20040218#18-FEB-2004
        exp_last_mod=20040818#18-AUG-2004
        exp_nmodify=2
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)  


    def test_extractDateLastModified(self):
        obs_date=parser.extractDateLastModified(self.rec.annotations)
        exp_date=20201220
        self.assertEqual(obs_date, exp_date)  
    
    def test_extractGo(self):
        r=self.rec
        prot1={f.type:f for f in (r.features[1],r.features[2])}
        obs_go=parser.extractGO(r, prot1, self.db)

        exp_go=[]
        self.assertEqual(obs_go, exp_go)  

    def test_extractEC(self):
        obs_ec=parser.extractEC(self.rec, self.rec.features[1])
        exp_ec=""
        self.assertEqual(obs_ec, exp_ec)  

    def test_extractTaxaID(self):
        obs_taxid=parser.extractTaxaID(self.rec.features[0])
        exp_taxid="10493"
        self.assertEqual(obs_taxid, exp_taxid)

    def test_extractParentEdges(self):
        r=self.rec        
        obs_edges=parser.extractParentEdges(r, self.db)
        exp_edges=['AY548484'] 
        self.assertEqual(obs_edges, exp_edges)   

    # def test_extractEdges(self):
    #     r=self.rec
    #     fd={f.type:f for f in (r.features[1],r.features[2])}
    #     obs_edges=parser.extractRelatedEdges(fd)
    #     exp_edges={"YP_031579.1":[]}
    #     self.assertEqual(obs_edges, exp_edges)    

    def test_identifyProteins(self):
        obs_prots=parser.identifyProteins(self.rec)
        obs_nprots=len(obs_prots)
        exp_nprots=3
        self.assertEqual(obs_nprots, exp_nprots)

    def test_createTopLevelNode(self):
        r=self.rec  
        rec_type=self.rec_type
        seq_type=self.seq_type
        db=self.db
        obs=parser.createTopLevelNode(r, rec_type, seq_type, db)
        obs['n']['seq']=obs['n']['seq'][0:10]
        
        exp_node={'id': 'NC_005946.1', 
        'seq_version': 1, 
        'date_first_upload': 20040218, 
        'num_modified': 2, 
        'date_last_modified': 20201220, 
        'organism': 'Frog virus 3', 
        'taxa_id': '10493',
        'taxonomy': 'Viruses,Varidnaviria,Bamfordvirae,Nucleocytoviricota,Megaviricetes,Pimascovirales', 
        'n_products': 3, 
        'name': 'Frog virus 3, complete genome', 
        'go': '', 
        'ec': '', 
        'seq': "AAGCTTTAAC"
        #'parent': ['AY548484']
        }
        self.assertEqual(obs['n'], exp_node)
        self.assertEqual(obs['e'],[('NC_005946.1', 'AY548484')])



    def test_extractChildren(self):
        r=self.rec  
        rec_type=self.rec_type
        seq_type=self.seq_type
        db=self.db
        parent=parser.createTopLevelNode(r, rec_type, seq_type, db)
        obs=parser.extractChildren(r, parent['n'], seq_type, db)
        
        #NAmes
        obs_ids=[p['id'] for p in obs['n']]
        exp_ids=['YP_031579.1', 'YP_031580.1', 'YP_031581.1']
        self.assertEqual(obs_ids, exp_ids)

        obs_names=[p['name'] for p in obs['n']]
        exp_names=['putative replicating factor', 'putative myristylated membrane protein', 'orf229L protein-like protein']
        self.assertEqual(obs_names, exp_names)

        obs_go=[p['go'] for p in obs['n']]
        exp_go=[[], [], []]
        self.assertEqual(obs_go, exp_go)

        obs_ec=[p['ec'] for p in obs['n']]
        exp_ec=['', '', '']
        self.assertEqual(obs_ec, exp_ec)

        #obs_ec=[p['parent'] for p in obs['n']]
        #exp_ec=['NC_005946.1', 'NC_005946.1','NC_005946.1']
        #self.assertEqual(obs_ec, exp_ec)

        obs_edges=obs['e']
        exp_edges=[('YP_031579.1', 'NC_005946.1'), 
        ('YP_031580.1', 'NC_005946.1'), 
        ('YP_031581.1', 'NC_005946.1')]
        self.assertEqual(obs_edges, exp_edges)



    def test_rs_nuc_parseRecord(self):
        r=self.rec  
        seq_type=self.seq_type
        db=self.db
        obs_node, obs_edge=parser.parseRecord(r, db, seq_type)
        exp_edge=[('NC_005946.1', 'AY548484'), 
        ('YP_031579.1', 'NC_005946.1'), 
        ('YP_031580.1', 'NC_005946.1'), 
        ('YP_031581.1', 'NC_005946.1')]
        exp_nodes=4
        self.assertEqual(len(obs_node), exp_nodes)
        self.assertEqual(obs_edge, exp_edge)





