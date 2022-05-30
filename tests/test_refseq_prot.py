from sqlite3 import dbapi2
import unittest   # The test framework
from importlib import reload  # Python 3.4+
from netdbqual import parser
from Bio import SeqIO
from io import StringIO

class Test_Parser(unittest.TestCase):
    
    def setUp(self):
        self.rec_str="""
LOCUS       NP_415645                408 aa            linear   CON 09-MAR-2022
DEFINITION  peptidase T [Escherichia coli str. K-12 substr. MG1655].
ACCESSION   NP_415645
VERSION     NP_415645.1
DBLINK      BioProject: PRJNA57779
            BioSample: SAMN02604091
DBSOURCE    REFSEQ: accession NC_000913.3
KEYWORDS    RefSeq.
SOURCE      Escherichia coli str. K-12 substr. MG1655
  ORGANISM  Escherichia coli str. K-12 substr. MG1655
            Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales;
            Enterobacteriaceae; Escherichia.
REFERENCE   1  (residues 1 to 408)
  AUTHORS   Riley,M., Abe,T., Arnaud,M.B., Berlyn,M.K., Blattner,F.R.,
            Chaudhuri,R.R., Glasner,J.D., Horiuchi,T., Keseler,I.M., Kosuge,T.,
            Mori,H., Perna,N.T., Plunkett,G. III, Rudd,K.E., Serres,M.H.,
            Thomas,G.H., Thomson,N.R., Wishart,D. and Wanner,B.L.
  TITLE     Escherichia coli K-12: a cooperatively developed annotation
            snapshot--2005
  JOURNAL   Nucleic Acids Res. 34 (1), 1-9 (2006)
   PUBMED   16397293
  REMARK    Publication Status: Online-Only
REFERENCE   2  (residues 1 to 408)
  AUTHORS   Hayashi,K., Morooka,N., Yamamoto,Y., Fujita,K., Isono,K., Choi,S.,
            Ohtsubo,E., Baba,T., Wanner,B.L., Mori,H. and Horiuchi,T.
  TITLE     Highly accurate genome sequences of Escherichia coli K-12 strains
            MG1655 and W3110
  JOURNAL   Mol. Syst. Biol. 2, 2006 (2006)
   PUBMED   16738553
REFERENCE   3  (residues 1 to 408)
  AUTHORS   Blattner,F.R., Plunkett,G. III, Bloch,C.A., Perna,N.T., Burland,V.,
            Riley,M., Collado-Vides,J., Glasner,J.D., Rode,C.K., Mayhew,G.F.,
            Gregor,J., Davis,N.W., Kirkpatrick,H.A., Goeden,M.A., Rose,D.J.,
            Mau,B. and Shao,Y.
  TITLE     The complete genome sequence of Escherichia coli K-12
  JOURNAL   Science 277 (5331), 1453-1462 (1997)
   PUBMED   9278503
REFERENCE   4  (residues 1 to 408)
  AUTHORS   Arnaud,M., Berlyn,M.K.B., Blattner,F.R., Galperin,M.Y.,
            Glasner,J.D., Horiuchi,T., Kosuge,T., Mori,H., Perna,N.T.,
            Plunkett,G. III, Riley,M., Rudd,K.E., Serres,M.H., Thomas,G.H. and
            Wanner,B.L.
  TITLE     Workshop on Annotation of Escherichia coli K-12
  JOURNAL   Unpublished
  REMARK    Woods Hole, Mass., on 14-18 November 2003 (sequence corrections)
REFERENCE   5  (residues 1 to 408)
  AUTHORS   Glasner,J.D., Perna,N.T., Plunkett,G. III, Anderson,B.D.,
            Bockhorst,J., Hu,J.C., Riley,M., Rudd,K.E. and Serres,M.H.
  TITLE     ASAP: Escherichia coli K-12 strain MG1655 version m56
  JOURNAL   Unpublished
  REMARK    ASAP download 10 June 2004 (annotation updates)
REFERENCE   6  (residues 1 to 408)
  AUTHORS   Hayashi,K., Morooka,N., Mori,H. and Horiuchi,T.
  TITLE     A more accurate sequence comparison between genomes of Escherichia
            coli K12 W3110 and MG1655 strains
  JOURNAL   Unpublished
  REMARK    GenBank accessions AG613214 to AG613378 (sequence corrections)
REFERENCE   7  (residues 1 to 408)
  AUTHORS   Perna,N.T.
  TITLE     Escherichia coli K-12 MG1655 yqiK-rfaE intergenic region, genomic
            sequence correction
  JOURNAL   Unpublished
  REMARK    GenBank accession AY605712 (sequence corrections)
REFERENCE   8  (residues 1 to 408)
  AUTHORS   Rudd,K.E.
  TITLE     A manual approach to accurate translation start site annotation: an
            E. coli K-12 case study
  JOURNAL   Unpublished
REFERENCE   9  (residues 1 to 408)
  CONSRTM   NCBI Genome Project
  TITLE     Direct Submission
  JOURNAL   Submitted (08-MAR-2022) National Center for Biotechnology
            Information, NIH, Bethesda, MD 20894, USA
REFERENCE   10 (residues 1 to 408)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (30-JUL-2014) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Protein update by submitter
REFERENCE   11 (residues 1 to 408)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (15-NOV-2013) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Protein update by submitter
REFERENCE   12 (residues 1 to 408)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (26-SEP-2013) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Sequence update by submitter
REFERENCE   13 (residues 1 to 408)
  AUTHORS   Rudd,K.E.
  TITLE     Direct Submission
  JOURNAL   Submitted (06-FEB-2013) Department of Biochemistry and Molecular
            Biology, University of Miami Miller School of Medicine, 118 Gautier
            Bldg., Miami, FL 33136, USA
  REMARK    Sequence update by submitter
REFERENCE   14 (residues 1 to 408)
  AUTHORS   Rudd,K.E.
  TITLE     Direct Submission
  JOURNAL   Submitted (24-APR-2007) Department of Biochemistry and Molecular
            Biology, University of Miami Miller School of Medicine, 118 Gautier
            Bldg., Miami, FL 33136, USA
  REMARK    Annotation update from ecogene.org as a multi-database
            collaboration
REFERENCE   15 (residues 1 to 408)
  AUTHORS   Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (07-FEB-2006) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Protein updates by submitter
REFERENCE   16 (residues 1 to 408)
  AUTHORS   Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (10-JUN-2004) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
  REMARK    Sequence update by submitter
REFERENCE   17 (residues 1 to 408)
  AUTHORS   Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (13-OCT-1998) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
REFERENCE   18 (residues 1 to 408)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (02-SEP-1997) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
REFERENCE   19 (residues 1 to 408)
  AUTHORS   Blattner,F.R. and Plunkett,G. III.
  TITLE     Direct Submission
  JOURNAL   Submitted (16-JAN-1997) Laboratory of Genetics, University of
            Wisconsin, 425G Henry Mall, Madison, WI 53706-1580, USA
COMMENT     PROVISIONAL REFSEQ: This record has not yet been subject to final
            NCBI review. The reference sequence is identical to AAC74211.
            Changes to proteins and annotation made on March 7, 2022.  Current
            U00096 annotation updates are derived from EcoCyc
            https://ecocyc.org/.  Suggestions for updates can be sent to
            biocyc-support@ai.sri.com. These updates are being generated from a
            collaboration  that includes EcoCyc, the University of Wisconsin,
            UniProtKB/Swiss-Prot, and the National Center for Biotechnology
            Information (NCBI).
            COMPLETENESS: full length.
            Method: conceptual translation.
FEATURES             Location/Qualifiers
     source          1..408
                     /organism="Escherichia coli str. K-12 substr. MG1655"
                     /strain="K-12"
                     /sub_strain="MG1655"
                     /db_xref="taxon:511145"
     Protein         1..408
                     /product="peptidase T"
                     /EC_number="3.4.11.4"
                     /calculated_mol_wt=44792
     Region          1..408
                     /region_name="PRK05469"
                     /note="tripeptide aminopeptidase PepT"
                     /db_xref="CDD:235484"
     Site            order(78,140,173..174,196,379)
                     /site_type="metal-binding"
                     /note="metal binding site [ion binding]"
                     /db_xref="CDD:349887"
     CDS             1..408
                     /gene="pepT"
                     /locus_tag="b1127"
                     /gene_synonym="ECK1113"
                     /coded_by="NC_000913.3:1185844..1187070"
                     /transl_table=11
                     /db_xref="UniProtKB/Swiss-Prot:P29745"
                     /db_xref="ASAP:ABE-0003800"
                     /db_xref="ECOCYC:EG11549"
                     /db_xref="GeneID:946333"
CONTIG      join(WP_000359434.1:1..408)
ORIGIN
        1 mdkllerfln yvsldtqska gvrqvpsteg qwkllhllke qleemglinv tlsekgtlma
       61 tlpanvpgdi paigfishvd tspdcsgknv npqivenyrg gdialgigde vlspvmfpvl
      121 hqllgqtlit tdgktllgad dkagiaeimt alavlqqkki phgdirvaft pdeevgkgak
      181 hfdvdafdar waytvdgggv gelefenfna asvnikivgn nvhpgtakgv mvnalslaar
      241 ihaevpades pemtegyegf yhlasmkgtv eradmhyiir dfdrkqfear krkmmeiakk
      301 vgkglhpdcy ielviedsyy nmrekvvehp hildiaqqam rdcdiepelk pirggtdgaq
      361 lsfmglpcpn lftggynyhg khefvtlegm ekavqvivri aeltaqrk
//
"""
        self.rec=SeqIO.parse(StringIO(self.rec_str), "gb").__iter__().__next__()
        self.db='refseq'
        self.seq_type="protein"
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
        exp_pseudo="Escherichia coli str. K-12 substr. MG1655"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_extractTaxonomy(self):
        obs_pseudo=parser.extractTaxonomy(self.rec.annotations)
        exp_pseudo='Bacteria,Proteobacteria,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia' 
        self.assertEqual(obs_pseudo, exp_pseudo)  

    #Figure out when thi was first submitted
    def test_extractDateModified(self):
        obs=parser.extractDateModified(self.rec, "")
        exp_upload=19970116
        exp_last_mod=20220308
        exp_nmodify=11
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)   


    def test_extractDateLastModified(self):
        obs_date=parser.extractDateLastModified(self.rec.annotations)
        exp_date=20220309
        self.assertEqual(obs_date, exp_date)  

  
    def test_extractGo(self):
        prot1={f.type:f for f in (self.rec.features[1],self.rec.features[4])}
        obs_go=parser.extractGO(self.rec, prot1, self.db)
        #obs_go=parser.extractGO(self.rec, self.rec_prot_rec_type)
        exp_go=[]
        self.assertEqual(obs_go, exp_go)  

    def test_extracEC(self):
        r=self.rec
        prot1={f.type:f for f in (r.features[1],r.features[4])}
        obs_ec=parser.extractEC(r, prot1['Protein'])
        #obs_ec=parser.extractEC(self.rec, self.rec_prot_rec_type)
        exp_ec="3.4.11.4"
        self.assertEqual(obs_ec, exp_ec)  

    def test_extractTaxaID(self):
        obs_taxid=parser.extractTaxaID(self.rec.features[0])
        exp_taxid="511145"
        self.assertEqual(obs_taxid, exp_taxid)

    def test_extractParentEdges(self):
        r_annots=self.rec
        obs_edges=parser.extractParentEdges(r_annots,self.db)
        exp_edges=['AAC74211'] 
        self.assertEqual(obs_edges, exp_edges)        


    def test_createTopLevelNode(self):
            r=self.rec  
            rec_type=self.rec_type
            seq_type=self.seq_type
            db=self.db
            obs=parser.createTopLevelNode(r, rec_type, seq_type, db)
            obs['n']['seq']=obs['n']['seq'][0:10]
            
            exp_node={'id': 'NP_415645.1', 
            'seq_version': 1,
             'date_first_upload': 19970116, 
             'num_modified': 11, 
             'date_last_modified': 20220309,
             'organism': 'Escherichia coli str. K-12 substr. MG1655', 
             'taxa_id': '511145',
             'taxonomy': 'Bacteria,Proteobacteria,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia', 
             'n_products': 1, 
             'name': 'peptidase T', 
             'go': [],
             'ec': '3.4.11.4', 
             'seq': 'MDKLLERFLN'
             #'parent': ['AAC74211']
             }

            self.assertEqual(obs['n'], exp_node)

            exp_edges=[('NP_415645.1', 'AAC74211')]
            self.assertEqual(obs['e'], exp_edges)

