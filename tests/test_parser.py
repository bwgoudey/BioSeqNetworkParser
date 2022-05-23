from sqlite3 import dbapi2
import unittest   # The test framework
from importlib import reload  # Python 3.4+
from netdbqual import parser
from Bio import SeqIO
from io import StringIO

class Test_Parser(unittest.TestCase):
    
    def setUp(self):
        self.record_str="""
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
        self.rec_prot_top=SeqIO.parse(StringIO(self.record_str), "gb").__iter__().__next__()
        self.rec_prot_top_db='genbank'
        self.rec_prot_top_seq_type="protein"
        self.rec_prot_rec_type="top"

        self.rec_prot_top_rs_str="""
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
        self.rec_prot_top_rs=SeqIO.parse(StringIO(self.rec_prot_top_rs_str), "gb").__iter__().__next__()
        self.rec_prot_top_rs_db='refseq'
        self.rec_prot_top_rs_seq_type="protein"
        self.rec_prot_top_rs_rec_type="top"

        self.rec_nuc_top_rs_str="""
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
        self.rec_nuc_top_rs=SeqIO.parse(StringIO(self.rec_nuc_top_rs_str), "gb").__iter__().__next__()
        self.rec_nuc_top_rs_db='refseq'
        self.rec_nuc_top_rs_seq_type="nucleotide"
        self.rec_nuc_rs_rec_type="top"

        self.rec_uniprot_str="""
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
        self.rec_uniprot=SeqIO.parse(StringIO(self.rec_uniprot_str), "gb").__iter__().__next__()
        self.rec_uniprot_db='uniprot'
        self.rec_uniprot_seq_type="protein"
        self.rec_uniprot_rec_type="top"
        self.rec_uniprot_dbsource_str=parser.processUniProtDBsource(self.rec_uniprot.annotations['db_source'])


      





    # What sort of record? Its a protein, described at top-level 
    # (i.e whole record is about this protein))
    def test_rs_prot_top_determineRecordType(self):
        obs_rec_type=parser.determineRecordType(self.rec_prot_top)
        exp_rec_type="top"
        self.assertEqual(obs_rec_type, exp_rec_type)        


    # Is it a real protein? Discard Psuedo proteins
    def test_gb_prot_is_pseudo(self):
        obs_pseudo=parser.isPseudo(self.rec_prot_top, self.rec_prot_rec_type, self.rec_prot_top_seq_type)
        exp_pseudo=False
        self.assertEqual(obs_pseudo, exp_pseudo)  


    # What organism are we looking at
    def test_gb_prot_extractOrganism(self):
        obs_pseudo=parser.extractOrganism(self.rec_prot_top.annotations)
        exp_pseudo="Staphylococcus epidermidis"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_gb_prot_extractTaxonomy(self):
        obs_pseudo=parser.extractTaxonomy(self.rec_prot_top.annotations)
        exp_pseudo='Bacteria,Firmicutes,Bacilli,Bacillales,Staphylococcaceae,Staphylococcus'
        
        obs_pseudo_short=parser.extractTaxonomy({'taxonomy':['Bacteria', 'Firmicutes', 'Bacilli']})
        exp_pseudo_short='Bacteria,Firmicutes,Bacilli,,,'
        
        self.assertEqual(obs_pseudo, exp_pseudo)  
        self.assertEqual(obs_pseudo_short, exp_pseudo_short)  


    def test_gpkDateStrToNum(self):
        obs_date=parser.gpkDateStrToNum("31-MAR-2020")
        exp_date=20200331
        self.assertEqual(obs_date, exp_date)

        obs_date=parser.gpkDateStrToNum("30-NOV-2020")
        exp_date=20201130
        self.assertEqual(obs_date, exp_date)


    #Figure out when thi was first submitted
    def test_gb_prot_extractDateModified(self):
        obs=parser.extractDateModified(self.rec_prot_top, self.rec_prot_rec_type)
        exp_upload=20220331
        exp_last_mod=20220331
        exp_nmodify=1
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)   


    def test_gb_prot_extractDateLastModified(self):
        obs_date=parser.extractDateLastModified(self.rec_prot_top.annotations)
        exp_date=20220407
        self.assertEqual(obs_date, exp_date)  

    def test_gb_prot_extractGo(self):
       
        r=self.rec_prot_top
        prot1={f.type:f for f in (r.features[1],r.features[2])}
        obs_go=parser.extractGO(self.rec_prot_top, prot1, self.rec_prot_top_db)
        exp_go=sorted(['GO:0006508', 'GO:0045148'])
        self.assertEqual(obs_go, exp_go)  

    def test_gb_prot_extractEC(self):
        r=self.rec_prot_top
        prot1={f.type:f for f in (r.features[1],r.features[2])}
        obs_ec=parser.extractEC(prot1['Protein'])
        exp_ec="3.4.11.4"
        self.assertEqual(obs_ec, exp_ec)  

    def test_gb_prot_extractTaxaID(self):
        obs_taxid=parser.extractTaxaID(self.rec_prot_top.features[0])
        exp_taxid="1282"
        self.assertEqual(obs_taxid, exp_taxid)

    def test_gb_prot_extractParentEdges(self):
        obs_edges=parser.extractParentEdges(self.rec_prot_top, self.rec_prot_top_db)
        exp_edges=['CP094865.1']
        self.assertEqual(obs_edges, exp_edges)

    # def test_gb_prot_extractEdges(self):
    #     obs_edges=parser.extractEdges(self.rec_prot_top, self.rec_prot_top_db)
    #     exp_edges=[(self.rec_prot_top.id, 'CP094865.1')]
    #     self.assertEqual(obs_edges, exp_edges)

    def test_gb_prot_top_createTopLevelNode(self):
        r=self.rec_prot_top  
        rec_type=self.rec_prot_rec_type
        seq_type=self.rec_prot_top_seq_type
        db=self.rec_prot_top_db
        obs=parser.createTopLevelNode(r, rec_type, seq_type, db)
        obs['n']['seq']=obs['n']['seq'][0:10]
        
        exp_node={'id': 'UOH70005.1', 
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
        'seq': 'MKKQIIERLT'
        #'parent': ['CP094865.1']
        }

        self.assertEqual(obs['n'], exp_node)
        self.assertEqual(obs['e'], [('UOH70005.1', 'CP094865.1')])


        ##########
         # What sort of record? Its a protein, described at top-level 
    # (i.e whole record is about this protein))
    def test_rs_prot_top_determineRecordType(self):
        obs_rec_type=parser.determineRecordType(self.rec_prot_top_rs)
        exp_rec_type="top"
        self.assertEqual(obs_rec_type, exp_rec_type)        


    # Is it a real protein? Discard Psuedo proteins
    def test_rs_prot_top_is_pseudo(self):
        obs_pseudo=parser.isPseudo(self.rec_prot_top_rs, self.rec_prot_rec_type, self.rec_prot_top_rs_seq_type)
        exp_pseudo=False
        self.assertEqual(obs_pseudo, exp_pseudo)  


    # What organism are we looking at
    def test_rs_prot_top_extractOrganism(self):
        obs_pseudo=parser.extractOrganism(self.rec_prot_top_rs.annotations)
        exp_pseudo="Escherichia coli str. K-12 substr. MG1655"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_rs_prot_top_extractTaxonomy(self):
        obs_pseudo=parser.extractTaxonomy(self.rec_prot_top_rs.annotations)
        exp_pseudo='Bacteria,Proteobacteria,Gammaproteobacteria,Enterobacterales,Enterobacteriaceae,Escherichia' 
        self.assertEqual(obs_pseudo, exp_pseudo)  

    #Figure out when thi was first submitted
    def test_rs_prot_top_extractDateModified(self):
        obs=parser.extractDateModified(self.rec_prot_top_rs, self.rec_prot_rec_type)
        exp_upload=19970116
        exp_last_mod=20220308
        exp_nmodify=11
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)   


    def test_rs_prot_top_extractDateLastModified(self):
        obs_date=parser.extractDateLastModified(self.rec_prot_top_rs.annotations)
        exp_date=20220309
        self.assertEqual(obs_date, exp_date)  

  
    def test_rs_prot_top_extractGo(self):
        r=self.rec_prot_top_rs
        prot1={f.type:f for f in (r.features[1],r.features[4])}
        obs_go=parser.extractGO(self.rec_prot_top_rs, prot1, self.rec_prot_top_db)
        #obs_go=parser.extractGO(self.rec_prot_top_rs, self.rec_prot_rec_type)
        exp_go=[]
        self.assertEqual(obs_go, exp_go)  

    def test_rs_prot_top_extracEC(self):
        r=self.rec_prot_top_rs
        prot1={f.type:f for f in (r.features[1],r.features[4])}
        obs_ec=parser.extractEC(prot1['CDS'])
        #obs_ec=parser.extractEC(self.rec_prot_top_rs, self.rec_prot_rec_type)
        exp_ec=""
        self.assertEqual(obs_ec, exp_ec)  

    def test_rs_prot_top_extractTaxaID(self):
        obs_taxid=parser.extractTaxaID(self.rec_prot_top_rs.features[0])
        exp_taxid="511145"
        self.assertEqual(obs_taxid, exp_taxid)

    def test_rs_prot_top_extractParentEdges(self):
        r_annots=self.rec_prot_top_rs
        obs_edges=parser.extractParentEdges(r_annots,self.rec_prot_top_rs_db)
        exp_edges=['AAC74211'] 
        self.assertEqual(obs_edges, exp_edges)        


    def test_rs_prot_top_createTopLevelNode(self):
            r=self.rec_prot_top_rs  
            rec_type=self.rec_prot_top_rs_rec_type
            seq_type=self.rec_prot_top_rs_seq_type
            db=self.rec_prot_top_rs_db
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




    def test_ExtractRefSeqParentEdge(self):
        x="REVIEWED REFSEQ: This record has been curated by NCBI staff. The\nreference sequence was derived from BY793553.1, AK134931.1,\nBF018223.1 and BG807829.1."        
        
        obs_refs=parser.extractRefSeqParentEdge(x)
        exp_refs=["BY793553.1", "AK134931.1","BF018223.1", "BG807829.1"]
        self.assertEqual(obs_refs, exp_refs)

        x="The reference sequence was derived from D16669.1."
        obs_refs=parser.extractRefSeqParentEdge(x)
        exp_refs=["D16669.1"]
        self.assertEqual(obs_refs, exp_refs)

        x='PROVISIONAL REFSEQ: This record has not yet been subject to final\nNCBI review. The reference sequence was derived from AY548484.\nOn Jul 12, 2012 this sequence version replaced NC_012637.1.\nCOMPLETENESS: full length.'
        obs_refs=parser.extractRefSeqParentEdge(x)
        exp_refs=["AY548484"]
        self.assertEqual(obs_refs, exp_refs)

    ##########
         # What sort of record? Its a protein, described at top-level 
    # (i.e whole record is about this protein))
    def test_rs_nuc_top_determineRecordType(self):
        obs_rec_type=parser.determineRecordType(self.rec_nuc_top_rs)
        exp_rec_type="top"
        self.assertEqual(obs_rec_type, exp_rec_type)        


    # Is it a real protein? Discard Psuedo proteins
    def test_rs_nuc_top_is_pseudo(self):
        obs_pseudo=parser.isPseudo(self.rec_nuc_top_rs, self.rec_nuc_rs_rec_type, self.rec_nuc_top_rs_seq_type)
        exp_pseudo=False
        self.assertEqual(obs_pseudo, exp_pseudo)  


    # What organism are we looking at
    def test_rs_nuc_top_extractOrganism(self):
        obs_pseudo=parser.extractOrganism(self.rec_nuc_top_rs.annotations)
        exp_pseudo="Frog virus 3"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_rs_nuc_top_extractTaxonomy(self):
        obs_pseudo=parser.extractTaxonomy(self.rec_nuc_top_rs.annotations)
        exp_pseudo='Viruses,Varidnaviria,Bamfordvirae,Nucleocytoviricota,Megaviricetes,Pimascovirales' 
        self.assertEqual(obs_pseudo, exp_pseudo)  

    #Figure out when thi was first submitted
    def test_rs_nuc_top_extractDateModified(self):
        obs=parser.extractDateModified(self.rec_nuc_top_rs, self.rec_nuc_rs_rec_type)
        exp_upload=20040218#18-FEB-2004
        exp_last_mod=20040818#18-AUG-2004
        exp_nmodify=2
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)  


    def test_rs_nuc_top_extractDateLastModified(self):
        obs_date=parser.extractDateLastModified(self.rec_nuc_top_rs.annotations)
        exp_date=20201220
        self.assertEqual(obs_date, exp_date)  
    
    def test_rs_nuc_top_extractGo(self):
        r=self.rec_nuc_top_rs
        prot1={f.type:f for f in (r.features[1],r.features[2])}
        obs_go=parser.extractGO(r, prot1, self.rec_prot_top_db)

        exp_go=[]
        self.assertEqual(obs_go, exp_go)  

    def test_rs_nuc_top_extractEC(self):
        obs_ec=parser.extractEC(self.rec_nuc_top_rs.features[1])
        exp_ec=""
        self.assertEqual(obs_ec, exp_ec)  

    def test_rs_nuc_top_extractTaxaID(self):
        obs_taxid=parser.extractTaxaID(self.rec_nuc_top_rs.features[0])
        exp_taxid="10493"
        self.assertEqual(obs_taxid, exp_taxid)

    def test_rs_nuc_top_extractParentEdges(self):
        r=self.rec_nuc_top_rs        
        obs_edges=parser.extractParentEdges(r, self.rec_nuc_top_rs_db)
        exp_edges=['AY548484'] 
        self.assertEqual(obs_edges, exp_edges)   

    # def test_rs_nuc_top_extractEdges(self):
    #     r=self.rec_nuc_top_rs
    #     fd={f.type:f for f in (r.features[1],r.features[2])}
    #     obs_edges=parser.extractRelatedEdges(fd)
    #     exp_edges={"YP_031579.1":[]}
    #     self.assertEqual(obs_edges, exp_edges)    

    def test_rs_nuc_top_identifyProteins(self):
        obs_prots=parser.identifyProteins(self.rec_nuc_top_rs)
        obs_nprots=len(obs_prots)
        exp_nprots=3
        self.assertEqual(obs_nprots, exp_nprots)

    def test_rs_nuc_top_createTopLevelNode(self):
        r=self.rec_nuc_top_rs  
        rec_type=self.rec_nuc_rs_rec_type
        seq_type=self.rec_nuc_top_rs_seq_type
        db=self.rec_nuc_top_rs_db
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



    def test_rs_nuc_top_extractChildren(self):
        r=self.rec_nuc_top_rs  
        rec_type=self.rec_nuc_rs_rec_type
        seq_type=self.rec_nuc_top_rs_seq_type
        db=self.rec_nuc_top_rs_db
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
        r=self.rec_nuc_top_rs  
        seq_type=self.rec_nuc_top_rs_seq_type
        db=self.rec_nuc_top_rs_db
        obs_node, obs_edge=parser.parseRecord(r, db, seq_type)
        exp_edge=[('NC_005946.1', 'AY548484'), 
        ('YP_031579.1', 'NC_005946.1'), 
        ('YP_031580.1', 'NC_005946.1'), 
        ('YP_031581.1', 'NC_005946.1')]
        exp_nodes=4
        self.assertEqual(len(obs_node), exp_nodes)
        self.assertEqual(obs_edge, exp_edge)




        ##########
         # What sort of record? Its a protein, described at top-level 
    # (i.e whole record is about this protein))
    def test_uniprot_top_determineRecordType(self):
        obs_rec_type=parser.determineRecordType(self.rec_uniprot)
        exp_rec_type="top"
        self.assertEqual(obs_rec_type, exp_rec_type)        


    # Is it a real protein? Discard Psuedo proteins
    def test_uniprot_top_is_pseudo(self):
        obs_pseudo=parser.isPseudo(self.rec_uniprot, self.rec_prot_rec_type, self.rec_uniprot_seq_type)
        exp_pseudo=False
        self.assertEqual(obs_pseudo, exp_pseudo)  


    # What organism are we looking at
    def test_uniprot_top_extractOrganism(self):
        obs_pseudo=parser.extractOrganism(self.rec_uniprot.annotations)
        exp_pseudo="Listeria welshimeri serovar 6b str. SLCC5334"
        self.assertEqual(obs_pseudo, exp_pseudo)  

    # What are the 4 highest levels of taxa. Might help with plotting. 
    def test_uniprot_top_extractTaxonomy(self):
        obs_pseudo=parser.extractTaxonomy(self.rec_uniprot.annotations)
        exp_pseudo='Bacteria,Firmicutes,Bacilli,Bacillales,Listeriaceae,Listeria' 
        self.assertEqual(obs_pseudo, exp_pseudo)  

    #Figure out when thi was first submitted
    def test_uniprot_top_extractDateModified(self):
        db_str=parser.processUniProtDBsource(self.rec_uniprot.annotations['db_source'])
        obs=parser.extractDateModified(self.rec_uniprot, self.rec_uniprot_rec_type, db_str)
        exp_upload=20080115
        exp_last_mod=20220223
        exp_nmodify=-1
        self.assertEqual(obs["first_upload"], exp_upload)  
        self.assertEqual(obs["last_modify"], exp_last_mod)  
        self.assertEqual(obs["num_modify"], exp_nmodify)   


    def test_uniprot_top_extractDateLastModified(self):
        obs_date=parser.extractDateLastModified(self.rec_uniprot.annotations)
        exp_date=20220223
        self.assertEqual(obs_date, exp_date)  

  
    def test_uniprot_top_extractGo(self):
        r=self.rec_uniprot
        prot1={f.type:f for f in (r.features[1],r.features[4])}
        obs_go=parser.extractGO(self.rec_uniprot, prot1, self.rec_uniprot_db)
        #obs_go=parser.extractGO(self.rec_uniprot, self.rec_prot_rec_type)
        exp_go=['GO:0005737', 'GO:0008237', 'GO:0008270', 'GO:0043171', 'GO:0045148']
        self.assertEqual(obs_go, exp_go)  

    def test_uniprot_top_extracEC(self):
        r=self.rec_uniprot
        prot1={f.type:f for f in (r.features[2],r.features[4])}
        obs_ec=parser.extractEC(prot1['Protein'])
        #obs_ec=parser.extractEC(self.rec_uniprot, self.rec_prot_rec_type)
        exp_ec="3.4.11.4"
        self.assertEqual(obs_ec, exp_ec)  

    def test_uniprot_top_extractTaxaID(self):
        obs_taxid=parser.extractTaxaID(self.rec_uniprot.features[0])
        exp_taxid="386043"
        self.assertEqual(obs_taxid, exp_taxid)

    def test_uniprot_top_extractParentEdges(self):
        db_str=parser.processUniProtDBsource(self.rec_uniprot.annotations['db_source'])
        r=self.rec_uniprot
        obs_edges=parser.extractParentEdges(r,self.rec_uniprot_db, db_str)
        exp_edges=['AM263198.1', 'CAK21216.1', 'WP_011702575.1']
        self.assertEqual(obs_edges, exp_edges)        


    def test_uniprot_top_createTopLevelNode(self):
            r=self.rec_uniprot  
            rec_type=self.rec_uniprot_rec_type
            seq_type=self.rec_uniprot_seq_type
            db=self.rec_uniprot_db
            dbsrc_str=self.rec_uniprot_dbsource_str
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
        r=self.rec_uniprot  
        db=self.rec_uniprot_db
        seq_type=self.rec_uniprot_seq_type
        obs_node, obs_edge=parser.parseRecord(r, db, seq_type)
        exp_edge=[('A0AJN4.1', 'AM263198.1'), ('A0AJN4.1', 'CAK21216.1'), ('A0AJN4.1', 'WP_011702575.1')]
        exp_nodes=1
        self.assertEqual(len(obs_node), exp_nodes)
        self.assertEqual(obs_edge, exp_edge)



