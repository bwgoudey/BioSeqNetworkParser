import unittest   # The test framework
from netdbqual import rec_parser
from importlib import reload  # Python 3.4+

class Test_ParseInference(unittest.TestCase):

    def test_homology_with_qualifier(self):
        inf="non-experimental evidence, no additional details recorded"
        note="similar to GB:AAL97515.1 (AE010016) percent identity 98 in 407 aa"

        obs=rec_parser.parseInference(inf, note)
        exp={'acc':'AAL97515',
            'seq_version':1,
            'identity': .98,
            'db':'g', 
            'type':'h',
            'model':""
            }
        self.assertEqual(obs, exp)

    def test_no_result(self):
        inf="non-experimental evidence, no additional details recorded"
        note="catalyzes the release of the N-terminal amino acid from a tripeptide"

        obs=rec_parser.parseInference(inf, note)
        exp={}
        self.assertEqual(obs, exp)

    #Here, YPO1631 refers to a  locus_tag within a protein https://www.ncbi.nlm.nih.gov/protein/CAL20276
    def test_ImpossibleLocus(self):
        inf="non-experimental evidence, no additional details recorded"
        note="""Similar to: Y. pestis YPO1631 pepT; peptidase T
                     (100% evalue=0); E. coli Z1832 pepT; putative peptidase T
                     (77.1% evalue=0)"""

        obs=rec_parser.parseInference(inf, note)
        exp=exp={'acc':'YPO1631',
            'seq_version':"",
            'identity': 1,
            'db':'g', 
            'type':'h',
            'model':""}
        self.assertEqual(obs, exp)


    def test_multiple_inference(self):
        inf="""non-experimental evidence, no additional details recorded"
                     (100% evalue=0); E. coli Z1832 pepT; putative peptidase T
                     (77.1% evalue=0)"""
        inf2="protein motif:HMMPfam:IPR011650"
        inf3="protein motif:HMMTigr:IPR010161"
        inf4="protein motif:ScanRegExp:IPR001261"
        inf5="similar to AA sequence:INSD:AAL20156.1"
        note="""KEGG: stm:STM1227 4.3e-211 pepT; peptidase T
                     K01258;
                     COG: COG2195 Di- and tripeptidases;
                     Psort location: Cytoplasmic, score:9.97"""

        obs=rec_parser.parseInference(inf, note)
        exp=exp={'acc':'YPO1631',
            'seq_version':"",
            'identity': 1,
            'db':'g', 
            'type':'h',
            'model':""}
        self.assertEqual(obs, exp)


    def test_ab_initio(self):
            inf="ab initio prediction:AMIGene:2.0"
            obs=rec_parser.parseInference(inf)
            exp=exp={'acc':'',
                'seq_version':"",
                'identity': "",
                'db':'', 
                'type':'d',
                'model':"AMIGene:2.0"}
            self.assertEqual(obs, exp)        

    

    def test_ab_initio(self):
        inf="protein motif:TFAM:TIGR01882"
        note="""KEGG: gtn:GTNG_1659 peptidase T;
                        TIGRFAM: peptidase T;
                        PFAM: peptidase dimerisation domain protein; peptidase
                        M20"""
        obs=rec_parser.parseInference(inf)
        exp=exp={'acc':':TIGR01882',
                'seq_version':"",
                'identity': "",
                'model':"tfam",
                'db':'tfam', 
                'type':'hmm'}
        self.assertEqual(obs, exp)

    def test_refseq_homology(self):
        inf="""COORDINATES: similar to AA
                        sequence:RefSeq:WP_014656343.2"""
        note="""Derived by automated computational analysis using
                gene prediction method: Protein Homology."""
        obs=rec_parser.parseInference(inf, note)
        exp={'acc':'WP_014656343',
                'seq_version':2,
                'identity': "",
                'db':'r', 
                'type':'h',
                'model':""}
        self.assertEqual(obs, exp)

    def test_refseq_abinit(self):
        inf="""COORDINATES: ab initio prediction:GeneMarkS+"""
        note="""Derived by automated computational analysis using
                     gene prediction method: GeneMarkS+."""
        obs=rec_parser.parseInference(inf, note)
        exp={'acc':'',
                'seq_version':"",
                'identity': "",
                'db':'r', 
                'type':'h',
                "model":"GeneMarkS+"}
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()