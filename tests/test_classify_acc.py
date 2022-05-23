import unittest   # The test framework
import netdbqual.classify_acc as ca
from importlib import reload  # Python 3.4+

class Test_ClassifyAcc(unittest.TestCase):
    reload(ca)

    #  genbank     contig     AABR07056638
    #  genbank     nucleotide CP093132    
    #  genbank     protein    CAK21216    
    #  pdb         protein    1FNO_       
    #  pir         protein    M62725      
    #  refseq      nucleotide NC_000913   
    #  refseq      protein    NP_460197   
    #  refseq.anon protein    WP_011702575

    def test_gb_nuc(self):
        acc="AABR07056638"
        exp_db="genbank"
        exp_type="contig"
        obs_db, obs_type=ca.classify_acc(acc)
        self.assertEqual(exp_db, obs_db)
        self.assertEqual(exp_type, obs_type)

    def test_gb_prot_8(self):
        acc="UOI52910"
        exp_db="genbank"
        exp_type="protein"
        obs_db, obs_type=ca.classify_acc(acc)
        self.assertEqual(exp_db, obs_db)
        self.assertEqual(exp_type, obs_type)


    def test_gb_nuc(self):
        acc="CP093132"
        exp_db="genbank"
        exp_type="nucleotide"
        obs_db, obs_type=ca.classify_acc(acc)
        self.assertEqual(exp_db, obs_db)
        self.assertEqual(exp_type, obs_type)

    def test_gb_prot(self):
        acc="CAK21216"
        exp_db="genbank"
        exp_type="protein"
        obs_db, obs_type=ca.classify_acc(acc)
        self.assertEqual(exp_db, obs_db)
        self.assertEqual(exp_type, obs_type)

    def test_pdb_protein(self):
        acc="1FNO_"
        exp_db="pdb"
        exp_type="protein"
        obs_db, obs_type=ca.classify_acc(acc)
        self.assertEqual(exp_db, obs_db)
        self.assertEqual(exp_type, obs_type)

    def test_pir_protein(self):
        acc="M62725"
        exp_db="pir"
        exp_type="protein"
        obs_db, obs_type=ca.classify_acc(acc)
        self.assertEqual(exp_db, obs_db)
        self.assertEqual(exp_type, obs_type)

    def test_refseq_nuc(self):
        acc="NC_000913"
        exp_db="refseq"
        exp_type="nucleotide"
        obs_db, obs_type=ca.classify_acc(acc)
        self.assertEqual(exp_db, obs_db)
        self.assertEqual(exp_type, obs_type)

    def test_refseq_prot(self):
        acc="NP_460197"
        exp_db="refseq"
        exp_type="protein"
        obs_db, obs_type=ca.classify_acc(acc)
        self.assertEqual(exp_db, obs_db)
        self.assertEqual(exp_type, obs_type)

    def test_refseq_anon_prot(self):
        acc="WP_011702575"
        exp_db="refseq.anon"
        exp_type="protein"
        obs_db, obs_type=ca.classify_acc(acc)
        self.assertEqual(exp_db, obs_db)
        self.assertEqual(exp_type, obs_type)


if __name__ == '__main__':
    unittest.main()