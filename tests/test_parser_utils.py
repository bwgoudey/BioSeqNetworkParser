import unittest   # The test framework
from importlib import reload  # Python 3.4+
from netdbqual import rec_parser


class Test_Parser_utils(unittest.TestCase):
    
    def test_gpkDateStrToNum(self):
        obs_date=rec_parser.gpkDateStrToNum("31-MAR-2020")
        exp_date=20200331
        self.assertEqual(obs_date, exp_date)

        obs_date=rec_parser.gpkDateStrToNum("30-NOV-2020")
        exp_date=20201130
        self.assertEqual(obs_date, exp_date)


    def test_ExtractRefSeqParentEdge(self):
        x="REVIEWED REFSEQ: This record has been curated by NCBI staff. The\nreference sequence was derived from BY793553.1, AK134931.1,\nBF018223.1 and BG807829.1."        
        
        obs_refs=rec_parser.extractRefSeqParentEdge(x)
        exp_refs=["BY793553.1", "AK134931.1","BF018223.1", "BG807829.1"]
        self.assertEqual(obs_refs, exp_refs)

        x="The reference sequence was derived from D16669.1."
        obs_refs=rec_parser.extractRefSeqParentEdge(x)
        exp_refs=["D16669.1"]
        self.assertEqual(obs_refs, exp_refs)

        x='PROVISIONAL REFSEQ: This record has not yet been subject to final\nNCBI review. The reference sequence was derived from AY548484.\nOn Jul 12, 2012 this sequence version replaced NC_012637.1.\nCOMPLETENESS: full length.'
        obs_refs=rec_parser.extractRefSeqParentEdge(x)
        exp_refs=["AY548484"]
        self.assertEqual(obs_refs, exp_refs)
