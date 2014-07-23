'''
Created on Jul 22, 2014

@author: siakhnin
'''
import unittest

import os, sys

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(CURRENT_DIR))

from variant_counter import VcfRecord

class VcfRecordTest(unittest.TestCase):


    def setUp(self):
        self.vcfRecord = VcfRecord(None,None,None)


    def tearDown(self):
        self.vcfRecord = None


    def test_isSnp_same_allele(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["A"]
        self.assertFalse(self.vcfRecord.isSNV())

    def test_isSnp_diff_alleles(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["G"]
        self.assertTrue(self.vcfRecord.isSNV())
        
    def test_isSnp_diff_lens_alt_bigger(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["AG"]
        self.assertFalse(self.vcfRecord.isSNV())
        
    def test_isSnp_diff_lens_ref_bigger(self):
        self.vcfRecord.ref = "AG"
        self.vcfRecord.alt = ["A"]
        self.assertFalse(self.vcfRecord.isSNV())
    
    def test_isSnp_multi_alt(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["T", "C"]
        self.assertTrue(self.vcfRecord.isSNV())
        
    def test_isSnp_multi_alt_ref_bigger(self):
        self.vcfRecord.ref = "AT"
        self.vcfRecord.alt = ["T", "C"]
        self.assertFalse(self.vcfRecord.isSNV())
        
    def test_isSnp_multi_alt_alt_bigger(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["AT", "C"]
        self.assertFalse(self.vcfRecord.isSNV())
        
    def test_isInsertion_same_allele(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["A"]
        self.assertFalse(self.vcfRecord.isInsertion())
        
    def test_isInsertion_dif_alleles(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["G"]
        self.assertFalse(self.vcfRecord.isInsertion())
        
    def test_isInsertion_single_nt(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["AG"]
        self.assertTrue(self.vcfRecord.isInsertion())
        
    def test_isInsertion_multi_nt(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["AGCCGT"]
        self.assertTrue(self.vcfRecord.isInsertion())
        
    def test_isInsertion_not_aligned(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["TAG"]
        self.assertFalse(self.vcfRecord.isInsertion())
        
    def test_isInsertion_single_nt_multi_value(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["AG", "AT"]
        self.assertTrue(self.vcfRecord.isInsertion())
        
    def test_isInsertion_multi_nt_multi_value(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["AGCCGT", "ATCG"]
        self.assertTrue(self.vcfRecord.isInsertion())
        
    def test_isInsertion_not_aligned_multi_value(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["TAG", "GACG"]
        self.assertFalse(self.vcfRecord.isInsertion())
        
    def test_isInsertion_multi_ref(self):
        self.vcfRecord.ref = "ACTG"
        self.vcfRecord.alt = ["ACTGAGC"]
        self.assertTrue(self.vcfRecord.isInsertion())
        
    def test_isInsertion_not_aligned_multi_ref(self):
        self.vcfRecord.ref = "AT"
        self.vcfRecord.alt = ["TAG"]
        self.assertFalse(self.vcfRecord.isInsertion())
        
    def test_isDeletion_same_allele(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["A"]
        self.assertFalse(self.vcfRecord.isDeletion())
        
    def test_isDeletion_dif_alleles(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["G"]
        self.assertFalse(self.vcfRecord.isDeletion())
        
    def test_isDeletion_single_nt(self):
        self.vcfRecord.ref = "AG"
        self.vcfRecord.alt = ["A"]
        self.assertTrue(self.vcfRecord.isDeletion())
        
    def test_isDeletion_multi_nt(self):
        self.vcfRecord.ref = "AGCCGT"
        self.vcfRecord.alt = ["A"]
        self.assertTrue(self.vcfRecord.isDeletion())
        
    def test_isDeletion_not_aligned(self):
        self.vcfRecord.ref = "TAG"
        self.vcfRecord.alt = ["A"]
        self.assertFalse(self.vcfRecord.isDeletion())
        
    def test_isDeletion_single_nt_multi_value(self):
        self.vcfRecord.ref = "AGT"
        self.vcfRecord.alt = ["A", "T"]
        self.assertFalse(self.vcfRecord.isDeletion())
        
    def test_isDeletion_multi_nt_multi_value(self):
        self.vcfRecord.ref = "AGTCG"
        self.vcfRecord.alt = ["AG", "AGT"]
        self.assertTrue(self.vcfRecord.isDeletion())
        
    def test_isDeletion_not_aligned_multi_value(self):
        self.vcfRecord.ref = "AT"
        self.vcfRecord.alt = ["T", "G"]
        self.assertFalse(self.vcfRecord.isDeletion())
        
    def test_getMutationType_same_allele(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["A"]
        self.assertEqual(self.vcfRecord.getMutationType(), "MNV")
        
    def test_getMutationType_dif_alleles(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["G"]
        self.assertNotEqual(self.vcfRecord.getMutationType(), "MNV")
        
    def test_getMutationType_multi_nt_same_number(self):
        self.vcfRecord.ref = "TCGAGT"
        self.vcfRecord.alt = ["AGCCGT"]
        self.assertEqual(self.vcfRecord.getMutationType(), "MNV")
    
    def test_getMutationType_multi_nt_dif_number(self):
        self.vcfRecord.ref = "TCGAGTTG"
        self.vcfRecord.alt = ["AGCCGT"]
        self.assertEqual(self.vcfRecord.getMutationType(), "MNV")
            
    def test_getMutationType_multi_nt_dif_number_single_ref(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["TAG"]
        self.assertEqual(self.vcfRecord.getMutationType(), "MNV")
        
    def test_getMutationType_multi_nt_dif_number_single_alt(self):
        self.vcfRecord.ref = "TAG"
        self.vcfRecord.alt = ["A"]
        self.assertEqual(self.vcfRecord.getMutationType(), "MNV")    
    
    def test_getMutationType_single_nt_multi_value(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["TAG", "GAAT"]
        self.assertEqual(self.vcfRecord.getMutationType(), "MNV")
        
    def test_getMutationType_snv(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["T"]
        self.assertEqual(self.vcfRecord.getMutationType(), "SNV")
        
    def test_getMutationType_ins(self):
        self.vcfRecord.ref = "A"
        self.vcfRecord.alt = ["AG", "ACT"]
        self.assertEqual(self.vcfRecord.getMutationType(), "INS")
        
    def test_getMutationType_del(self):
        self.vcfRecord.ref = "TGA"
        self.vcfRecord.alt = ["T", "TG"]
        self.assertEqual(self.vcfRecord.getMutationType(), "DEL") 
    
    def suite(self):
        return unittest.TestLoader().loadTestsFromTestCase(VcfRecordTest)   
    
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()