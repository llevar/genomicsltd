'''
Created on Jul 22, 2014

@author: siakhnin
'''
import unittest

import os, sys
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(CURRENT_DIR))

from variant_counter import VcfMutationCounter
from variant_counter import VcfRecord

class VcfMutationCounterTest(unittest.TestCase):


    def setUp(self):
        self.vcfMutationCounter = VcfMutationCounter()


    def tearDown(self):
        self.vcfMutationCounter = None


    def test_addMutationsToCount_snv(self):
        vcf_records = []
        vcf_record = VcfRecord("A", ["C"], ["ENSG00000101278"])
        
        vcf_records.append(vcf_record)
        
        self.vcfMutationCounter.addMutationsToCount(vcf_records)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["SNV"], 1)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["INS"], 0)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["DEL"], 0)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["MNV"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["SNV"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["INS"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["DEL"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["MNV"], 0)
        
    
    def test_addMutationsToCount_ins(self):
        vcf_records = []
        vcf_record = VcfRecord("ATG", ["ATGGC"], ["ENSG00000101278"])
        
        vcf_records.append(vcf_record)
        
        self.vcfMutationCounter.addMutationsToCount(vcf_records)
        
        
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["SNV"], 0)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["INS"], 1)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["DEL"], 0)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["MNV"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["SNV"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["INS"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["DEL"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["MNV"], 0)
    
    def test_addMutationsToCount_del(self):
        vcf_records = []
        vcf_record = VcfRecord("ATG", ["AT", "A"], ["ENSG00000101278"])
        
        vcf_records.append(vcf_record)
        
        self.vcfMutationCounter.addMutationsToCount(vcf_records)
        
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["SNV"], 0)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["INS"], 0)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["DEL"], 1)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["MNV"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["SNV"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["INS"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["DEL"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["MNV"], 0)
    
    def test_addMutationsToCount_mnv(self):
        vcf_records = []
        vcf_record = VcfRecord("ATG", ["CAT", "GGGAA"], ["ENSG00000101278"])
        
        vcf_records.append(vcf_record)
        
        self.vcfMutationCounter.addMutationsToCount(vcf_records)
        
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["SNV"], 0)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["INS"], 0)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["DEL"], 0)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["MNV"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["SNV"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["INS"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["DEL"], 0)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["MNV"], 1)
    
    def test_addMutationsToCount_multi(self):
        vcf_records = []
        vcf_record = VcfRecord("ATG", ["CAT", "GGGAA"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("ATG", ["ATGGA", "ATGATG"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("GGG", ["G"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("T", ["G", "C"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        
        self.vcfMutationCounter.addMutationsToCount(vcf_records)
        
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["SNV"], 1)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["INS"], 1)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["DEL"], 1)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["MNV"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["SNV"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["INS"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["DEL"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["MNV"], 1)
    
    def test_addMutationsToCount_multi_many(self):
        vcf_records = []
        vcf_record = VcfRecord("ATG", ["CAT", "GGGAA"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("ATGG", ["CGCAT", "A"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("ATG", ["ATGGA", "ATGATG"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("GGG", ["G"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("T", ["G", "C"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("T", ["A"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("T", ["G", "C"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("G", ["T"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        vcf_record = VcfRecord("C", ["G"], ["ENSG00000101278"])
        vcf_records.append(vcf_record)
        
        self.vcfMutationCounter.addMutationsToCount(vcf_records)
        
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["SNV"], 5)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["INS"], 1)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["DEL"], 1)
        self.assertEquals(self.vcfMutationCounter.global_mutation_tallies["MNV"], 2)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["SNV"], 5)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["INS"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["DEL"], 1)
        self.assertEqual(self.vcfMutationCounter.mutations_by_gene["ENSG00000101278"]["MNV"], 2)
    
    def suite(self):
        return unittest.TestLoader().loadTestsFromTestCase(VcfMutationCounterTest)   
            
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()