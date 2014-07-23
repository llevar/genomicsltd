'''
Created on Jul 22, 2014

@author: siakhnin
'''
import unittest
import os, sys

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(CURRENT_DIR))

from variant_counter import VcfParser

class VcfParserTest(unittest.TestCase):


    def setUp(self):
        self.vcf_parser = VcfParser()


    def tearDown(self):
        self.vcf_parser = None


    def test_parseVariant_snv(self):
        variant_string = "20 65900 rs6053810 G A 1649 PASS BRF=0.04;FR=1.0000;HP=2;HapScore=1;MGOF=8;MMLQ=37;MQ=61.04;NF=23;NR=22;PP=1649;QD=37.4171380566;SC=TGACACTTACGTTGCTTCTAA;SbPval=0.58;Source=Platypus;TC=45;TCF=23;TCR=22;TR=45;WE=65908;WS=65890;Gene=ENSG00000178591;SYMBOL=DEFB125;CSN=1;BIOTYPE=protein_coding;Transcript=CAN;Consequence=upstream_gene_variant;DISTANCE=2451;AF=0.7277;U10K=5599/7242;WGS500=261/274;GR=0.18;PH=0.284;PS=0.031 GT:GL:GOF:GQ:NR:NV 1/1:-168.4,-12.15,0.0:8:99:45:45"

        variant_record = self.vcf_parser.parseVariant(variant_string)
        
        self.assertIsNotNone(variant_record)
        self.assertListEqual(variant_record.gene, ["ENSG00000178591"])
        self.assertEquals(variant_record.ref, "G")
        self.assertListEqual(variant_record.alt, ["A"])
        
    def test_parseVariant_ins(self):
        variant_string = "20 67500 rs112142516 T TTGGTATCTAG 2640 PASS BRF=0.0;FR=1.0000;HP=2;HapScore=1;MGOF=0;MMLQ=37;MQ=95.71;NF=18;NR=35;PP=2640;QD=20;SC=CCTGATTTCCTTGGTATTAAA;SbPval=0.63;Source=Platypus;TC=58;TCF=21;TCR=37;TR=53;WE=67514;WS=67490;Gene=ENSG00000178591;SYMBOL=DEFB125;CSN=1;BIOTYPE=protein_coding;Transcript=CAN;Consequence=upstream_gene_variant;DISTANCE=850;U10K=5210/7242;WGS500=246/274;GR=0.199;PH=0.103;PS=0.003 GT:GL:GOF:GQ:NR:NV 1/1:-274.0,-15.3,0.0:0:99:58:53"

        variant_record = self.vcf_parser.parseVariant(variant_string)
        
        self.assertIsNotNone(variant_record)
        self.assertListEqual(variant_record.gene, ["ENSG00000178591"])
        self.assertEquals(variant_record.ref, "T")
        self.assertListEqual(variant_record.alt, ["TTGGTATCTAG"])
        
    def test_parseVariant_del(self):
        variant_string = "20 72104 rs11087028 TA T 1633 PASS BRF=0.02;FR=1.0000;HP=4;HapScore=1;MGOF=1;MMLQ=37;MQ=60.04;NF=19;NR=20;PP=1633;QD=42.7435897436;SC=TTTAAGTCTGTAAAACCATAC;SbPval=0.5;Source=Platypus;TC=46;TCF=23;TCR=23;TR=39;WE=72113;WS=72094;Gene=ENSG00000178591;SYMBOL=DEFB125;CSN=1;BIOTYPE=protein_coding;Transcript=CAN;Consequence=intron_variant;AF=0.44;U10K=3662/7242;WGS500=200/274;GR=0.468;PH=0.413;PS=0.021 GT:GL:GOF:GQ:NR:NV 1/1:-166.7,-12.64,0.0:1:99:46:39"

        variant_record = self.vcf_parser.parseVariant(variant_string)
        
        self.assertIsNotNone(variant_record)
        self.assertListEqual(variant_record.gene, ["ENSG00000178591"])
        self.assertEquals(variant_record.ref, "TA")
        self.assertListEqual(variant_record.alt, ["T"])
        
    def test_parseVariant_mnv(self):
        variant_string = "20 238530 . TATCTGGTTGC GATCTGGTTGT 1105 PASS BRF=0.0;FR=0.5000;HP=2;HapScore=2;MGOF=11;MMLQ=37;MQ=60.0;NF=10;NR=7;PP=1105;QD=68.1451691095;SC=AACACTAGCATATCTGGTTGC;SbPval=0.61;Source=Platypus;TC=43;TCF=25;TCR=18;TR=17;WE=238548;WS=238520;Gene=ENSG00000186458;SYMBOL=DEFB132;CSN=1;BIOTYPE=protein_coding;Transcript=CAN;Consequence=intron_variant;GR=1.66;PH=0.167;PS=0.000 GT:GL:GOF:GQ:NR:NV 1/0:-116.45,0.0,-173.85:11:99:43:17"

        variant_record = self.vcf_parser.parseVariant(variant_string)
        
        self.assertIsNotNone(variant_record)
        self.assertListEqual(variant_record.gene, ["ENSG00000186458"])
        self.assertEquals(variant_record.ref, "TATCTGGTTGC")
        self.assertListEqual(variant_record.alt, ["GATCTGGTTGT"])

    def test_parseVariant_not_variant(self):
        variant_string = "##platypusOptions={'shareBAMFiles': 0, 'assembleBadReads': 1, 'bamFiles': ['GELS00000000010/Normal.bam'], 'assemblyRegionSize': 1500, 'minReads': 2, 'qualBinSize': 1, 'ploidy': 2, 'refFile': '/well/bsg/scratch/rimmer/Genomes/hs37d5.fa', 'maxHaplotypes': 50, 'filterVarsByCoverage': 1, 'maxSize': 1500, 'skipDifficultWindows': 0, 'parseNCBI': 0, 'skipRegionsFile': None, 'noCycles': 0, 'assembleAll': 1, 'minPosterior': 5, 'filterDuplicates': 1, 'abThreshold': 0.001, 'minFlank': 10, 'bufferSize': 100000, 'useEMLikelihoods': 0, 'logFileName': 'LogFiles/GELS00000000010_10.log', 'nCPU': 1, 'filterReadsWithUnmappedMates': 1, 'qdThreshold': 10, 'maxVariants': 8, 'scThreshold': 0.95, 'filterReadsWithDistantMates': 1, 'maxReads': 5000000, 'badReadsWindow': 11, 'genIndels': 1, 'useIndelErrorModel': 0, 'minMapQual': 20, 'maxVarDist': 15, 'printEvery': 10, 'maxGOF': 30, 'minVarDist': 9, 'rlen': 100, 'minGoodQualBases': 20, 'maxEMIterations': 100, 'filterReadPairsWithSmallInserts': 1, 'nInd': 1, ..."
        variant_record = self.vcf_parser.parseVariant(variant_string)
        
        self.assertIsNone(variant_record)
        
    def test_parseVcfFile_no_variants(self):
        vcf_records = self.vcf_parser.parseVcfFile("vcf_parser_test_no_variants.vcf")
        self.assertListEqual(vcf_records, [])
        
    def test_parseVcfFile_with_variants(self):
        vcf_records = self.vcf_parser.parseVcfFile("vcf_parser_tests_with_variants.vcf")
        
        self.assertEquals(len(vcf_records), 2)
        
        vcf_record = vcf_records[0]
        self.assertIsNotNone(vcf_record)
        self.assertListEqual(vcf_record.gene, ["ENSG00000101298"])
        self.assertEquals(vcf_record.ref, "C")
        self.assertListEqual(vcf_record.alt, ["CTG"])
        
        vcf_record = vcf_records[1]
        self.assertIsNotNone(vcf_record)
        self.assertListEqual(vcf_record.gene, ["ENSG00000101298"])
        self.assertEquals(vcf_record.ref, "TTG")
        self.assertListEqual(vcf_record.alt, ["T", "TTGTGTG"])
        
    def suite(self):
        return unittest.TestLoader().loadTestsFromTestCase(VcfParserTest)   
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()