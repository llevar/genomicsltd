from unittest import TestSuite
from vcf_record_tests import VcfRecordTest
from vcf_parser_tests import VcfParserTest
from vcf_mutation_counter_tests import VcfMutationCounterTest

test_cases = (VcfRecordTest, VcfParserTest, VcfMutationCounterTest)

def load_tests(loader, tests, pattern):
    suite = TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite