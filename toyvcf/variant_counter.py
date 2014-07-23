'''
Created on Jul 22, 2014

This program takes a vcf filename as a parameter and for each variant with a PASS Filter value and an existing Gene prints:

1. The total number of variants (mutations) represented in the file together with a breakdown of this total by mutation type (SNV/INS/DEL/MNV).
2. A list of genes within which at least one variant has been found, the total number of variants found in that gene, and a breakdown by variant type.

Variants that do not have a Gene specified are excluded from this list!

@author: siakhnin
'''

import sys
import re

def main():
    ''' Main function that runs the program '''
    
    #Parse vcf file and extract mutations.
    vcf_parser = VcfParser()
    vcf_record_list = vcf_parser.parseVcfFile(sys.argv[1])
    
    #Count extracted mutations by gene and type.
    vcf_mutation_counter = VcfMutationCounter()
    vcf_mutation_counter.addMutationsToCount(vcf_record_list)
    
    #Print mutation counts by type
    vcf_mutation_counter.printGlobalMutationTally()
    
    #Print mutation types by gene and by type
    vcf_mutation_counter.printMutationsByGene()
    
class VcfRecord:
    '''
    Class VcfRecord represents one variant as described in a vcf file. 
    
        Currently supported fields are:
            ref - Reference allele
            alt - Alternative alleles
            gene - Gene IDs overlapping the variant
    '''
    ref = None
    alt = []
    gene = []
    
    def __init__(self, ref, alt, gene):
        self.ref = ref
        self.alt = alt
        self.gene = gene
    
    def isSNV(self):
        ''' Returns True if reference and all alternative bases are 1 nucleotide long and are different. '''
            
        for alt_bases in self.alt:
            if not(len(self.ref) == len(alt_bases) == 1 and self.ref != alt_bases):
                return False
        return True
    
    
    def isInsertion(self):
        '''
            Returns True if length of all alternative bases is greater than the 
            length of the reference bases and reference is a left-aligned substring of the alternative
        '''
        for alt_bases in self.alt:
            if not(len(self.ref) < len(alt_bases) and re.match(self.ref, alt_bases)):
                return False
        return True

    def isDeletion(self):
        '''
            Returns True if length of all alternative bases is smaller than the 
            length of the reference bases and alternative is a left-aligned substring of the reference
        '''
    
        for alt_bases in self.alt:
            if not(len(self.ref) > len(alt_bases) and re.match(alt_bases, self.ref)):
                return False
        return True
    
    def getMutationType(self):
        '''
            Returns which mutation class this VcfRecord encodes
        '''    
    
        if self.isSNV():
            return "SNV"
        elif self.isInsertion():
            return "INS"
        elif self.isDeletion():
            return "DEL"
        else:
            return "MNV"


class VcfParser:
    '''
        Class VcfParser represents a vcf file parser. 
        It can parse a vcf file and extract information into VcfRecords.
    '''
    
    def __init__(self):
        ''' Initialize the parser with a parse pattern for vcf entries. '''
        self.parse_pattern = re.compile(r"""
            .*?                          #match Chromosome
            \s.[0-9]*                    #match Position
            \s.*?                        #match ID
            \s(?P<ref>[A|C|G|T|N]+)      #capture REF bases
            \s(?P<alt>.*?)               #capture ALT bases
            \s.*?PASS                    #match FILTER PASS 
            \s.*?Gene=(?P<gene>.*?);     #capture INFO Gene name
            """, re.VERBOSE)
    
    def parseVcfFile(self, filename):
        ''' Returns a list of VcfRecords parsed from the file named filename. '''    
    
        input_vcf = open(filename, "r")
        vcf_record_list = []
        
        for line in input_vcf:
            #Attempt to parse the current line as a VcfRecord
            vcf_record = self.parseVariant(line)
            
            #If current line was successfully parsed add it to the list of vcf records
            if vcf_record != None:
                vcf_record_list.append(vcf_record)
                
        return vcf_record_list
    
    
    def parseVariant(self, variant_string):
        '''
            Returns a VcfRecord represented by the variant_string passed in as a parameter
            or None if variant_string does not represent a VcfRecord.
            Currently supported fields are:
                ref - Reference allele
                alt - Alternative alleles
                gene - Gene IDs overlapping the variant
            
        '''
        ref = alt = gene = None
        my_match = self.parse_pattern.match(variant_string)
        
        #If match was successful create and return the VcfRecord using extracted values
        if my_match != None:
            ref = my_match.group("ref")
            alt = my_match.group("alt").split(",")
            gene = my_match.group("gene").split(",")
        
            return VcfRecord(ref, alt, gene)
    
class VcfMutationCounter:
    '''
        Class VcfMutationCounter maintains counts of mutations by type and gene ID.
    '''
    # Declare mutation types here to ensure consistent naming.
    snv_key = "SNV"
    ins_key = "INS"
    del_key = "DEL"
    mnv_key = "MNV"
    
    
    def __init__(self):
        self.global_mutation_tallies = {self.snv_key:0, self.ins_key:0, self.del_key:0, self.mnv_key:0}
        self.mutations_by_gene = {}
    
   
    def addMutationsToCount(self, vcf_record_list):
        '''
            Takes a list of VcfRecords as a parameter and adds the encoded mutations
            to the lists of mutation counts being maintained by this VcfMutationCounter.
        '''
        for vcf_record in vcf_record_list:
            
            #Determine the type of this mutation
            mutation_class = vcf_record.getMutationType()
            
            #Increment the count of mutations for the appropriate type
            self.global_mutation_tallies[mutation_class] = self.global_mutation_tallies.setdefault(mutation_class,0) + 1
            
            #For all genes overlapping this mutation, increment the count of appropriate 
            #mutations for that gene.
            for gene_name in vcf_record.gene:
                
                #If this gene had no prior mutations recorded, initialize a list of mutation counts
                if self.mutations_by_gene.get(gene_name) == None:
                    self.mutations_by_gene[gene_name] = {self.snv_key:0, self.ins_key:0, self.del_key:0, self.mnv_key:0}
                
                gene_mutation_tallies = self.mutations_by_gene.setdefault(gene_name,{})
                gene_mutation_tallies[mutation_class] = gene_mutation_tallies.setdefault(mutation_class,0) + 1
                
    def printGlobalMutationTally(self):
        '''
            Prints the tally of mutations for this VcfMutationCounter by mutation type
        '''
        print "\n============== Total Variant Counts =============="
        print "All Types: {}".format(sum(self.global_mutation_tallies.values()))
        for mutation_type in self.global_mutation_tallies.keys():
            print "{}s: {}".format(mutation_type, str(self.global_mutation_tallies[mutation_type]))
    
    def printMutationsByGene(self):
        '''
            Prints the tally of mutations for this VcfMutationCounter by gene ID and mutation type
        '''
        print "\n\n============== Variant Counts By Gene =============="
        format_template = "{:<20}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}"
        print format_template.format("Gene ID", "All Types", self.snv_key, self.ins_key, self.del_key, self.mnv_key)
        for gene_name in self.mutations_by_gene.keys():
            print format_template.format(gene_name, 
                                          sum(self.mutations_by_gene[gene_name].values()),
                                          self.mutations_by_gene[gene_name].get(self.snv_key),
                                          self.mutations_by_gene[gene_name].get(self.ins_key),
                                          self.mutations_by_gene[gene_name].get(self.del_key),
                                          self.mutations_by_gene[gene_name].get(self.mnv_key))
           
if __name__ == '__main__':
    main()