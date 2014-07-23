This program takes a vcf filename as a parameter and for each variant with a PASS Filter value and an existing Gene prints:

1. The total number of variants (mutations) represented in the file together with a breakdown of this total by mutation type (SNV/INS/DEL/MNV).
2. A list of genes within which at least one variant has been found, the total number of variants found in that gene, and a breakdown by variant type.

Variants that do not have a Gene specified are excluded from this list!


Running the program
--------------------

python variant_counter.py testData.vcf


Running the unit tests
-----------------------

python -m unittest discover -s <path_to_tests> -p '*_tests.py'



Written by Sergei Iakhnin on July 22, 2014
