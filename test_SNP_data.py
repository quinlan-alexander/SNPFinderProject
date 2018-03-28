import unittest
import arv

from find_snp import SNPediaData
from find_snp import Twenty3andMeData
from find_snp import MyVariantInfo
from find_snp import SNPFinder

class TestFindSNPs(unittest.TestCase):

    def setUp(self):
        self.snp_data = SNPediaData()
        self.twenty3_and_me_data = Twenty3andMeData('MTHFR')
        self.myvariant_info = MyVariantInfo()
        genome = arv.load('genome_Barry_Alexander_v3_Full_20180311122349.txt')
        self.snpFinder = SNPFinder('foo', genome, 'MTHFR')

    def test_gets_SNP_data(self):
        self.snp_data.get_data('rs1801133')
        self.assertEqual('(C;C)', self.snp_data.get_common_allele())

    def test_honors_orientation(self):
        self.snp_data.get_data('rs1801133')
        self.assertEqual('(G;G)', self.snp_data.get_common_allele('plus'))

    def test_orientations_match(self):
        self.snp_data.get_data('rs1801133')
        self.assertEqual('(C;C)', self.snp_data.get_common_allele('minus'))

    def test_gene_lookup(self):
        self.assertEqual("MTHFR", self.twenty3_and_me_data.get_gene_name())

    def test_gene_snp_collection(self):
        self.assertEqual(['rs12023469', 'rs15854'],self.twenty3_and_me_data.get_snp_collection()[0:2])

    def test_myvariant_data(self):
        self.myvariant_info.get_data('rs58991260')
        self.assertEqual('G',self.myvariant_info.get_common_allele())
        self.assertEqual('A', self.myvariant_info.get_alt_allele())

    def test_ethnicity_trait(self):
        self.assertEqual('light-skinned, European',self.snpFinder.subject_ethnicity())
