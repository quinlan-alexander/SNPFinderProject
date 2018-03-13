import unittest

from find_snp import SNPediaData
from find_snp import Twenty3andMeData

class TestORFs(unittest.TestCase):

    def setUp(self):
        self.snp_data = SNPediaData()
        self.twenty3_and_me_data = Twenty3andMeData()

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
        self.assertEqual('{}', self.twenty3_and_me_data.get_gene('MTHFR'))
