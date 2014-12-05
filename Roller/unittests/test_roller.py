import unittest
import Roller
import numpy as np

class TestRoller(unittest.TestCase):
    def setUp(self):
        self.roller = Roller.Roller('../../data/emt/compressed_katrina_data.txt', 5, None, "time", " ")

    def test_next(self):
        self.roller.next()
        window = self.roller.get_window_raw()
        time_slice = window['time'].unique()
        correct_window = [1,2,3]
        self.assertTrue(np.array_equal(correct_window, time_slice))

    def test_get_only_genes(self):
        only_genes = self.roller.get_window()
        header = only_genes.columns.values
        correct_header = ['AP1','AP2','AP3','AP4', 'AR', 'Bcat', 'Brachyury', 'cmyc', 'CRE', 'E2F','ELK1', 'ER', 'ETS1', 'FOXA', 'FOXO3A', 'GATA1', 'GATA2', 'GATA3','GLI', 'GR', 'HIF1', 'HNF1A', 'HOXA1', 'HSE', 'KLF4', 'LHX8','MEF2', 'MNX', 'MNX1', 'MYB', 'NANOG', 'NFAT', 'NFKb', 'NOBOX','Notch1', 'OCT', 'P53', 'PAX1', 'PEA3', 'PR', 'PTTG', 'RAR','RUNX1', 'RUNX2', 'SMAD1', 'SMAD3', 'SOX', 'SP1', 'SRF', 'STAT1','STAT3', 'STAT4', 'STAT5', 'VDR', 'WT1', 'YY1']
        self.assertTrue(np.array_equal(correct_header, header))

if __name__ == '__main__':
    unittest.main()
