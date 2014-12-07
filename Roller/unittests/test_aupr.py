import unittest
from Roller.util.Evaluator import Evaluator
import numpy as np
from random import randint
import numpy.testing as npt
import random
from sets import Set
import pdb

"""this test generally checks to see if aupr calculations are robust and if trivial cases evaluate as expected"""

#todo: this test should probably work with multiple gold standards instead of just the dream one. please incorporate test methods of different standards.

class TestAUPR(unittest.TestCase):
    def setUp(self):
        gold_standard_file = "../../data/dream4/insilico_size10_1_goldstandard.tsv"
        self.evaluator = Evaluator(gold_standard_file, sep = '\t')

    def test_perfect_pr(self):
        prediction = self.evaluator.gs_data
        precision, recall, aupr = self.evaluator.calc_pr(prediction)
        #lazy unit testing. check if first 14 are not 1.0s in recall, check if others are all 1.0s
        check = Set([1.0 for x in xrange(14)])
        not_ones =Set(recall[0:14]) - check
        self.assertEqual(len(not_ones),14)
        #one more lazy unit test. check if there is one point where P = R = 1.0
    def test_if_perfect_case(self):
        prediction = self.evaluator.gs_data
        precision, recall, aupr = self.evaluator.calc_pr(prediction)
        point = Set(precision).intersection(Set(recall))
        real_point = 1.0
        self.assertIn(real_point, point)

if __name__ == '__main__':
    unittest.main()
