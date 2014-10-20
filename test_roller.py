import unittest
from Roller import Roller
import numpy as np

class TestRoller(unittest.TestCase):
    def setUp(self):
        self.roller = Roller('compressed_katrina_data.txt')

    def test_next(self):
        self.roller.next()
        window = self.roller.get_window()
        time_slice = window['time'].unique()
        correct_window = [1,2,3]
        self.assertTrue(np.array_equal(correct_window, time_slice))

if __name__ == '__main__':
    unittest.main()
