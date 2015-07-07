__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Roller import Roller
import numpy as np
import pandas as pd

class tdRoller(Roller):
    """
    A roller that incorporates time delays
    """

