__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

class Window(object):
    def __init__(self, dataframe, alpha=None):
        self.dataframe = dataframe
        self.alpha = alpha
        self.data = dataframe.values

