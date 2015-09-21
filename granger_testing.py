__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import pandas as pd
import sys
import statsmodels as sm
import numpy as np
from statsmodels.tsa.api import VAR
from statsmodels.tsa.base.datetools import dates_from_str

# some example data
mdata = sm.datasets.macrodata.load_pandas().data

# prepare the dates index
dates = mdata[['year', 'quarter']].astype(int).astype(str)
quarterly = dates["year"] + "Q" + dates["quarter"]

quarterly = dates_from_str(quarterly)

mdata = mdata[['realgdp','realcons','realinv']]
mdata.index = pd.DatetimeIndex(quarterly)
data = np.log(mdata).diff().dropna()

# make a VAR model
model = VAR(data)
print(model.select_order(15))
results = model.fit(3)

results.test_causality('realgdp',['realinv', 'realcons'], kind='f')