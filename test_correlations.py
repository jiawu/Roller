__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import u_functions as ufun
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys

t = np.linspace(0,2*np.pi,10000)
s = np.cos(t)

phases = np.linspace(0,2*np.pi,100)
x = ufun.simple(t).signal()
print np.cov(x,x)
z = ufun.simple(t, phase=0.5*np.pi).signal()

#plt.plot(t,x,t,z)
#plt.show()
#plt.plot(x,z)
#plt.show()
pearsons = []
spearmans = []
for phase in phases:
    y = ufun.simple(t, phase=phase).signal()
    coef, p = stats.pearsonr(x,y)
    pearsons.append(coef)
    coef, p = stats.spearmanr(x,y)
    spearmans.append(coef)

plt.plot(phases,pearsons, phases, spearmans)
plt.show()

freqs = np.linspace(1,10,100)
pearsons2 = []
spearmans2 = []

for freq in freqs:
    w = ufun.simple(t,freq=freq).signal()
    coef, p = stats.pearsonr(x,w)
    pearsons2.append(coef)
    coef, p = stats.spearmanr(x,w)
    spearmans2.append(coef)

plt.plot(freqs,pearsons2, freqs, spearmans2)
plt.show()

a1=1
b1=1
c1=0
d1=0
a2=a1
b2=b1
c2=c1
d2=d1

cov = a1*np.sin(b1*t+c1)*a2*np.sin(b2*t+c2)+a1*d2*np.sin(b1*t+c1)+d1*a2*np.sin(b2*t+c2)+d1*d2
print np.mean(cov)