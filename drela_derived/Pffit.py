import pandas as pd
import numpy as np
from numpy import logspace, log, log10
import matplotlib.pyplot as plt
from gpfit.fit import fit


# Fitting for Pf, pressure drop parameter
df = pd.read_csv('PfFit.csv')

Re = np.array(df['Re'])
Pf = np.array(df['Pf'])
Remax = np.amax(Re)
Pfmax = np.amax(Pf)

logRe = log(Re/Remax)
logPf = log(Pf/Pfmax)

K=2
Type = 'SMA'
cstrt, rms_error = fit(logRe,logPf,K,Type)

