import pandas as pd
import numpy as np
from numpy import logspace, log, log10
import matplotlib.pyplot as plt
from gpfit.fit import fit


# Fitting for etah, effectiveness parameter
df = pd.read_csv('etahFit.csv')

Re = np.array(df['Re'])
etah = np.array(df['etah'])
Remax = np.amax(Re)
etahmax = np.amax(etah)

logRe = log(Re/Remax)
logetah = log(etah/etahmax)

K=1
Type = 'SMA'
cstrt, rms_error = fit(logRe,logetah,K,Type)

plt.plot(log(Re), log((0.799*(Re/Remax)**-0.0296)**10*etahmax))
plt.plot(log(Re),log(etah))
plt.show()