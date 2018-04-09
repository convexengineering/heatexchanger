import layer
import cellplot
import imp
from matplotlib.pyplot import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from gpkit import Model, units
from relaxed_constants import relaxed_constants
from gpkit.constraints.bounded import Bounded

from cellplot import gen_plots
from writetotext import genHXData

# Initializing SP single-layer HX model
imp.reload(layer)
Nw, Na = 5, 5
m = layer.Layer(Na, Nw)

# Model input parameters
m.substitutions.update({
                        # Geometric parameters,
                        m.n_fins:      5.,
                        m.x_dim:       5.*units('cm'),
                        m.y_dim:       10.*units('cm'),
                        m.z_dim:       1.*units('cm')
                        })

# Objective function
m.cost = (m.D_air+m.D_wat)/m.Q
#m = Model(m.cost,Bounded(m))
#m = relaxed_constants(m)

# Solving HX problem
sol = m.localsolve(verbosity=2)
#post_process(sol)
print sol('Q')

# Generating 2D plots
gen_plots(m, sol, Nw, Na)

# Writing complete solution file sol.txt
with open("sol.txt", "w") as f:
    f.write(sol.table())

# Generating ESP data file
genHXData(m, sol)
