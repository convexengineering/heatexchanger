import layer
import cellplot
from materials import *

import imp
from matplotlib.pyplot import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from gpkit import Model, units
from relaxed_constants import relaxed_constants
from gpkit.constraints.bounded import Bounded

from cellplot import gen_plots

# Initializing SP single-layer HX model
imp.reload(layer)
Ncold, Nhot = 3, 3
m = layer.Layer(Ncold, Nhot)

# Model input parameters
m.substitutions.update({
                        # Geometric parameters,
                        m.n_fins:      5.,
                        m.x_dim:       5.*units('cm'),
                        m.y_dim:       10.*units('cm'),
                        m.z_dim:       1.*units('cm')
                        })

# Objective function
m.cost = 1/m.Q
#m = Model(m.cost,Bounded(m))
#m = relaxed_constants(m)

# Solving HX problem
sol = m.localsolve(verbosity=2)
#post_process(sol)
print sol('Q')

# Generating 2D plots
gen_plots(m, sol, Ncold, Nhot)

# Writing complete solution file sol.txt
with open("sol.txt", "w") as f:
    f.write(sol.table())
