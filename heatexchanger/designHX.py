import layer
import cellplot
import rectpipe
from materials import *

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
imp.reload(rectpipe)
Ncoldpipes, Nhotpipes = 4,4
coldFluid = Air()
hotFluid = Water()
material = StainlessSteel()
n_fins = 5
m = layer.Layer(Ncoldpipes, Nhotpipes, coldFluid, hotFluid, material)

# Model input parameters
m.substitutions.update({
                        # Geometric parameters,
                        m.n_fins:         n_fins,
                        m.max_porosity:   0.7,
                        m.x_dim:          5.*units('cm'),      # max length of cold flow
                        m.y_dim:          10.*units('cm'),     # max length of hot flow
                        m.z_dim:          1.*units('cm'),      # max height of layer

                        # Inlet flow parameters
                        m.coldpipes.v_in:    20.*np.ones(Ncoldpipes)*units('m/s'),
                        m.hotpipes.v_in:     1.*np.ones(Nhotpipes)*units('m/s'),
                        m.coldpipes.P_in:    101000*np.ones(Ncoldpipes)*units('Pa'),
                        m.coldpipes.P_out:   101000*np.ones(Ncoldpipes)*units('Pa'),
                        m.hotpipes.P_in:     101000*np.ones(Nhotpipes)*units('Pa'),
                        m.hotpipes.P_out:    101000*np.ones(Nhotpipes)*units('Pa'),

                        # Heat transfer parameters
                        #m.T_min_cold:         300.*units('K'),
                        m.T_max_hot:           450.*units('K'),
                        m.hotpipes.T_in:       500.*np.ones(Nhotpipes)*units('K'),
                        m.coldpipes.T_in:      303.*np.ones(Ncoldpipes)*units('K'),
                        })

# Objective function
m.cost = (m.D_hot+m.D_cold)/m.Q

#m = Model(m.cost,Bounded(m))
#m = relaxed_constants(m)

# Solving HX problem
sol = m.localsolve(verbosity=2)
#post_process(sol)
print sol('Q')

# Generating 2D plots
gen_plots(m, sol, Ncoldpipes, Nhotpipes)

# Writing complete solution file sol.txt
with open("sol.txt", "w") as f:
    f.write(sol.table())

# Generating ESP data file
genHXData(m, sol)
