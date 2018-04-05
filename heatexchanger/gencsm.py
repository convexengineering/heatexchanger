import numpy as np
from scipy.interpolate import interp2d


def gencsm(m, sol):
    nu = m.Nwaterpipes
    nwater = nu
    nv = 1
    nw = m.Nairpipes
    nair = nw
    nParams = 8

    # Creating corner coordinates
    x = [sum(sol(m.waterpipes.w)[0:i].magnitude) for i in range(nu+1)]
    y = [sum(sol(m.airpipes.w)[0:i].magnitude) for i in range(nw+1)]
    xy = np.array([(x[i], y[j]) for j in range(nw+1) for i in range(nu+1)])
    xycent = np.array([[(x[i+1] + x[i])/2, (y[j+1]+y[j])/2]
                       for j in range(nw) for i in range(nu)])
    z = sol(m.c.z_hot) + sol(m.c.z_cld) + 2*sol(m.c.t_plate)

    hxVals = [sol(m.c.t_plate).magnitude,
              sol(m.c.t_hot).magnitude,
              sol(m.c.t_cld).magnitude,
              sol(m.c.z_hot).magnitude/sol(m.c.z_cld).magnitude]
    intlist = [interp2d(xycent[:, 0], xycent[:, 1], hxVals[i], kind='linear')
               for i in range(len(hxVals))]

    f = open('HX.csm', 'w')
    # TODO: where to put Nair, Nwater?
    f.write("""# HeateXchanger
# autogenerated CSM file

#           ^  y,u
#           |
#           +-----------------------+
#          /:                      /|
#         / :                     / |
#        +-----------------------+  |
#  ----> |  :                    |  |
#  cold  |  + . . . . . . . . . . . + -> x,w
#  ----> | '                     | /
#        |'                      |/
#        +-----------------------+
#       /         ^   ^
#      v  z,v     |hot|

# knot locations for tile placement
dimension uknots   1 7 1
despmtr   uknots   0.00;0.10;0.22;0.36;0.54;0.74;1.00

dimension vknots   1 3 1
despmtr   vknots   0.0;0.5;1.0

dimension wknots   1 5 1
despmtr   wknots   0.00;0.14;0.34;0.62;1.00

# flow quantities
""")
    for name, val in m.design_parameters.items():
        if val in m.substitutions:
            val = m.substitutions[val]
        f.write("despmtr   %s   %s\n" % (name, val))
    f.write("""
# duct definition (regular hexahedron)
dimension corners  8 3 0
set       corners  "0.0;   0.0;   0.0;   \\
                    0.0;   0.0;   z_dim; \\
                    0.0;   y_dim; 0.0;   \\
                    0.0;   y_dim; z_dim; \\
                    x_dim; 0.0;   0.0;   \\
                    x_dim; 0.0;   z_dim; \\
                    x_dim; y_dim; 0.0;   \\
                    x_dim; y_dim; z_dim;"

udparg    hex       corners   corners
udparg    hex       uknots    vknots     # u and v switched because
udparg    hex       vknots    uknots     # of way trivariate is made
udprim    hex       wknots    wknots
store duct

# tile the configuration
restore duct
udparg tile filename  $$/demo_tile.csm
udparg tile tablename <<
    7   3   5   4
0.00   0.10   0.22   0.36   0.54   0.74   1.00
0.00   0.50   1.00
0.00   0.14   0.34   0.62   1.00

""")

    f.write('thkPlate' + '\n')
    f.write('v' + '\n')
    f.write('thkHot' + '\n')
    f.write('w' + '\n')
    f.write('thkCold' + '\n')
    f.write('u' + '\n')
    f.write('hot2cold' + '\n')
    f.write('.' + '\n\n')

    for w in range(nw+1):
        for v in range(nv+1):
            for u in range(nu+1):
                for k in range(len(intlist)):
                    f.write("%.4f" % intlist[k](xy[u+w*u, 0], xy[u+w*u, 1])[0])
                    f.write(' ')
                f.write('\n')
            f.write('\n')

    f.write(""">>
udparg tile nutile    1
udparg tile nvtile    1
udparg tile nwtile    1
udprim tile

assert @@numbodys 1

restore duct
attribute _name $duct
""")

    f.close()
