import numpy as np
from scipy.interpolate import interp2d


def gencsm(m, sol, ID):
    nu = m.Ncoldpipes
    nhot = nu
    nv = 2
    nw = m.Nhotpipes
    ncold = nw
    nParams = 8

    # Creating corner coordinates
    x = [sum(sol(m.coldpipes.w)[0:i].to("m").magnitude) for i in range(nu+1)]
    y = [sum(sol(m.hotpipes.w)[0:i].to("m").magnitude) for i in range(nw+1)]
    xy = np.array([(x[i], y[j]) for j in range(nw+1) for i in range(nu+1)])
    xycent = np.array([[(x[i+1] + x[i])/2, (y[j+1]+y[j])/2]
                       for j in range(nw) for i in range(nu)])

    hxVals = [10*sol(m.cells.t_plate).to("m").magnitude,
              10*sol(m.cells.t_hot).to("m").magnitude,
              10*sol(m.cells.t_cld).to("m").magnitude,
              (sol(m.cells.z_hot)/sol(m.cells.z_cld)).magnitude]
    intlist = [interp2d(xycent[:, 0], xycent[:, 1], hxVals[i], kind='linear')
               for i in range(len(hxVals))]

    f = open('HX.csm', 'w')
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

# flow quantities
""")
    for name, val in m.design_parameters.items():
        if val in m.substitutions:
            val = m.substitutions[val]
        f.write("despmtr   %s   %s\n" % (name, val))

    f.write("""
# knot locations for tile placement
dimension uknots   1 %i 1
despmtr   uknots   %s

dimension vknots   1 3 1
despmtr   vknots   0.0;0.5;1.0

dimension wknots   1 %i 1
despmtr   wknots   %s
""" % (1+nu, ";".join(["%.2f" % (w/max(x)) for w in x]),
       1+nw, ";".join(["%.2f" % (w/max(y)) for w in y])))
    f.write("""
# duct definition (regular hexahedron)
dimension corners  8 3 0
set       corners  "0.0;     0.0;     0.0;     \\
                    0.0;     0.0;     z_width; \\
                    0.0;     y_width; 0.0;     \\
                    0.0;     y_width; z_width; \\
                    x_width; 0.0;     0.0;     \\
                    x_width; 0.0;     z_width; \\
                    x_width; y_width; 0.0;     \\
                    x_width; y_width; z_width;"

udparg    hex       corners   corners
udparg    hex       uknots    vknots     # u and v switched because
udparg    hex       vknots    uknots     # of way trivariate is made
udprim    hex       wknots    wknots
store duct

# tile the configuration
restore duct
udparg tile filename  $$/demo_tile.csm
udparg tile tablename <<
    %i   3   %i   4
%s
0.00   0.50   1.00
%s

""" % (1+nu, 1+nw,
       "   ".join(["%.2f" % (w/max(x)) for w in x]),
       "   ".join(["%.2f" % (w/max(y)) for w in y])))

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
dump HX_%03i.egads

restore duct
attribute _name $duct
""" % ID)

    f.close()
