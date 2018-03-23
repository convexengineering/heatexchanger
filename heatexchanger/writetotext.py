import numpy as np
from scipy.interpolate import interp2d

def genHXData(m,sol):
    nu = m.Nwaterpipes
    nwater = nu
    nv = 1
    nw = m.Nairpipes
    nair = nw
    nParams = 8

    # Creating corner coordinates
    x = [sum(sol(m.waterpipes.w)[0:i].magnitude) for i in range(nu+1)]
    y = [sum(sol(m.airpipes.w)[0:i].magnitude) for i in range(nw+1)]
    xy = np.array([(x[i],y[j]) for j in range(nw+1) for i in range(nu+1)])
    xycent = np.array([[(x[i+1]+ x[i])/2,(y[j+1]+y[j])/2] for j in range (nw) for i in range(nu)])
    z = sol(m.c.z_hot) + sol(m.c.z_cld) + 2*sol(m.c.t_plate)

    f = open('hxOut.txt','w')

    f.write(str(nu) + ' ' + str(nv) + ' ' + str(nw) + ' ' + str(nParams))
    f.write('\n')
    for i in range(len(x)):
        f.write(str(x[i]/max(x)) + ' ')
    f.write('\n')
    f.write('0. 1.')
    f.write('\n')
    for i in range(len(y)):
        f.write(str(y[i]/max(y)) + ' ')
    f.write('\n\n')

    hxVals = [sol(m.c.t_plate).magnitude,
               sol(m.c.t_hot).magnitude,
               sol(m.c.t_cld).magnitude,
               sol(m.c.z_hot).magnitude/sol(m.c.z_cld).magnitude]
    intlist = [interp2d(xycent[:,0],xycent[:,1],hxVals[i],kind='linear') for i in range(len(hxVals))]

    for w in range(nw+1):
        for v in range(nv+1):
            for u in range(nu+1):
                for k in range(len(intlist)):
                    f.write(str(intlist[k](xy[u+w*u,0],xy[u+w*u,1])[0]))
                    f.write('\n')
                f.write('\n')
    f.close()
