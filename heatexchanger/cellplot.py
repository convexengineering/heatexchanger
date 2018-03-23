from matplotlib.pyplot import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def nsf(num, n=1):
    """n-Significant Figures"""
    numstr = ("{0:.%ie}" % (n - 1)).format(num)
    return float(numstr)


def plot_cells(m, Z, cm=cm.RdBu_r, verbosity=0, zscale=None, zoff=None):
    "Plots a given array for every heat-exchange cell"
    f, a = subplots(figsize=(12, 12))
    sol = m.solution
    Nwaterpipes = m.Nwaterpipes
    Nairpipes = m.Nairpipes
    Zmin, Zmax = Z.min(), Z.max()
    dpos = 0
    for i in range(Nwaterpipes):
        wpos = 0
        for j in range(Nairpipes):
            wcell = sol(m.airpipes.w)[j].magnitude
            dcell = sol(m.waterpipes.w)[i].magnitude
            if not zscale:
                z = ((Z[j, i] - Zmin) / (Zmax - Zmin)).magnitude
            else:
                z = (Z[j, i].magnitude - zoff) / zscale
            r = Rectangle((dpos, wpos), dcell, wcell, facecolor=cm(z))
            if verbosity:
                labelx, labely = dpos, wpos + wcell / 2
                label = "%.3g" % Z[j, i].magnitude
                if verbosity > 1:
                    label = "[%i,%i] : " % (i, j) + label
                else:
                    labelx += 0.4 * dcell
                a.text(labelx, labely, label)
            wpos += wcell
            a.add_patch(r)
        dpos += dcell
    ylim([0, wpos])
    xlim([0, dpos])
    a.set_frame_on(False)
    a.set_xlabel("width traveled by air [m]")
    a.set_ylabel("depth traveled by water [m]")
    return f, a


def hist_cells(m, Z, cm=cm.RdBu_r, verbosity=0, zscale=None, zoff=None):
    "Plots a given array for every heat-exchange cell"
    sol = m.solution
    Nwaterpipes = m.Nwaterpipes
    Nairpipes = m.Nairpipes
    Zmin, Zmax = Z.min(), Z.max()
    dpos = 0

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # x =

    # Construct arrays for the anchor positions of the 16 bars.
    # Note: np.meshgrid gives arrays in (ny, nx) so we use 'F' to flatten xpos,
    # ypos in column-major order. For numpy >= 1.7, we could instead call meshgrid
    # with indexing='ij'.
    xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25)
    xpos = xpos.flatten('F')
    ypos = ypos.flatten('F')
    zpos = np.zeros_like(xpos)

    ylim([0, wpos])
    xlim([0, dpos])
    a.set_frame_on(False)
    a.set_xlabel("width traveled by air [m]")
    a.set_ylabel("depth traveled by water [m]")
    return f, arun


if __name__ == "__main__":
    from layer import Layer
    Nw, Na = 5, 5
    m = Layer(Nw, Na)
    # m.substitutions.update({
    #     'V_tot':1*units('cm^3'),
    #     'Q'    :4*units('W')
    #     })
    penalties = (m.waterpipes.dP_scale.prod()*m.airpipes.dP_scale.prod()*m.waterpipes.dT.prod()*m.airpipes.dT.prod())**-1
    m.cost = penalties*1*m.Q**-1*(1*m.waterpipes.D.sum()+ 1*m.airpipes.D.sum())
    #m = Model(m.cost,Bounded(m))
    #m = relaxed_constants(m)
    sol = m.localsolve(verbosity=4)
    #post_process(sol)
    print sol('Q')

    # Liquid temperature
    f, a = plot_cells(m, sol(m.c.T_hot), cm=cm.Reds,
                      zscale=1 / Nw / Na, zoff=0.625 / Nw / Na, verbosity=2)
    a.set_title("Liquid Temperature [K]")
    f.savefig("plots/T_liq.png")

    # Air temperature
    f, a = plot_cells(m, sol(m.c.T_cld), cm=cm.Blues,
                      zscale=1 / Nw / Na, zoff=0.625 / Nw / Na, verbosity=2)
    a.set_title("Air Temperature [K]")
    f.savefig("plots/T_air.png")

    # Heat transfer
    Q = sol(m.Q).magnitude
    f, a = plot_cells(m, sol(m.c.dQ), cm=cm.Reds,
                      zscale=1 / Nw / Na, zoff=0.625 / Nw / Na, verbosity=2)
    a.set_title("Heat Transfer (%.2f Watts total)" % Q)
    f.savefig("plots/dQ.png")

    # Water velocity in each cell
    f, a = plot_cells(m, sol(m.waterpipes.v_avg), cm=cm.Blues,
                      zscale=1 / Nw / Na, zoff=0.625 / Nw / Na, verbosity=2)
    a.set_title("Average velocity in water cell")
    f.savefig("plots/waterV.png")

    # Air velocity in each cell
    f, a = plot_cells(m, sol(m.airpipes.v_avg).T, cm=cm.Reds,
                      zscale=1 / Nw / Na, zoff=0.625 / Nw / Na, verbosity=2)
    a.set_title("Average velocity in air cell")
    f.savefig("plots/airV.png")

    # Wall temperature
    f, a = plot_cells(m, sol(m.c.T_r), cm=cm.Reds,
                      zscale=1 / Nw / Na, zoff=0.625 / Nw / Na, verbosity=2)
    a.set_title("Mean wall temperature (K)")
    f.savefig("plots/Tr.png")
