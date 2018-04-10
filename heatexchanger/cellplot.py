from matplotlib.pyplot import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import imp
import layer

def nsf(num, n=1):
    """n-Significant Figures"""
    numstr = ("{0:.%ie}" % (n - 1)).format(num)
    return float(numstr)


def plot_cells(m, Z, cm=cm.RdBu_r, verbosity=0, zscale=None, zoff=None):
    "Plots a given array for every heat-exchange cell"
    f, a = subplots(figsize=(12, 12))
    sol = m.solution
    Nhotpipes = m.Nhotpipes
    Ncoldpipes = m.Ncoldpipes
    Zmin, Zmax = Z.min(), Z.max()
    dpos = 0
    for i in range(Nhotpipes):
        wpos = 0
        for j in range(Ncoldpipes):
            wcell = sol(m.coldpipes.w)[j].magnitude
            dcell = sol(m.hotpipes.w)[i].magnitude
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
    a.set_xlabel("width traveled by cold fluid [cm]")
    a.set_ylabel("depth traveled by hot fluid [cm]")
    return f, a


def hist_cells(m, Z, cm=cm.RdBu_r, verbosity=0, zscale=None, zoff=None):
    "Plots a given array for every heat-exchange cell"
    sol = m.solution
    Nhotpipes = m.Nhotpipes
    Ncoldpipes = m.Ncoldpipes
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
    a.set_xlabel("width traveled by cold fluid [cm]")
    a.set_ylabel("depth traveled by hot fluid [cm]")
    return f, arun

def gen_plots(m, sol, Ncld, Nhot):
    f, a = plot_cells(m, sol(m.c.T_hot), cm=cm.Reds,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Hot fluid temperature [K]")
    f.savefig("plots/T_hot.png")

    # Air temperature
    f, a = plot_cells(m, sol(m.c.T_cld), cm=cm.Blues,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Cold fluid temperature [K]")
    f.savefig("plots/T_cld.png")

    # Heat transfer
    Q = sol(m.Q).magnitude
    f, a = plot_cells(m, sol(m.c.dQ), cm=cm.Reds,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Heat transfer (%.2f Watts total)" % Q)
    f.savefig("plots/dQ.png")

    # Water velocity in each cell
    f, a = plot_cells(m, sol(m.hotpipes.v_avg), cm=cm.Blues,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Average velocity in hot cell (m/s)")
    f.savefig("plots/v_hot.png")

    # Air velocity in each cell
    f, a = plot_cells(m, sol(m.coldpipes.v_avg).T, cm=cm.Reds,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Average velocity in cold cell (m/s)")
    f.savefig("plots/v_cld.png")

    # Wall temperature
    f, a = plot_cells(m, sol(m.c.T_r), cm=cm.Reds,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Mean wall temperature (K)")
    f.savefig("plots/Tr.png")

    # Hot cell height
    f, a = plot_cells(m, sol(m.c.z_hot), cm=cm.Reds,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Hot cell height (cm)")
    f.savefig("plots/z_hot.png")

    # Cold cell height
    f, a = plot_cells(m, sol(m.c.z_cld), cm=cm.Blues,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Cold cell height (cm)")
    f.savefig("plots/z_cld.png")

    # Total cell height
    f, a = plot_cells(m, sol(m.c.z_cld)+sol(m.c.z_hot)+sol(m.c.t_plate), cm=cm.Blues,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Total cell height (cm)")
    f.savefig("plots/z.png")

    # Hot fin thickness
    f, a = plot_cells(m, sol(m.c.t_hot), cm=cm.Reds,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Hot fin thickness (cm)")
    f.savefig("plots/t_hot.png")

    # Cold fin thickness
    f, a = plot_cells(m, sol(m.c.t_cld), cm=cm.Blues,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Cold fin thickness (cm)")
    f.savefig("plots/t_cld.png")

    # Velocity area (to confirm mass flow rates) on hot side
    f, a = plot_cells(m, sol(m.hotpipes.v)*sol(m.hotpipes.A), cm=cm.Reds,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Velocity*area of hot flow (cm)")
    f.savefig("plots/vA_hot.png")

    # Velocity area (to confirm mass flow rates) on cold side
    f, a = plot_cells(m,(sol(m.coldpipes.v)*sol(m.coldpipes.A)).T, cm=cm.Reds,
                      zscale=1 / Nhot / Ncld, zoff=0.625 / Nhot / Ncld, verbosity=2)
    a.set_title("Velocity*area of cold flow (m^3/s)")
    f.savefig("plots/vA_cld.png")

