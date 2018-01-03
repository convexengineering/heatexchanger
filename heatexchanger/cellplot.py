from matplotlib.pyplot import *


def plot_cells(m, Z, cm=cm.RdBu_r, verbosity=0):
    "Plots a given array for every heat-exchange cell"
    f, a = subplots(figsize=(12, 12))
    sol = m.solution
    Zmin, Zmax = Z.min(), Z.max()
    dpos = 0
    for i in range(m.Nwaterpipes):
        wpos = 0
        for j in range(m.Nairpipes):
            wcell = sol(m.airpipes.w)[j].magnitude
            dcell = sol(m.waterpipes.w)[i].magnitude
            z = ((Z[j, i] - Zmin)/(Zmax-Zmin)).magnitude
            r = Rectangle((dpos, wpos), dcell, wcell, facecolor=cm(z))
            if verbosity:
                labelx, labely = dpos, wpos + wcell/2
                label = "%.3g" % Z[j, i].magnitude
                if verbosity > 1:
                    label = "[%i,%i] : " % (i, j) + label
                else:
                    labelx += 0.4*dcell
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


if __name__ == "__main__":
    from layer import Layer
    m = Layer(10, 15)
    m.cost = 1/m.Q
    sol = m.localsolve()

    f, a = plot_cells(m, sol(m.waterpipes.T)[1:], verbosity=1)
    a.set_title("Liquid Temperature [K]")
    f.savefig("T_liq.png")

    f, a = plot_cells(m, sol(m.c.dQ), cm=cm.Reds)
    a.set_title("Heat Transfer [W]")
    f.savefig("dQ.png")
