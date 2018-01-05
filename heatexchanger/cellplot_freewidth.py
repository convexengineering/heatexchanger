from matplotlib.pyplot import *


def plot_cells(m, Z, cm=cm.RdBu_r, verbosity=0, zscale=None, zoff=None):
    "Plots a given array for every heat-exchange cell"
    f, a = subplots(figsize=(12, 12))
    sol = m.solution
    Zmin, Zmax = Z.min(), Z.max()
    dpos = 0
    for i in range(m.Nwaterpipes):
        wpos = 0
        for j in range(m.Nairpipes):
            wcell = sol(m.airpipes.w[j, i]).magnitude
            dcell = sol(m.waterpipes.w[i, j]).magnitude
            if not zscale:
                z = ((Z[j, i] - Zmin)/(Zmax-Zmin)).magnitude
            else:
                z = (Z[j, i].magnitude - zoff)/zscale
            r = Rectangle((dpos, wpos), dcell, wcell, facecolor=cm(z))
            if verbosity:
                labelx, labely = dpos, wpos + wcell/2
                label = "%.3g" % Z[j, i].magnitude
                if verbosity > 1:
                    label = "[%i,%i] : " % (i, j) + label
                else:
                    labelx += 0.4*dcell
                a.text(labelx, labely, label)
            wpos += sol(m.airpipes.w).magnitude[j, :].max()
            a.add_patch(r)
        dpos += sol(m.waterpipes.w).magnitude[i, :].max()
    ylim([0, wpos])
    xlim([0, dpos])
    a.set_frame_on(False)
    a.set_xlabel("width traveled by air [m]")
    a.set_ylabel("depth traveled by water [m]")
    return f, a


if __name__ == "__main__":
    from layer_freewidth import Layer
    Nw, Na = 5, 5
    m = Layer(Nw, Na)
    m.cost = 1/m.Q
    sol = m.localsolve()
    print sol(m.c.dQ).min(), 0.625/Nw/Na

    f, a = plot_cells(m, sol(m.c.T_hot), verbosity=0)
    a.set_title("Liquid Temperature [K]")
    f.savefig("T_liq.png")

    Q = sol(m.Q).magnitude
    f, a = plot_cells(m, sol(m.c.dQ), cm=cm.Reds, zscale=1/Nw/Na, zoff=0.625/Nw/Na)
    a.set_title("Heat Transfer (%.2f Watts total)" % Q)
    f.savefig("dQ.png")
