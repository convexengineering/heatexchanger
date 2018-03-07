from matplotlib.pyplot import *

def nsf(num, n=1):
    """n-Significant Figures"""
    numstr = ("{0:.%ie}" % (n-1)).format(num)
    return float(numstr)


def plot_cells(m, Z, cm=cm.RdBu_r, verbosity=0, zscale=None, zoff=None):
    "Plots a given array for every heat-exchange cell"
    f, a = subplots(figsize=(12, 12))
    sol = m.solution
    Nwaterpipes = m.original.Nwaterpipes
    Nairpipes = m.original.Nairpipes
    Zmin, Zmax = Z.min(), Z.max()
    dpos = 0
    for i in range(Nwaterpipes):
        wpos = 0
        for j in range(Nairpipes):
            wcell = sol(m.original.airpipes.w)[j].magnitude
            dcell = sol(m.original.waterpipes.w)[i].magnitude
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
    Nw, Na = 5, 5
    print sol(m.original.c.dQ).min(), 0.625/Nw/Na

    f, a = plot_cells(m, sol(m.original.c.T_hot), cm=cm.Reds, zscale=1/Nw/Na, zoff=0.625/Nw/Na, verbosity=2)
    a.set_title("Liquid Temperature [K]")
    f.savefig("plots/T_liq.png")

    Q = sol(m.original.Q).magnitude
    f, a = plot_cells(m, sol(m.original.c.dQ), cm=cm.Reds, zscale=1/Nw/Na, zoff=0.625/Nw/Na, verbosity=2)
    a.set_title("Heat Transfer (%.2f Watts total)" % Q)
    f.savefig("plots/dQ.png")

    # Water drag in each cell
    waterD = sum(sum(sol(m.original.waterpipes.D_seg).magnitude))
    f, a = plot_cells(m, sol(m.original.waterpipes.D_seg), cm=cm.Blues, zscale=1/Nw/Na, zoff=0.625/Nw/Na, verbosity=2)
    a.set_title("Drag force due to each water cell (%.2f Watts total)" % waterD)
    f.savefig("plots/waterD.png")

    # Air drag in each cell
    airD = sum(sum(sol(m.original.airpipes.D_seg).magnitude))
    f, a = plot_cells(m, sol(m.original.airpipes.D_seg), cm=cm.Reds, zscale=1/Nw/Na, zoff=0.625/Nw/Na, verbosity=2)
    a.set_title("Drag force due to each air cell (%.2f N total)" % airD)
    f.savefig("plots/airD.png")

    # Wall temperature
    f, a = plot_cells(m, sol(m.original.c.T_r), cm=cm.Reds, zscale=1/Nw/Na, zoff=0.625/Nw/Na, verbosity=2)
    a.set_title("Mean wall temperature (K)")
    f.savefig("plots/Tr.png")




