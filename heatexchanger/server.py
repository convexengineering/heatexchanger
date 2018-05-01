from SimpleWebSocketServer import SimpleWebSocketServer, WebSocket
import json
from layer import Layer
from materials import *
from gencsm import gencsm
from shutil import copyfile
from designHX import designHX

EXIT = [False]
ID = 0
LASTSOL = [None]


def genfiles(m, sol):
    global ID
    gensoltxt(m, sol, ID)
    gencsm(m, sol, ID)
    copyfile("HX.csm", "HX_%03i.csm" % ID)
    ID += 1


def gensoltxt(m, sol, ID):
    with open("sol_%03i.txt" % ID, "w") as f:
        for var in sorted(m.varkeys, key=str):
            f.write("%s [%s]\t\t%f\n" % (var, var.unitstr(dimless="-"),
                                         sol["variables"][var]))


class HXGPServer(WebSocket):

    def handleMessage(self):
        print "< received", repr(self.data)
        try:
            self.data = json.loads(self.data)
            print self.data

            Ncoldpipes = self.data["Cold_Channels"]
            Nhotpipes = self.data["Hot_Channels"]
            if (Ncoldpipes, Nhotpipes) == LASTSOL[0][0]:
                x0 = LASTSOL[0][1]["variables"]
            else:
                x0 = None

            coldfluid_model = Air
            hotfluid_model = Water
            material_model = StainlessSteel
            m = designHX(Ncoldpipes, Nhotpipes, coldfluid_model, hotfluid_model, material_model)
            m.cost = 1/m.Q

            for name, value in self.data.items():
                try:
                    key = m.design_parameters[name]
                    m.substitutions[key] = value
                except KeyError as e:
                    print repr(e)

            sol = m.localsolve(x0=x0)
            LASTSOL[0] = ((Ncoldpipes, Nhotpipes), sol)
            genfiles(m, sol)

            self.send({"status": "optimal",
                       "msg": ("Successfully optimized."
                               " Optimal heat transfer: %.1f watts "
                               % sol["variables"][m.Q])})
        except Exception as e:
            self.send({"status": "unknown", "msg": "The last solution"
                      " raised an exception; tweak it and send again."})
            print type(e), e

    def send(self, msg):
        print "> sent", repr(msg)
        self.sendMessage(unicode(json.dumps(msg)))

    def handleConnected(self):
        print self.address, "connected"

    def handleClose(self):
        print self.address, "closed"
        EXIT[0] = True


if __name__ == "__main__":
    # TODO: uncomment to produce the initial CSM file before serving
    Ncoldpipes, Nhotpipes = 3,3
    coldfluid_model = Air
    hotfluid_model = Water
    material_model = StainlessSteel
    m = designHX(Ncoldpipes, Nhotpipes, coldfluid_model, hotfluid_model, material_model)
    sol = m.localsolve(verbosity=2)
    LASTSOL[0] = ((m.Ncoldpipes,m.Nhotpipes), sol)
    genfiles(m, sol)
    server = SimpleWebSocketServer('', 8000, HXGPServer)
    while not EXIT[0]:
        server.serveonce()
    print "Python server has exited."
