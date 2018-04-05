from SimpleWebSocketServer import SimpleWebSocketServer, WebSocket
import json
from layer import Layer
from gencsm import gencsm


def gensoltxt(m, sol):
    with open("sol.txt", "w") as f:
        for var in sorted(m.varkeys, key=str):
            f.write("%s [%s]\t\t%f\n" % (var, var.unitstr(dimless="-"),
                                         sol["variables"][var]))


class HXGPServer(WebSocket):

    def handleMessage(self):
        print "< received", repr(self.data)
        try:
            self.data = json.loads(self.data)
            print self.data

            try:
                Nairpipes = self.data["Air_Channels"]
                Nwaterpipes = self.data["Water_Channels"]
            except Exception:
                Nairpipes, Nwaterpipes = 5, 5

            for key in gpkit.MODELNUM_LOOKUP:
                gpkit.MODELNUM_LOOKUP[key] = 0
            m = Layer(Nairpipes, Nwaterpipes)
            m.cost = (m.D_air+m.D_wat)/m.Q

            for name, value in self.data.items():
                try:
                    key = m.design_parameters[name]
                    m.substitutions[key] = value
                except KeyError as e:
                    print repr(e)

            sol = m.localsolve()
            gensoltxt(m, sol)
            gencsm(m, sol)

            self.send({"status": "optimal",
                       "msg": "Successfully optimized!"})
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


if __name__ == "__main__":
    # TODO: uncomment to produce the initial CSM file before serving
    m = Layer(5, 5)
    m.cost = (m.D_air+m.D_wat)/m.Q
    sol = m.localsolve()
    gensoltxt(m, sol)
    gencsm(m, sol)
    # SimpleWebSocketServer('', 8000, HXGPServer).serveforever()
