from SimpleWebSocketServer import SimpleWebSocketServer, WebSocket
import json
from layer import Layer
from time import sleep


class HXGPServer(WebSocket):

    def handleMessage(self):
        print "< received", repr(self.data)
        try:
            sleep(1)
            print self.data

            m = Layer(Nairpipes=self.data["Nairpipes"],
                      Nwaterpipes=self.data["Nwaterpipes"])

            conversion_dictionary = {
                "Rbar": m.water.rho
            }


            subs = {}
            for key, value in self.data.items():
                if key in conversion_dictionary:
                    subs[conversion_dictionary[key]] = value
                else:
                    print "Key not found conversion_dictionary:", key

            m.localsolve()

            self.send({"status": "optimal",
                       "msg": "Successfully optimized!"})
        except Exception, e:
            self.send({"status": "unknown", "msg": "The last solution"
                      " raised an exception; tweak it and send again."})
            print type(e), e

    def send(self, msg):
        self.sendMessage(unicode(json.dumps(msg)))

    def handleConnected(self):
        print self.address, "connected"

    def handleClose(self):
        print self.address, "closed"


SimpleWebSocketServer('', 8000, HXGPServer).serveforever()
