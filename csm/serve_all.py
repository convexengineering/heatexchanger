import os
import sys
import webbrowser
from subprocess import Popen, call
from time import sleep


def wait_until_file_exists(filename):
    while not os.path.exists(filename):
        sleep(1)
    print "found", filename


try:
    os.remove("HX.csm")
    os.remove("HX_000.egads")
except OSError:
    pass
Popen(["python", os.sep.join(["..", "heatexchanger", "server.py"])])
wait_until_file_exists("HX.csm")
Popen(["serveCSM", "HX.csm"])
wait_until_file_exists("HX_000.egads")
sleep(5)
webbrowser.open_new(os.sep.join(["..", "ESP", "ESP-localhost7681.html"]))
