#!/bin/sh
cd heatexchanger
python server.py > server.log &  # note: deletes sol.txt and HX.egads
while [ ! -f sol.txt ]; do sleep 1; done  # poll every second
serveCSM HX.csm  > csm.log &
while [ ! -f HX.egads ]; do sleep 1; done  # poll every second
python -m webbrowser -t ../ESP/ESP-localhost7681.html

