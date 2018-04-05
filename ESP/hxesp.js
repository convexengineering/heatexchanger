var url = "localhost:8000" //prompt("Enter hostname:port for the gp server", "localhost:8000")

window.gp = {
  optcount: 0,

  dom: {},

  sol: {},

  esp: {build_outdated: false,
        update: function() {
            gp.dom.buildButton.disabled = false
            browserToServer("getCsmFile|")
            window.oldactivateBuildButton()
          },
        },

  awaiting_response: false,
  last_sent: null,
  sendpmtrs: function() {
    gp.dom.optimizeButton.innerText = "Optimizing..."
    gp.dom.optimizeButton.style.backgroundColor = "#DA70D6"
    gp.optcount++

    for (var i=0; i < pmtr.length; i++) {
      console.log(i, pmtr[i].name, pmtr[i].value[0], gp.sol[pmtr[i].name])
      gp.sol[pmtr[i].name] = pmtr[i].value[0]
    }

    console.log(gp.sol)
    gp.websocket.send(JSON.stringify(gp.sol))
    // gp.websocket.send("sol")
  },

  websocket: new WebSocket("ws://"+url+"/")
}

// gp.websocket.onopen = gp.sendpmtrs
gp.websocket.onmessage = function(evt) {
  data = JSON.parse(evt.data);
  console.log("Data received:", data)
  postMessage("GP: " + data.msg)
  if (data.status == "optimal") {
    gp.dom.optimizeButton.innerText = "Optimized"
    gp.dom.optimizeButton.style.backgroundColor = null
    gp.esp.update()
  } else {
    gp.dom.optimizeButton.innerText = "Error"
    gp.dom.optimizeButton.style.backgroundColor = "#FF3F3F"
  }
}

window.oldactivateBuildButton = activateBuildButton
window.activateBuildButton = function() {
  gp.dom.optimizeButton.innerText = "Press to Optimize"
  gp.dom.optimizeButton.style.backgroundColor = "#3FFF3F";
}

gp.dom.buttonForm = document.getElementById("butnfrm")
gp.dom.buildButton = document.getElementById("buildButton")
gp.dom.buildButton.disabled = true
gp.dom.optimizeButton = document.createElement("button")
gp.dom.optimizeButton.id = "optButton"
gp.dom.optimizeButton.type = "button"
gp.dom.optimizeButton.innerText = "Optimized"
gp.dom.optimizeButton.onclick = gp.sendpmtrs
gp.dom.buttonForm.insertBefore(gp.dom.optimizeButton, gp.dom.buildButton)
