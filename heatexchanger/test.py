from layer import Layer

m = Layer(2, 3)
m.cost = 1/m.Q
sol = m.localsolve(verbosity=1)
print sol("Q")
print sol(m.D_cold)
print sol(m.D_hot)
