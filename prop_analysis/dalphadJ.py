from sympy import *

pitch,D,J = symbols("pitch D J")

# from blade element theory, at the 75% radius:
beta75 = atan(pitch/(2*pi*D*0.5*0.75))
phi75 = atan(J/(pi*0.75))
alpha75 = beta75 - phi75
J_sub = solve(alpha75, J)[0]

dalphadJ = diff(alpha75,J).subs({J:J_sub})

print(dalphadJ)
