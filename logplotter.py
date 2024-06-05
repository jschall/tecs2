from model import *
import pandas as pd
import matplotlib.pyplot as plt
from math import degrees

subs = {
    prop_J0T  : 1.04,
    prop_J0P  : 0.92,
    prop_dCTdJ: -0.09,
    prop_dCPdJ: -0.09,
    prop_maxCT: 0.07,
    prop_minCT: -0.015,
    prop_D: 0.6604,
    motor_kv: 207,
    motor_R: 0.0195,
    SKE: 0.5*Symbol("TAS")**2,
    mass:15,
    plane_wingspan: 6,
    plane_chord: 0.334,
    G: 9.81,
}

pitch = 0.5207
D = .6604
beta75 = atan(pitch/(2*pi*D*0.5*0.75))
phi75 = atan(J/(pi*0.75))
alpha75 = beta75 - phi75

df = pd.read_csv('data.csv')

prop_specific_power_out_truth = Symbol("STEdot")+flight_specific_power_required
prop_thrust_truth = prop_specific_power_out_truth/TAS * mass
prop_CT_truth = prop_thrust_truth/(rho*prop_n**2*prop_D**4)

STEdot_predicted = lambdify([Symbol("t"),rho,Symbol("TAS"),prop_omega,Symbol("STEdot"),bank_rad,plane_CD0],simplify(STEdot.subs(subs).evalf()))

prop_CT_truth = lambdify([Symbol("t"),rho,Symbol("TAS"),prop_omega,Symbol("STEdot"),bank_rad,plane_CD0],simplify(prop_CT_truth.subs(subs).evalf()))
J = lambdify([Symbol("t"),rho,Symbol("TAS"),prop_omega,Symbol("STEdot"),bank_rad],simplify(J.subs(subs).evalf()))
data = df.to_numpy()
alpha75 = lambdify([Symbol("t"),rho,Symbol("TAS"),prop_omega,Symbol("STEdot"),bank_rad],simplify(alpha75.subs(subs).evalf()))

fig, axes = plt.subplots(1, 1, sharex='all', sharey='all')

t_list = []
J_list = []
TAS_list = []
prop_CT_truth1 = []
prop_CT_truth2 = []
prop_CT_truth3 = []
prop_omega_list = []
STEdot_list = []
alpha75_list = []
STEdot_predicted_list = []
for row in data:
    t,rho,TAS,prop_omega,STEdot,bank_rad = tuple(row)
    prop_omega = abs(prop_omega)
    if prop_omega >= 10 and TAS > 12:
        t_list.append(t)
        J_list.append(J(*row))
        prop_CT_truth1.append(prop_CT_truth(*row,0.014))
        prop_omega_list.append(prop_omega)
        TAS_list.append(TAS)
        alpha75_list.append(degrees(alpha75(*row)))
        STEdot_list.append(STEdot)
        STEdot_predicted_list.append(STEdot_predicted(*row,0.014))

axes.set_title("Measured CT assuming CD0 = 0.015\nColor=TAS")
#axes.scatter(J_list,prop_CT_truth1,c=TAS_list,alpha=0.1)

plt.plot(t_list, STEdot_list)
plt.plot(t_list,STEdot_predicted_list)

fig.tight_layout()
plt.show()

#rho,TAS,prop_omega,STEdot,bank_rad = tuple(df.to_numpy())

#def testfunc(*args):
    #print(len(args))
    #print(type(args[0]))

#prop_CT_truth(*df.to_numpy().transpose())
