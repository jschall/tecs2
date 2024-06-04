from model import *
import pandas as pd
import matplotlib.pyplot as plt


subs = {
    prop_J0T  : 1,
    prop_J0P  : 0.92,
    prop_dCTdJ: -0.09,
    prop_dCPdJ: -0.09,
    prop_D: 0.6604,
    motor_kv: 207,
    motor_R: 0.0195,
    SKE: 0.5*Symbol("TAS")**2,
    mass:15,
    plane_wingspan: 6,
    plane_chord: 0.334,
    G: 9.81,
}

df = pd.read_csv('data.csv')

prop_specific_power_out_truth = Symbol("STEdot")+flight_specific_power_required
prop_thrust_truth = prop_specific_power_out_truth/TAS * mass
prop_CT_truth = prop_thrust_truth/(rho*prop_n**2*prop_D**4)

prop_CT_truth = lambdify([rho,Symbol("TAS"),prop_omega,Symbol("STEdot"),bank_rad,plane_CD0],simplify(prop_CT_truth.subs(subs).evalf()))
J = lambdify([rho,Symbol("TAS"),prop_omega,Symbol("STEdot"),bank_rad],simplify(J.subs(subs).evalf()))
data = df.to_numpy()

fig, axes = plt.subplots(1, 3, sharex='all', sharey='all')

J_list = []
TAS_list = []
prop_CT_truth1 = []
prop_CT_truth2 = []
prop_CT_truth3 = []
prop_omega_list = []
for row in data:
    rho,TAS,prop_omega,STEdot,bank_rad = tuple(row)
    if prop_omega > 10:
        J_list.append(J(*row))
        prop_CT_truth1.append(prop_CT_truth(*row,0.014))
        prop_CT_truth2.append(prop_CT_truth(*row,0.008))
        prop_CT_truth3.append(prop_CT_truth(*row,0.003))
        prop_omega_list.append(prop_omega)
        TAS_list.append(TAS)


axes[0].set_title("Predicted CT assuming CD0 = 0.014\nColor=TAS")
axes[0].scatter(J_list,prop_CT_truth1,c=TAS_list,alpha=0.2)
axes[1].set_title("Predicted CT assuming CD0 = 0.008\nColor=TAS")
axes[1].scatter(J_list,prop_CT_truth2,c=TAS_list,alpha=0.2)
axes[2].set_title("Predicted CT assuming CD0 = 0.003\nColor=TAS")
axes[2].scatter(J_list,prop_CT_truth3,c=TAS_list,alpha=0.2)


plt.show()

#rho,TAS,prop_omega,STEdot,bank_rad = tuple(df.to_numpy())

#def testfunc(*args):
    #print(len(args))
    #print(type(args[0]))

#prop_CT_truth(*df.to_numpy().transpose())
