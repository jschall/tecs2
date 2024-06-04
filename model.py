from sympy import *
from math import radians

####### State vector (x) #######
# SKE is specific kinetic energy, a function of true airspeed.
# SPEdot is rate of change of specific potential energy, a function of climb rate.
# prop_omega is the rotation rate of the propeller in radians per second.
# Vmotor is the voltage in to the motor. This is a state in the controller because we want to apply a limit to the rate of change of throttle.
SKE = Symbol("SKE", real=True)
SPEdot = Symbol("SPEdot", real=True)
prop_omega = Symbol("prop_omega", real=True, positive=True)
Vmotor = Symbol("Vmotor", real=True, positive=True)
x = Matrix([SKE,SPEdot,prop_omega,Vmotor])

####### Input vector (u) #######
# thetadot is the rate of change of pitch angle. It can be related to vertical acceleration through the centripetal force equation, under the assumption that alpha is not varying.
# Vmotordot is the rate of change of Vmotor
thetadot = Symbol("thetadot", real=True)
Vmotordot = Symbol("Vmotordot", real=True)
u = Matrix([thetadot, Vmotordot])

####### Parameters #######
# Motor parameters
# kv: motor kv in V/RPM
# R: motor resistance
motor_kv = Symbol("motor_kv", real=True, positive=True)
motor_R = Symbol("motor_R", real=True, positive=True)

# Propeller parameters
# J0T: advance ratio at which thrust coefficient is zero. Typical value 1.25*pitch/diameter, but varies up to 30%.
# J0P: advance ratio at which power coefficient is zero
# dCTdJ: slope of thrust coefficient WRT advance ratio
# dCPdJ: slope of power coefficient WRT advance ratio
prop_D = Symbol("prop_D", real=True, positive=True)
prop_Ixx = Symbol("prop_Ixx", real=True, positive=True)
prop_J0T   = Symbol("prop_J0T", real=True, positive=True)
prop_J0P   = Symbol("prop_J0P", real=True, positive=True) #
prop_dCTdJ = Symbol("prop_dCTdJ", real=True, negative=True)
prop_dCPdJ = Symbol("prop_dCPdJ", real=True, negative=True)

# Plane parameters
# wingspan: wingspan in meters
# chord: mean chord length in meters
# CD0: drag coefficient of airplane at alpha=0 without propeller. Defined as drag force = 0.5*rho*V**2*S*CD0, where S is wing area. Typical range 0.01 for a very efficient glider to 0.05 for a draggy airframe.
# mass: mass of airplane in kg
# bank_rad: the current angle of bank in radians
plane_wingspan = Symbol("plane_wingspan", real=True, positive=True)
plane_chord = Symbol("plane_chord", real=True, positive=True)
plane_CD0 = Symbol("plane_CD0", real=True, positive=True)
mass = Symbol("mass", real=True, positive=True)
bank_rad = Symbol("bank_rad", real=True)

# rho: air density in kg/m**3
rho = Symbol("rho", real=True, positive=True)

# G: gravity in m/s/s
G = Symbol("G", real=True, positive=True)

# dt: delta time in seconds
dt = Symbol("dt", real=True, positive=True)

TAS = sqrt(2*SKE)


####### Propeller and electric motor model: #######
# https://m-selig.ae.illinois.edu/props/propDB.html

# Propeller speed in revolutions per second
prop_n = prop_omega/(2*pi)

# Propeller advance ratio
J = TAS/(prop_n*prop_D)

# Model CT as a function of J
prop_CT = (J-prop_J0T)*prop_dCTdJ # Thrust coefficient
prop_CP = (J-prop_J0P)*prop_dCPdJ # Power coefficient
prop_CQ = prop_CP/(2*pi) # Torque coefficient

# eta is efficiency
prop_eta = prop_CT*J/prop_CP

# Propeller torque in N*m
prop_torque = -prop_CQ*rho*prop_n**2*prop_D**5

# Propeller thrust in N
prop_thrust = prop_CT * rho * prop_n**2 * prop_D**4

# Propeller specific power output in W/kg
prop_specific_power_out = TAS * prop_thrust / mass

# Propeller specific shaft power in W/kg
prop_specific_power_in = prop_CP*rho*prop_n**3*prop_D**5 / mass

# Motor torque in N*m
motor_flux_linkage = 60/(3**Rational(1,2)*motor_kv*pi)

motor_torque = motor_flux_linkage*(Vmotor - motor_flux_linkage*prop_omega)/(motor_R*Rational(3,2))

# Propeller angular acceleration in rad/s/s
prop_omega_dot = (motor_torque+prop_torque) / prop_Ixx

####### Power curve model: #######
# https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node97.html

plane_S = plane_wingspan*plane_chord
plane_AR = plane_wingspan/plane_chord

flight_specific_power_required = (Rational(1,2) * rho * TAS**3 * plane_S * plane_CD0 + (mass*G/cos(bank_rad))**2/(Rational(1,2)*rho*TAS*plane_S * (pi*E*plane_AR)))/mass

# Rate of change of specific total energy in J/kg
STEdot = prop_specific_power_out - flight_specific_power_required

# Rate of change of specific kinetic energy
SKE_dot = STEdot - SPEdot

# Vertical acceleration
SPEdot_dot = thetadot * TAS * G

def solveforvmotor():
    soln = solve(motor_torque+prop_torque, prop_omega)
    #print(simplify(soln[1]))

    prop_omega_expr = simplify(soln[1])

    subx,expr = cse(prop_thrust.subs({prop_omega:prop_omega_expr})-Symbol("prop_thrust"),ignore=(Vmotor,))

    soln = solve(expr, Vmotor, simplify=False)
    for s in soln:
        print("\n\n")
        pprint(s[0].xreplace(dict(subx)).xreplace(dict(subx)).xreplace(dict(subx)).xreplace(dict(subx)))
        print("\n\n")

def generate_mavexplorer_strings():
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
        plane_CD0: .014,
        G: 9.81,
    }

    pitch = 0.5207
    D = .6604
    beta75 = atan(pitch/(2*pi*D*0.5*0.75))
    phi75 = atan(J/(pi*0.75))
    alpha75 = beta75 - phi75
    prop_specific_power_out_truth = Symbol("STEdot")+flight_specific_power_required
    prop_thrust_truth = prop_specific_power_out_truth/TAS * mass
    prop_CT_truth = prop_thrust_truth/(rho*prop_n**2*prop_D**4)

    stedot_truth_graphstring = "TEC3.KED+TEC3.PED"+"#stedot_truth"
    stedot_graphstring = str(simplify(STEdot.subs(subs).evalf())).replace(" ","")+"#stedot_predicted"
    prop_power_graphstring = str(simplify(prop_specific_power_out.subs(subs).evalf())).replace(" ","")+"#prop_power_out_predicted"
    prop_CT_truth_graphstring = str(simplify(prop_CT_truth.subs(subs).evalf())).replace(" ","")+"#prop_CT_truth"
    prop_power_in_graphstring = str(simplify(prop_specific_power_in.subs(subs).evalf())).replace(" ","")+"#prop_power_in_predicted"
    prop_power_in_truth_graphstring = "ESC[4].Curr*ESC[4].Volt/15#prop_power_in_truth"
    J_graphstring = str(simplify(J.subs(subs).evalf())).replace(" ","")+"#J"
    eta_graphstring = "min(max("+str(simplify(prop_eta.subs(subs).evalf())).replace(" ","")+",0),10)#eta_predicted"
    alpha75_graphstring = "degrees("+str(simplify(alpha75.subs(subs).evalf())).replace(" ","")+")#alpha75"

    #print("graph %s %s %s %s %s" % (prop_power_truth_graphstring, stedot_graphstring, prop_power_in_graphstring, prop_power_in_truth_graphstring, alpha75_graphstring))
    print(J_graphstring)
    print(prop_CT_truth_graphstring)







#generate_mavexplorer_strings()

#

###pprint(soln)
##print(str((soln[1][0]/(2*pi)*60).evalf()).replace(" ",""))

#subs[prop_omega] = soln[1][0]

#print(str(STEdot.subs(subs).evalf()).replace(" ",""))
#print(str(eta.subs(subs).evalf()).replace(" ",""))

#subs.update({prop_omega: soln[1][1], Vmotor:soln[1][0]})

#print(J.subs(subs).evalf(), eta.subs(subs).evalf())


#f = Matrix([SKE+dt*SKE_dot, SPEdot+dt*SPEdot_dot, prop_omega+dt*prop_omega_dot])

#p,q,n = (len(u), len(x), len(x))

#A = f.jacobian(x)
#B = f.jacobian(u)
