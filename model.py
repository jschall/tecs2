from sympy import *
from math import radians
# states:
SKE = Symbol("SKE", real=True)
SPEdot = Symbol("SPEdot", real=True)
prop_omega = Symbol("prop_omega", real=True)
Vmotor = Symbol("Vmotor", real=True)
x = Matrix([SKE,SPEdot,prop_omega])

# inputs:
thetadot = Symbol("thetadot", real=True)
Vmotordot = Symbol("Vmotordot", real=True)
u = Matrix([thetadot, Vmotor])

# parameters:
motor_kv = Symbol("motor_kv", real=True, positive=True)
motor_R = Symbol("motor_R", real=True, positive=True)
prop_CT = Symbol("prop_CT", real=True, positive=True)
prop_D = Symbol("prop_D", real=True, positive=True)
prop_Ixx = Symbol("prop_Ixx", real=True, positive=True)
prop_efficiency = Symbol("prop_eta", real=True, positive=True)
prop_pitch = Symbol("prop_pitch", real=True, positive=True)

prop_J0T   = Symbol("prop_J0T", real=True, positive=True) # advance ratio at which thrust coefficient is zero
prop_J0P   = Symbol("prop_J0P", real=True, positive=True) # advance ratio at which power coefficient is zero
prop_dCTdJ = Symbol("prop_dCTdJ", real=True) # slope of thrust coefficient WRT advance ratio
prop_dCPdJ = Symbol("prop_dCPdJ", real=True) # slope of power coefficient WRT advance ratio

rho = Symbol("rho", real=True, positive=True)
mass = Symbol("mass", real=True, positive=True)
G = Symbol("G", real=True, positive=True)
plane_CD0 = Symbol("plane_CD0", real=True, positive=True)
plane_AR = Symbol("plane_AR", real=True, positive=True)
plane_S = Symbol("plane_S", real=True, positive=True)
bank_rad = Symbol("bank_rad", real=True, positive=True)
dt = Symbol("dt", real=True, positive=True)

TAS = sqrt(2*SKE)


####### Propeller and electric motor model: #######
# https://m-selig.ae.illinois.edu/props/propDB.html
# TODO design an equation to predict CT curve as a function of propeller measurements (diameter, pitch, chord, blades)

# Propeller speed in revolutions per second
prop_n = prop_omega/(2*pi)

# Propeller advance ratio
J = TAS/(prop_n*prop_D)

# Model CT as a function of J
prop_CT = (J-prop_J0T)*prop_dCTdJ # Thrust coefficient
prop_CP = (J-prop_J0P)*prop_dCPdJ # Power coefficient
prop_CQ = prop_CP/(2*pi) # Torque coefficient

# Propeller torque in N*m
prop_torque = -prop_CQ*rho*prop_n**2*prop_D**5

# Propeller thrust in N
prop_thrust = prop_CT * rho * prop_n**2 * prop_D**4

# Propeller specific power output in W/kg
prop_specific_power_out = TAS * prop_thrust / mass

# Motor torque in N*m
motor_flux_linkage = 60/(3**0.5*motor_kv*pi)

motor_torque = motor_flux_linkage*(Vmotor - motor_flux_linkage*prop_omega)/(motor_R*1.5)

#motor_torque = 30*(pi*Vmotor*motor_kv - 60*prop_omega)/(pi**2*motor_R*motor_kv**2)

# Propeller angular acceleration in rad/s/s
prop_omega_dot = (motor_torque+prop_torque) / prop_Ixx

####### Power curve model: #######
# https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node97.html

flight_specific_power_required = (0.5 * rho * TAS**3 * plane_S * plane_CD0 + (mass*G/cos(bank_rad))**2/(0.5*rho*TAS*plane_S * (pi*E*plane_AR)))/mass

# Rate of change of specific total energy in J/kg
STEdot = prop_specific_power_out - flight_specific_power_required

# Rate of change of specific kinetic energy
SKE_dot = STEdot - SPEdot

# Vertical acceleration
SPEdot_dot = thetadot * TAS * G

soln = solve([motor_torque+prop_torque], [prop_omega])

soln = solve(Eq(STEdot.subs({prop_omega:soln[1][0]}),Symbol("STEdot")), Vmotor)
print(srepr(soln))

subs = {
    prop_J0T  : 0.9,
    prop_J0P  : 0.92,
    prop_dCTdJ: -0.11,
    prop_dCPdJ: -0.09,
    prop_D: 0.6604,
    motor_kv: 220,
    motor_R: 0.0195,
    SKE: 0.5*Symbol("TECS.sp")**2,
    #SKE:0.5*17**2,
    mass:15,
    plane_S: 2.005,
    plane_AR: 16.6,
    bank_rad:radians(1)*Symbol("ATT.Roll"),
    #bank_rad:radians(0),
    plane_CD0: .006,
    rho:(1.0/Symbol("CTUN.E2T"))**0.5*1.225,
    #rho:1.225,
    G: 9.81,
    Vmotor: Symbol("VESC.VI")*Symbol("VESC.Duty"),
    #prop_omega:Symbol("ESC[4].RPM")/60*2*pi
}

eta = prop_CT*J/prop_CP

##pprint(soln)
#print(str((soln[1][0]/(2*pi)*60).evalf()).replace(" ",""))

subs[prop_omega] = soln[1][0]

#print(str(STEdot.subs(subs).evalf()).replace(" ",""))
#print(str(eta.subs(subs).evalf()).replace(" ",""))

#subs.update({prop_omega: soln[1][1], Vmotor:soln[1][0]})

#print(J.subs(subs).evalf(), eta.subs(subs).evalf())


#f = Matrix([SKE+dt*SKE_dot, SPEdot+dt*SPEdot_dot, prop_omega+dt*prop_omega_dot])

#p,q,n = (len(u), len(x), len(x))

#A = f.jacobian(x)
#B = f.jacobian(u)
