from sympy import *
from math import radians

####### State vector (x) #######
# SKE is specific kinetic energy, a function of true airspeed.
# SPEdot is rate of change of specific potential energy, a function of climb rate.
# prop_omega is the rotation rate of the propeller in radians per second.
# Vmotor is the voltage in to the motor. This is a state in the controller because we want to apply a limit to the rate of change of throttle.
SKE = Symbol("SKE", real=True)
SPEdot = Symbol("SPEdot", real=True)
prop_omega = Symbol("prop_omega", real=True)
Vmotor = Symbol("Vmotor", real=True)
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
prop_dCTdJ = Symbol("prop_dCTdJ", real=True)
prop_dCPdJ = Symbol("prop_dCPdJ", real=True)

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
bank_rad = Symbol("bank_rad", real=True, positive=True)

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

plane_S = plane_wingspan*plane_chord
plane_AR = plane_wingspan/plane_chord

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
    plane_wingspan: 6,
    plane_chord: 0.334,
    bank_rad:radians(1)*Symbol("ATT.Roll"),
    #bank_rad:radians(0),
    plane_CD0: .016,
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
