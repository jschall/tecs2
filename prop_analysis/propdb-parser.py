import glob
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
import re
import os
from math import *

import pandas as pd

def line(x,m,b):
    return m*x+b

def two_lines(x,ptx,pty,b1,b2):
    if x < ptx:
        return pty+b1*(x-ptx)
    else:
        return pty+b2*(x-ptx)

def one_line(x,ptx,pty,b1,b2):
    return pty+b2*(x-ptx)

def zero_intercept(ptx,pty,b1,b2):
    line1_zero = ptx-pty/b1
    line2_zero = ptx-pty/b2

    if abs(two_lines(line1_zero,ptx,pty,b1,b2)) < abs(two_lines(line2_zero,ptx,pty,b1,b2)):
        return line1_zero
    return line2_zero

fns = glob.glob('UIUC-propDB/**/data/*.txt')

propdata = {
    "pitch/D":[],
    "pitch":[],
    "D":[],
    "c75":[],
    "c75/R":[],
    "CT0":[],
    "dCTdJ":[],
    "etamax":[],
    "RPM/D":[],
    "dalphadJ0*C75/D":[],
    "ctmin":[]
    }

maxrpms = {}

for fn in fns:
    try:
        name,diameter_inches,pitch_inches,rpm = re.match(r"^([^_]+)_([0-9.]+)x([0-9.]+).+([0-9]{4})\.txt", os.path.split(fn)[-1]).groups()
        rpm = float(rpm)

        maxrpm = maxrpms.get((name,diameter_inches,pitch_inches),0)
        if rpm > maxrpm:
            maxrpms[(name,diameter_inches,pitch_inches)] = rpm
    except:
        continue

index = []

for fn in fns:
    with open(fn) as f:
        if "volume-2" in fn:
            continue

        #if "apc" not in fn:
            #continue

        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        rows = [_ for _ in reader]

        try:
            rows = rows[rows.index(['J','CT','CP','eta'])+1:]
        except:
            continue


        try:
            name,diameter_inches,pitch_inches,rpm = re.match(r"^([^_]+)_([0-9.]+)x([0-9.]+).+([0-9]{4})\.txt", os.path.split(fn)[-1]).groups()
            print(diameter_inches)
        except:
            print(fn, "can't be matched")
            continue

        try:
            geomfn = os.path.split(fn)[0]+"/"+name+"_"+diameter_inches+"x"+pitch_inches+"_geom.txt"
            D = float(diameter_inches)*0.0254
            pitch = float(pitch_inches)*0.0254
            rpm = float(rpm)
            if pitch > D*2:
                pitch /= 10
            if pitch < D*0.2:
                D /= 10

            #if maxrpms[(name,diameter_inches,pitch_inches)] > rpm:
                #continue

            c75 = None
            beta75 = None
            with open(geomfn) as geom_f:
                geom_reader = csv.reader(geom_f, delimiter=' ', skipinitialspace=True)
                geom_rows = [_ for _ in geom_reader][1:]

                for row in geom_rows:
                    if row[0] == '0.75':
                        c75 = float(row[1])*D/2
                        beta75 = float(row[2])
            if c75 is None:
                continue
        except:
            continue

        J = []
        CT = []
        CP = []
        eta = []

        prev_J = None
        for row in rows:
            if float(row[0]) == prev_J:
                continue
            prev_J = float(row[0])
            J.append(float(row[0]))
            CT.append(float(row[1]))
            CP.append(float(row[2]))
            eta.append(float(row[3]))

        J = np.asarray(J)
        CT = np.asarray(CT)
        CP = np.asarray(CP)
        eta = np.asarray(eta)

        if len(J) < 12:
            continue

        def objective_CT(params):
            global J, CT
            m,b = params
            CT_linefit = np.asarray([line(j, *params) for j in J[-6:]])
            return sum((CT_linefit-CT[-6:])**2)

        x0 = [0,0]
        res = minimize(objective_CT, x0)
        params=res.x

        #print(fn, res.fun, params)
        #print(len(J), J[0], J[-1])
        CT_linefit = np.asarray([line(j, *params) for j in J])

        pitch_calculated = 2*pi*tan(radians(beta75))*D*0.5*0.75

        beta75_calculated = atan(pitch/(2*pi*D*0.5*0.75))

        dalphadJ0 = -1.33333333333333/(pi*(1 + 1.77777777777778*pitch**2/(pi**2*D**2)))

        #plt.plot(D,pitch,marker='x')

        propdata["ctmin"].append(min(CT))

        if min(CT) < -max(CT)*0.2:
            plt.plot(J, CT)

        J_corrected = J*D/pitch
        J_corrected2 = J*D/pitch_calculated

        #print(J[list(eta).index(max(list(eta)))])
        #plt.plot(J_corrected2,eta/J,color='red',alpha=0.2)

        #plt.subplot(211)
        #plt.plot(J_corrected,CT,linestyle='',marker='x')
        #plt.plot(J_corrected,CT_linefit,label=fn)
        #plt.subplot(212)
        #plt.plot(J,CT,linestyle='',marker='x')
        #plt.plot(J,CT_linefit,label=fn)
        #print(fn)
        propdata["dalphadJ0*C75/D"].append(dalphadJ0*c75/D)
        propdata["pitch/D"].append(pitch/D)
        propdata["pitch"].append(pitch)
        propdata["D"].append(D)
        propdata["c75"].append(c75)
        propdata["c75/R"].append(c75/D*2)
        propdata["CT0"].append(-params[1]/params[0])
        propdata["dCTdJ"].append(params[0])
        propdata["etamax"].append(max(eta))
        propdata["RPM/D"].append(rpm/D/120)
        index.append(fn)

        if pitch > 0.3:
            print(fn)

df = pd.DataFrame(propdata, index=index)

print(df[df.ctmin == df.ctmin.min()])
print(df.min())
print()
print(df.mean())
print()
print(df.std())
print()

X = df[['pitch/D','c75/R']]
y = df['CT0']


#fig = plt.figure(figsize=(10, 10))
#ax = plt.axes(projection='3d')
#ax.scatter3D(df['dalphadJ0*C75/D'],df['pitch/D'],df['CT0'])
#plt.xlabel('dalphadJ0*C75/D')
#plt.ylabel('D')

#fig, axes = plt.subplots(nrows=2, ncols=3)
#df.plot.scatter(x='pitch/D', y='etamax', ax=axes[0,0])
#df.plot.scatter(x='D', y='etamax', ax=axes[0,1])
#df.plot.scatter(x='c75/R', y='etamax', ax=axes[0,2])

#df.plot.scatter(x='pitch/D', y='CT0', ax=axes[1,0])
#df.plot.scatter(x='D', y='CT0', ax=axes[1,1])
#df.plot.scatter(x='c75/R', y='CT0', ax=axes[1,2])

from sklearn import linear_model
regr = linear_model.HuberRegressor()
regr.fit(X, y)
print("score", regr.score(X,y))

print("CT0", regr.predict([[0.7884615384615384,.0067]]))

#print(regr.predict([[.254,.2794,.125]])) # apce11x10

X = df[['dalphadJ0*C75/D']]
y = df['dCTdJ']
regr = linear_model.HuberRegressor()
regr.fit(X, y)
print(regr.coef_)
print(regr.intercept_)

D = .6604
pitch = .5207
c75 = 0.022
dalphadJ0 = -1.33333333333333/(pi*(1 + 1.77777777777778*pitch**2/(pi**2*D**2)))

print("dalphadJ0", dalphadJ0)
print("score", regr.score(X,y))
print("dCTdJ", regr.predict([[dalphadJ0*0.022/D]]))
#print(regr.coef_*np.matrix([.53,.6604,.0067]).transpose()+regr.intercept_)
#print(regr.coef_.transpose())
#print(regr.intercept_)

#print(regr.predict([[.254,.2794,.125]]))# apce11x10


#plt.title("APCE only")
#plt.subplot(211)
#plt.ylabel("CT")
#plt.xlabel("J WRT pitch speed")
#plt.subplot(212)
#plt.ylabel("CT")
#plt.xlabel("J")
plt.show()
