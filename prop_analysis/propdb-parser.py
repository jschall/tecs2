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
    "pitch":[],
    "D":[],
    "c75":[],
    "CT0":[],
    "dCTdJ":[],
    "etamax":[]
    }

index = []

for fn in fns:
    with open(fn) as f:
        if "volume-2" in fn:
            continue

        #if "apce" not in fn:
            #continue

        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        rows = [_ for _ in reader]

        try:
            rows = rows[rows.index(['J','CT','CP','eta'])+1:]
        except:
            continue


        try:
            name,diameter_inches,pitch_inches,rpm = re.match(r"^([^_]+)_([0-9.]+)x([0-9.]+).+([0-9]{4})\.txt", os.path.split(fn)[-1]).groups()
            geomfn = os.path.split(fn)[0]+"/"+name+"_"+diameter_inches+"x"+pitch_inches+"_geom.txt"
            D = float(diameter_inches)*0.0254
            pitch = float(pitch_inches)*0.0254
            rpm = float(rpm)
            if pitch > D*2:
                pitch /= 10
            if pitch < D*0.2:
                D /= 10

            c75 = None
            beta75 = None
            with open(geomfn) as geom_f:
                geom_reader = csv.reader(geom_f, delimiter=' ', skipinitialspace=True)
                geom_rows = [_ for _ in geom_reader][1:]

                for row in geom_rows:
                    if row[0] == '0.75':
                        c75 = float(row[1])
                        beta75 = float(row[2])
            if c75 is None:
                continue
        except:
            #print(fn, "can't be matched")
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
            CT.append(float(row[2]))
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

        #plt.plot(D,pitch,marker='x')

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
        propdata["pitch"].append(pitch)
        propdata["D"].append(D)
        propdata["c75"].append(c75)
        propdata["CT0"].append(-params[1]/params[0]*D/pitch)
        propdata["dCTdJ"].append(params[0])
        propdata["etamax"].append(max(eta))
        index.append(fn)

        if pitch > 0.3:
            print(fn)

df = pd.DataFrame(propdata, index=index)

X = df[['pitch','D','c75']]
y = df['CT0']

#fig, axes = plt.subplots(nrows=2, ncols=3)

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')

ax.scatter3D(df['pitch'],df['c75'],df['dCTdJ'])
plt.xlabel('pitch')
plt.ylabel('c75')

#df.plot.scatter(x='pitch', y='dCTdJ', ax=axes[0,0])
#df.plot.scatter(x='D', y='dCTdJ', ax=axes[0,1])
#df.plot.scatter(x='c75', y='dCTdJ', ax=axes[0,2])

#df.plot.scatter(x='pitch', y='CT0', ax=axes[1,0])
#df.plot.scatter(x='D', y='CT0', ax=axes[1,1])
#df.plot.scatter(x='c75', y='CT0', ax=axes[1,2])

from sklearn import linear_model
regr = linear_model.HuberRegressor()
regr.fit(X, y)
print("score", regr.score(X,y))
#print(regr.coef_)
#print(regr.intercept_)
print("CT0", regr.predict([[.53,.6604,.0067]])*.53/.6604)

#print(regr.predict([[.254,.2794,.125]])) # apce11x10

y = df['dCTdJ']
regr = linear_model.HuberRegressor()
regr.fit(X, y)
print("score", regr.score(X,y))
print("dCTdJ", regr.predict([[.53,.6604,.0067]]))
#print(regr.predict([[.254,.2794,.125]]))# apce11x10


#plt.title("APCE only")
#plt.subplot(211)
#plt.ylabel("CT")
#plt.xlabel("J WRT pitch speed")
#plt.subplot(212)
#plt.ylabel("CT")
#plt.xlabel("J")
plt.show()
