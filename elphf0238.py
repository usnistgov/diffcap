#!/usr/bin/env python

## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # Author: Jonathan Guyer
 # E-mail: <guyer@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are 
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 # 
 # #############################################################################
 ##

"""Divalent cations with impurity concentrations of $10^{-18}$
"""

import os
import sys
import json
import fipy as fp

jsonfile = sys.argv[1]

with open(jsonfile, 'rb') as ff:
    params = json.load(ff)

if "sumatra_label" in params:
    os.chdir(os.path.join("Data", params["sumatra_label"]))

L = fp.Variable("3.2 nm")
nx = int(params['ncells'])
dx = L / nx
print "L =", L
mesh = fp.Grid1D(nx=nx, dx=dx/L)

Vs = fp.Variable("0.018 l/mol")
R = fp.Variable("1 Nav * kB")
T = fp.Variable("298 K")
Faraday = fp.Variable("1 Nav * e")


xi = fp.CellVariable(mesh=mesh, name=r"\xi", value=0., hasOld=True)

def p(xi):
    return xi**3 * (6 * xi**2 - 15 * xi + 10)

def g(xi):
    return (xi * (1 - xi))**2

def pPrime(xi):
    return 30. * g(xi)

def gPrime(xi):
    return 2. * xi * (1 - xi) * (1 - 2 * xi)
    
def pPrimePrime(xi):
    return 30. * gPrime(xi)
    
def gPrimePrime(xi):
    return 2. * (1 - 6 * xi + 6 * xi**2)



phi = fp.CellVariable(mesh=mesh, name=r"\phi", value=0., hasOld=True)

interstitials = [
    fp.CellVariable(mesh=mesh, name="$e^-$", value=0., hasOld=True)
]

substitutionals = [
    fp.CellVariable(mesh=mesh, name="$M^+$", value=0., hasOld=True),
    fp.CellVariable(mesh=mesh, name="$A^-$", value=0., hasOld=True)
]

N = fp.Variable(name="$N$")
m = fp.Variable(name="$m$")


interstitials[0].z   = -1
substitutionals[0].z = +2
substitutionals[1].z = -2
N.z                  = 0

interstitials[0].muSL   = -41.041
substitutionals[0].muSL = -2.9188
substitutionals[1].muSL =  37.429
N.muSL                  =  41.410

gamma = fp.Variable("0.2 J/m**2")
delta_xi = fp.Variable("3e-11 m")

barrier = ((3. * Vs * gamma / delta_xi)
           / (R * T)).inBaseUnits()
print "barrier", barrier
for j in substitutionals:
    j.W = barrier

N.W = barrier

for j in interstitials:
    j.W = 0. 

xi.kappa = (6 * delta_xi * gamma
            / (L**2 * R * T / Vs)).inBaseUnits()

from scipy import optimize

def equilibrium(Y, interstitials, substitutionals, solvent, YM):
    components = interstitials + substitutionals + [solvent]
    N = len(components)
    I = len(interstitials)
    
    YmL = 1.
    YmS = 1.
    for j, Yj in enumerate(interstitials):
        YmL += Y[j]
        YmS += Y[j + N]
        
    sumL = 1.
    sumS = 1.
    for j, Yj in enumerate(substitutionals + [solvent]):
        sumL -= Y[I + j]
        sumS -= Y[I + j + N]
    
    out = [Y[1] - YM]
    
    chargeL = 0.
    chargeS = 0.
    for j, Yj in enumerate(components):
        out.append(fp.numerix.exp(Yj.muSL + Yj.z * Y[2*N]) * Y[j + N] / YmS - Y[j] / YmL)
        chargeL += Yj.z * Y[j]
        chargeS += Yj.z * Y[j + N]
    
    out += [sumL, sumS, chargeL, chargeS]
    
    return out

components = interstitials + substitutionals + [N]

YM = (Vs * str(params["concentration"])).inBaseUnits()
Y0 = [3e-6, YM, YM, 1, 2, 1, 3e-6, 3e-6, 0]
Y0 = optimize.fsolve(equilibrium, Y0, xtol=1e-12, 
                     args=(interstitials, substitutionals, N, YM))

Galvani = fp.Variable(value=(fp.Variable(str(params["galvani_potential"])) * Faraday / (R*T)).inBaseUnits(), 
                      name=r'\Delta\phi^\circ')

m.YL = 0.
m.YS = 0.
for i, j in enumerate(components):
    j.YL = Y0[i]
    j.YS = Y0[i+len(components)]
    m.YL += j.YL
    m.YS += j.YS
    
    print j.name, j.YL, j.YS

phi.L = -Y0[2 * len(components)] - Galvani
phi.S = 0
print m.name, m.YL, m.YS

print "deltaV", Y0[2*len(components)]

for i in substitutionals:
    i.muSL = i.muSL - N.muSL
    i.z = i.z - N.z
    i.W = i.W - N.W

for i in components:
    i.muSL = i.muSL - i.z * Galvani

x = mesh.cellCenters[0]

for j in interstitials:
    j.mubarS = j.muSL + fp.numerix.log(j.YS / m.YS) + j.z * phi.S
    j.mubarL = fp.numerix.log(j.YL / m.YL) + j.z * phi.L
    j.mubarS.name = r"$\bar{\mu}^S_{%s}$" % j.name.strip("$")
    j.mubarL.name = r"$\bar{\mu}^L_{%s}$" % j.name.strip("$")
    
    print j.name, j.mubarS, j.mubarL
    
for j in substitutionals:
    j.mubarS = j.muSL + fp.numerix.log(j.YS / N.YS) + j.z * phi.S
    j.mubarL = fp.numerix.log(j.YL / N.YL) + j.z * phi.L
    j.mubarS.name = r"$\bar{\mu}^S_{%s}$" % j.name.strip("$")
    j.mubarL.name = r"$\bar{\mu}^L_{%s}$" % j.name.strip("$")

    print j.name, j.mubarS, j.mubarL
    
for j in substitutionals + interstitials:
    j.constrain(j.mubarS, where=mesh.facesLeft)
    j.constrain(j.mubarL, where=mesh.facesRight)

N.mubarL = fp.numerix.log(N.YL / m.YL) + N.z * phi.L
N.mubarL.name = r"$\bar{\mu}^L_{%s}$" % N.name.strip("$")

if "restart" in params:
    values = fp.numerix.loadtxt(params["restart"], skiprows=1, unpack=True)
    
    xi.setValue(values[1])
    phi.setValue(values[2])
    
    # skip over stored electrochemical potentials
    _start = len(interstitials + substitutionals) + 3
    concentration = values[_start:_start+len(components)] 
    
    Ym = 1.
    for Y in concentration[0:len(interstitials)]:
        Ym += Y
        
    # force initial electrochemical potentials to be consistent with initial concentrations
    for j, Y in zip(interstitials, concentration[0:len(interstitials)]):
        j.setValue(j.muSL * p(xi) + fp.numerix.log(Y / Ym) + j.z * phi)
      
    YN = concentration[-1]
    _start = len(interstitials)
    for j, Y in zip(substitutionals, concentration[_start:_start + len(substitutionals)]):
        j.setValue(j.muSL * p(xi) + fp.numerix.log(Y / YN) + j.z * phi)
else:
    xi.setValue(1)
    xi.setValue(0, where=x > 0.25)

    phi.setValue(phi.L, where=xi == 0)
    phi.setValue(phi.S, where=xi == 1)
                
    for j in substitutionals + interstitials:
        j.setValue(j.mubarS, where=xi == 1)
        j.setValue(j.mubarL, where=xi == 0)

for i in components:
    i.Q = i.muSL * p(xi) + i.W * g(xi) + i.z * phi
    i.dmudxi = i.muSL * pPrime(xi) + i.W * gPrime(xi)
    i.d2mudxi2 = i.muSL * pPrimePrime(xi) + i.W * gPrimePrime(xi)
    i.Z = fp.numerix.exp(i - i.Q)
    i.Z_SG = fp.ScharfetterGummelFaceVariable(i - i.Q)

N.Y = fp.CellVariable(mesh=mesh, value=1.)
N.Y_SG = fp.FaceVariable(mesh=mesh, value=1.)
for i in substitutionals:
    N.Y += i.Z
    N.Y_SG += i.Z_SG
    
N.Y = 1. / N.Y
N.Y_SG = 1. / N.Y_SG
N.Y.name = N.name
N.Y_SG.name = N.name

m.Y = fp.CellVariable(mesh=mesh, value=1.)
m.Y_SG = fp.FaceVariable(mesh=mesh, value=1.)
for i in interstitials:
    m.Y -= i.Z
    m.Y_SG -= i.Z_SG
    
m.Y = 1. / m.Y
m.Y_SG = 1. / m.Y_SG
m.Y.name = 'Y_m'
m.Y_SG.name = 'Y_m'

interstitialCharge = 0.
interstitialEnthalpy = 0.
for i in interstitials:
    i.Y = m.Y * i.Z
    i.Y.name = i.name
    i.Y_SG = m.Y_SG * i.Z_SG
    i.Y_SG.name = i.name
    i.dmudY = (m.Y - i.Y) / (m.Y * i.Y)
    i.dmudY_SG = (m.Y_SG - i.Y_SG) / (m.Y_SG * i.Y_SG)
    i.dmudphi = i.z
    interstitialCharge += i.Y * i.z
    interstitialEnthalpy += i.Y * i.dmudxi

substitutionalCharge = 0.
substitutionalEnthalpy = 0.
for i in substitutionals:
    i.Y = N.Y * i.Z
    i.Y.name = i.name
    i.Y_SG = N.Y_SG * i.Z_SG
    i.Y_SG.name = i.name
    i.dmudY = (N.Y + i.Y) / (N.Y * i.Y)
    i.dmudY_SG = (N.Y_SG + i.Y_SG) / (N.Y_SG * i.Y_SG)
    i.dmudphi = i.z
    substitutionalCharge += i.Y * i.z
    substitutionalEnthalpy += i.Y * i.dmudxi
    
for i in interstitials:
    i.dYdphi = -i.Y * (i.z + interstitialCharge)
    i.dYdxi = -i.Y * (i.dmudxi + interstitialEnthalpy)
    
for i in substitutionals:
    i.dYdphi = -i.Y * (i.z - substitutionalCharge)
    i.dYdxi = -i.Y * (i.dmudxi - substitutionalEnthalpy)

enthalpy = N.dmudxi
enthalpyPrime = N.d2mudxi2
for i in interstitials + substitutionals:
    enthalpy += i.Y * i.dmudxi
    enthalpyPrime += i.dYdxi * i.dmudxi + i.Y * i.d2mudxi2
    
N.mubar = N.muSL * p(xi) + N.W * g(xi) + fp.numerix.log(N.Y / m.Y) + N.z * phi
    
D0 = fp.Variable("1e-5 cm**2/s")
tau = (L**2 / D0).inBaseUnits()
xi.M = (fp.Variable("1e4 cm**3 / (J*s)") / (Vs / (tau * R * T))).inBaseUnits()


for i in interstitials + substitutionals:
    i.D = 1.

for i in interstitials:
    i.M = i.D * m.Y * i.Y / (m.Y - i.Y)
    i.M_SG = i.D * m.Y_SG * i.Y_SG / (m.Y_SG - i.Y_SG)
    i.dMdmu = i.M * (m.Y / (m.Y - i.Y) + i.Y)
    i.dMdmu_SG = i.M_SG * (m.Y_SG / (m.Y_SG - i.Y_SG) + i.Y_SG)

for i in substitutionals:
    i.M = i.D * N.Y * i.Y / (N.Y + i.Y)
    i.M_SG = i.D * N.Y_SG * i.Y_SG / (N.Y_SG + i.Y_SG)
    i.dMdmu = i.M * (N.Y / (N.Y + i.Y) - i.Y)
    i.dMdmu_SG = i.M_SG * (N.Y_SG / (N.Y_SG + i.Y_SG) - i.Y_SG)
    
dt = fp.Variable(value=(fp.Variable("7.5e-10 s") / tau).inBaseUnits())
print "dt", dt #, (dt * L**2 / dx**2).inBaseUnits()

permittivity0 = (fp.Variable("1 eps0") 
                 / (L * L * Faraday * Faraday / (Vs * R * T))).inBaseUnits()

def dielectric(xi):
    return 78.49

def dielectricPrime(xi):
    return 0.
    
def dielectricPrimePrime(xi):
    return 0.

transient = fp.Variable(value=0.)

xi.eq = (fp.TransientTerm(coeff=transient / xi.M) 
         == fp.DiffusionTerm(coeff=xi.kappa)
         - enthalpy 
         + (permittivity0 * dielectricPrime(xi) / 2.) * xi.grad.dot(xi.grad))
         
xi.J  = (fp.TransientTerm(coeff=transient / xi.M) 
         == fp.DiffusionTerm(coeff=xi.kappa)
         - fp.ImplicitSourceTerm(coeff=enthalpyPrime 
                                 - (permittivity0 * dielectricPrimePrime(xi) / 2.) 
                                 * xi.grad.dot(xi.grad)))

## invYm = 1. / m.Y
invYmGrad = 0.
for i in interstitials:
    invYmGrad -= i.Z_SG * (i - i.Q).faceGrad
    
## invYN = 1. / N.Y
invYNGrad = 0.
for i in substitutionals:
    invYNGrad += i.Z_SG * (i - i.Q).faceGrad
                              
for j in interstitials:
    j.eq = (fp.TransientTerm(coeff=transient) 
            == (fp.DiffusionTerm(coeff=j.D)
                - fp.PowerLawConvectionTerm(coeff=j.M_SG * j.dmudY.faceGrad)
                + fp.ImplicitSourceTerm(coeff=(j.M_SG 
                                               * j.dmudY.faceGrad).divergence)
                + transient * j.dmudxi * (xi - xi.old) / dt 
                + transient * j.dmudphi * (phi - phi.old) / dt))
            
    dmudYdMdmugradmu = j.dmudY_SG * j.dMdmu_SG * j.faceGrad
    j.J = (fp.TransientTerm(coeff=transient) 
           == (fp.DiffusionTerm(coeff=j.D)
          - fp.PowerLawConvectionTerm(coeff=j.M_SG * j.dmudY.faceGrad)
               + fp.ImplicitSourceTerm(coeff=(j.M_SG 
                                              * j.dmudY.faceGrad).divergence)
               + fp.PowerLawConvectionTerm(coeff=dmudYdMdmugradmu)
               - fp.ImplicitSourceTerm(coeff=dmudYdMdmugradmu.divergence)
               + fp.ImplicitSourceTerm(coeff=j.dmudY 
                                       * (j.dMdmu_SG * j.faceGrad).divergence)))
           
    for i in [i for i in interstitials if i is not j]:
        j.eq += (((i.D * i.Y_SG / (m.Y_SG - i.Y_SG)) 
                  * i.faceGrad).divergence
                 - (i.arithmeticFaceValue * i.M_SG * invYmGrad).divergence
                 + i * (i.M_SG * invYmGrad).divergence)
 
        dmudYdMdmugradmu_cross = (i.D * i.Y_SG * j.Y_SG / (m.Y_SG - i.Y_SG)) * i.faceGrad
        j.J += (fp.PowerLawConvectionTerm(coeff=dmudYdMdmugradmu_cross)
                - fp.ImplicitSourceTerm(coeff=dmudYdMdmugradmu_cross.divergence)
                + fp.ImplicitSourceTerm(coeff=(i.M_SG * j.Y_SG * i.faceGrad).divergence / m.Y))
                
                              
for j in substitutionals:
    j.eq = (fp.TransientTerm(coeff=transient) 
            == (fp.DiffusionTerm(coeff=j.D)
                  - fp.PowerLawConvectionTerm(coeff=j.M_SG * j.dmudY.faceGrad)
                  + fp.ImplicitSourceTerm(coeff=(j.M_SG * j.dmudY.faceGrad).divergence)
                  + transient * j.dmudxi * (xi - xi.old) / dt 
                  + transient * j.dmudphi * (phi - phi.old) / dt))
            
    dmudYdMdmugradmu = j.dmudY_SG * j.dMdmu_SG * j.faceGrad
    j.J = (fp.TransientTerm(coeff=transient) 
           ==  (fp.DiffusionTerm(coeff=j.D)
                - fp.PowerLawConvectionTerm(coeff=j.M_SG * j.dmudY.faceGrad)
                + fp.ImplicitSourceTerm(coeff=(j.M_SG * j.dmudY.faceGrad).divergence)
                + fp.PowerLawConvectionTerm(coeff=dmudYdMdmugradmu)
                - fp.ImplicitSourceTerm(coeff=dmudYdMdmugradmu.divergence)
                + fp.ImplicitSourceTerm(coeff=j.dmudY * (j.dMdmu_SG * j.faceGrad).divergence)))

    for i in [i for i in substitutionals if i is not j]:
        j.eq -= (((i.D * i.Y_SG / (N.Y_SG + i.Y_SG)) * i.faceGrad).divergence
                 - (i.arithmeticFaceValue * i.M_SG * invYNGrad).divergence
                 + i * (i.M_SG * invYNGrad).divergence)
 
        dmudYdMdmugradmu_cross = (i.D * i.Y_SG * j.Y_SG / (N.Y_SG + i.Y_SG)) * i.faceGrad
        j.J += (fp.PowerLawConvectionTerm(coeff=dmudYdMdmugradmu_cross)
                - fp.ImplicitSourceTerm(coeff=dmudYdMdmugradmu_cross.divergence)
                + fp.ImplicitSourceTerm(coeff=(i.M_SG * j.Y_SG * i.faceGrad).divergence / N.Y))
                 
charge = fp.CellVariable(mesh=mesh, value=N.z)
chargePrime = 0.
for i in interstitials + substitutionals:
    charge += i.z * i.Y
    chargePrime += i.z * i.dYdphi
    
charge.name = r"\rho"
    
phi.eq = fp.DiffusionTerm(coeff=permittivity0 * dielectric(xi)) == -charge

phi.J = fp.DiffusionTerm(coeff=permittivity0 * dielectric(xi)) == -fp.ImplicitSourceTerm(coeff=chargePrime)

phi.constrain(0., where=mesh.facesLeft)

relaxation = fp.Variable(value=params["relaxation"])
for var in [xi, phi] + interstitials + substitutionals:
    var.delta = fp.CellVariable(mesh=mesh, name="d" + var.name, value=0., hasOld=True)
    var.F = fp.CellVariable(mesh=mesh)
    var.J = var.J + relaxation * var.F / mesh.cellVolumes

for var in [xi] + interstitials + substitutionals:
    var.delta.constrain(0., where=mesh.facesLeft)
    var.delta.constrain(0., where=mesh.facesRight)
    
phi.delta.constrain(0., where=mesh.facesLeft)

surfaceEnergy =  (xi.kappa * xi.grad.dot(xi.grad)
                  - dielectric(xi) * permittivity0 
                  * phi.grad.dot(phi.grad)).cellVolumeAverage

surfaceCharge =  (p(xi) * charge).cellVolumeAverage

N.first = N.mubar - (xi.kappa * xi.grad.dot(xi.grad)
                     - dielectric(xi) * permittivity0 
                     * phi.grad.dot(phi.grad)) / 2
N.first.name = "$1^{st}$"

from matplotlibElPhFViewer import MatplotlibElPhFViewer
#viewer = MatplotlibElPhFViewer(phase=xi,
#                               potential=phi,
#                               components=[i.Y for i in interstitials + substitutionals 
#                                           + [N]],
#                               charge=charge,
#                               potentials=[i - i.mubarS for i in interstitials + substitutionals],
#                               limits={
#                                   'phasemax':1.,
#                                   'phasemin':0.,
#                               })
                               
def dampBrownLindsay(deltaV):
    # Voltage damping of Brown & Lindsay, 
    # Solid-State Electronics 19 (1976) 991
    
    dVsgn = fp.numerix.sign(deltaV)
    dVabs = fp.numerix.absolute(deltaV)
    deltaV = fp.numerix.where((1. < dVabs) & (dVabs < 3.7), 
                              dVsgn * dVabs**0.2, deltaV)
    deltaV = fp.numerix.where(dVabs >= 3.7, dVsgn * fp.numerix.log(dVabs), deltaV)
    
    return deltaV


maxmu = 1e6

outer = 0
while outer < int(params["outer_sweeps"]) and maxmu > 1e-30:
    for var in [xi, phi] + interstitials + substitutionals:
        var.updateOld()
        var.delta.setValue(0.)
        var.delta.updateOld()
        
    ress = []

    xi.F.setValue(xi.eq.justResidualVector(var=xi,
                                           dt=dt))
    res = xi.J.sweep(var=xi.delta, dt=dt)
    ress += [res]
    
    xi.setValue(xi + xi.delta())

    for i in [phi] + interstitials + substitutionals:
        for sweep in range(10):
            i.delta.setValue(0.)
            i.F.setValue(i.eq.justResidualVector(var=i, 
                                                 dt=dt))
            res = i.J.sweep(var=i.delta, dt=dt)
            
            i.setValue(i + dampBrownLindsay(i.delta()))
            
#             print res

        ress += [res]


    outer += 1
    
    maxmu = 0
    for var in interstitials + substitutionals:
        maxmu = max(maxmu, abs(var - var.mubarS).max().value)

    #viewer.plot()
    if outer == int(params["outer_sweeps"]) -1: 
	with open("test.txt", "a") as ff:
	    ff.write("{0} {1} {2}\n".format(Galvani(), surfaceEnergy(), surfaceCharge())) 
    else: 
	print Galvani(), surfaceEnergy(), surfaceCharge()

fp.TSVViewer(vars=[xi, phi] + interstitials + substitutionals + [j.Y for j in components]).plot(filename="output.tsv")
