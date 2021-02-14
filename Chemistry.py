#!/usr/bin/env python
# coding: utf-8

# In[60]:


from chempy import balance_stoichiometry
from sympy import *

A_N = 6.0221409e+23
R = 0.08205

H = 1.008
Li = 6.92
Fe = 55.845
Ne = 20.180
Cl = 35.45
C = 12.011
O = 15.999
Co = 58.933
P = 30.974
S = 32.06
N = 14.007
F = 18.998
Mg = 24.305
I = 126.90
Ca = 40.078
Cu = 63.546
Pb = 207.2
Ti = 47.867
Zn = 65.38
Na = 22.990
B = 10.81
K = 39.098
Br = 79.904
Ag = 107.87
Cr = 51.9961
Ba = 137.33
Al = 26.981539
Hg = 200.59
Ar = 39.95
Xe = 131.29
Kr = 83.798
He = 4.0026
acetate = C+H*3+C+O*2


def emper_form(per1,M1,per2,M2,per3=0,M3=0):
    if per3 == 0:
        moles = [(per1/100)/M1,(per2/100)/M2]
    else:
        moles = [(per1/100)/M1,(per2/100)/M2,(per3/100)/M3,]
    return [n/min(moles) for n in moles]

def atom2mole(a):
    return a/A_N

def mole2atom(n):
    return n*A_N

def gram2mole(m, M):
    return m/M

def mole2gram(n, M):
    return n*M

def grams_atoms(m, M):
    return (m/M)*A_N

def atoms_grams(a, M):
    return (a/A_N)*M

def molarity(M, N, V):
    m, n, v = symbols('m,n,v')
    return solve(Eq(m , n/v).subs({m:M, n:N, v:V}))[0]

def delta_E(n0, n1):
    return -R*h*c*(1/n1**2 - 1/n0**2)

def emission_spectrum(n0, n1):
    wl = abs(h*c*1e+9/delta_E(n0,n1))
    if wl <= 400:
        region = "ultraviolet"
    elif wl >= 780:
        region = "infrared"
    else:
        region = "visible"
    return (wl, region)

def En(n):
    return -2.179e-18/n**2

def C2K(C):
    return C+273.15

def K2C(K):
    return K - 273.15

def psi2atm(psi):
    return psi/14.7

def torr2atm(torr):
    return torr/760

def kpa2atm(kpa):
    return kpa/101.3

def bar2atm(bar):
    return bar/1.013

def atm2kpa(atm):
    return atm*101.3

def atm2bar(atm):
    return atm*1.013

def atm2torr(atm):
    return atm*760

def atm2psi(atm):
    return atm*14.7

def J2cal(J):
    return J/4.184

def cal2J(cal):
    return cal*4.184
    
def igas(p, v, N, t):
    P, V, n, T = sp.symbols('P,V,n,T')
    eq = sp.Eq(P*V , n*R*T).subs({P:p, V:v, n:N, T:t})
    return sp.solve(eq)[0]

def b_law(p1, v1, p2, v2):
    P1, V1, P2, V2 = sp.symbols('P1,V1,P2,V2')
    eq = sp.Eq(P1*V1 , P2*V2).subs({P1:p1,V1:v1,P2:p2,V2:v2})
    return sp.solve(eq)[0]

def c_law(p1, v1, t1, p2, v2, t2):    
    P1, T1, V1, P2, T2, V2 = sp.symbols('P1,T1,V1,P2,T2,V2')
    eq = sp.Eq((P1*V1)/T1 , (P2*V2)/T2).subs({P1:p1,T1:t1,V1:v1,P2:p2,T2:t2,V2:v2})
    return sp.solve(eq)[0]

def Vrms(T, M):
    return sp.sqrt((3*8.3145*T)/(M/1000))

def Vave(T, M):
    return sp.sqrt((8*8.3145*T)/(3.14159*(M/1000)))

def effusion(r1, m1, r2, m2):
    R1, M1, R2, M2 = sp.symbols('R1, M1, R2, M2')
    eq = sp.Eq(R1/R2 , sp.sqrt(M2/M1)).subs({R1:r1, R2:r2, M1:m1, M2:m2})
    return sp.solve(eq)[0]

def density(D, M, V):
    d, m, v = symbols('d,m,v')
    return solve(Eq(d, m/v).subs({d:D, m:M, v:V}))[0]

def dilution(M1, V1, M2, V2):
    m1, v1, m2, v2 = symbols('m1,v1,m2,v2')
    return solve(Eq(m1*v1 , m2*v2).subs({m1:M1, m2:M2, v1:V1, v2:V2}))[0]

def pH(H3O):
    return float(-log(H3O, 10))

def energy(E, W, H):
    e,w,h = symbols('e,w,h')
    return solve(Eq(e, w+h).subs({e:E,w:W,h:H}))[0]

