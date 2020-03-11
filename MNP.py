#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:28:38 2020

@author: Arthur
"""

from sympy import *
from sympy.physics.mechanics import *
from sympy.utilities.codegen import codegen
from pydy.codegen.ode_function_generators import generate_ode_function
from pydy.codegen.octave_code import OctaveMatrixGenerator

#coordenadas generalizadas
M1x , M1y , M1z  = dynamicsymbols('M1x  M1y  M1z')
M1vx, M1vy, M1vz = dynamicsymbols('M1vx M1vy M1vz')
M1xd,  M1yd,  M1zd  = dynamicsymbols('M1x  M1y  M1z',1)
M1vxd, M1vyd, M1vzd = dynamicsymbols('M1vx M1vy M1vz', 1)

#derivadas
M1phi, M1psi, M1sig = dynamicsymbols('M1phi M1psi M1sig')
M1vphi, M1vpsi, M1vsig = dynamicsymbols('M1vphi M1vpsi M1vsig')
M1phid, M1psid, M1sigd = dynamicsymbols('M1phi M1psi M1sig',1)
M1vphid, M1vpsid, M1vsigd = dynamicsymbols('M1vphi M1vpsi M1vsig',1)

#constantes
m, l, g, k, r, d = symbols('m l g k r d')

#sistema de coordenadas inercial
N = ReferenceFrame('N')
O = Point('O') #origem
O.set_vel(N,0)

ax, ay = symbols('ax ay')

Mref = ReferenceFrame('Mref')
Mref.orient(N, 'space', (M1phi, M1psi, M1sig), 'xyz')

P = Point('P')
P.set_pos(O, M1x*N.x + M1y*N.y + M1z*N.z)
P.set_vel(N, M1vx*N.x + M1vy*N.y + M1vz*N.z)

Ixx = (m*r**2)/2
Iyy = m*(3*r**2 + d**2)/12
I = inertia(N, Ixx, Iyy, Iyy)
inertia_tuple = (I, P)
M1 = RigidBody('M1', P, Mref, m, inertia_tuple)

# definindo as molas
# 4 pontos fixos:
A1 = O.locatenew('A1', ax*N.x + ay*N.y)
A2 = O.locatenew('A2',-ax*N.x + ay*N.y)
A3 = O.locatenew('A3',-ax*N.x - ay*N.y)
A4 = O.locatenew('A4', ax*N.x - ay*N.y)

Apts = [A1, A2, A3, A4]

for pt in Apts:
    pt.set_vel(N, 0)

bx, by, bz = symbols('bx by bz')

B1 = P.locatenew('B1', bx*Mref.x + by*Mref.y + bz*Mref.z)
B2 = P.locatenew('B2',-bx*Mref.x + by*Mref.y + bz*Mref.z)
B3 = P.locatenew('B3',-bx*Mref.x - by*Mref.y + bz*Mref.z)
B4 = P.locatenew('B4', bx*Mref.x - by*Mref.y + bz*Mref.z)

Bpts = [B1, B2, B3, B4]

for pt in Bpts:
    pt.v2pt_theory(P,N,Mref)

#forces
forces = []
for ap, bp in zip(Apts, Bpts):
    forces.append((bp, -k*(bp.pos_from(ap).magnitude() - l)*(bp.pos_from(ap).normalize())))

forces.append((P, -m*g*N.z))

#kinematic eqs
kinematic_differentials = [M1vx - M1xd,
                           M1vy - M1yd,
                           M1vz - M1zd,
                           M1vphi - M1phid,
                           M1vpsi - M1psid,
                           M1vsig - M1sigd]
coords = [M1x, M1y, M1z, M1phi, M1psi, M1sig]
velocs = [M1vx, M1vy, M1vz, M1vphi, M1vpsi, M1vsig]
kane = KanesMethod(N, q_ind=coords, u_ind=velocs, kd_eqs=kinematic_differentials)
(fr, frst) = kane.kanes_equations([M1], loads = forces)


#eqs = kane.mass_matrix.inv() * kane.forcing
#codegen(('rigid_pend_1',eqs),language='Octave', to_files = True)

#eqs = kane.rhs()
#eqs = eqs.subs([(M1y,0),(M1sig,0),(M1phi,0)])
#eqs = eqs.subs([(cos(M1psi),1-(1/2)*M1psi**2),(sin(M1psi),M1psi)])
#eqs = eqs.subs([(ay,0),(by,0)])
#eqs = simplify(eqs)
#eqs = eqs.subs([(bx, ax)])
#eqs = eqs.subs(bz,0)