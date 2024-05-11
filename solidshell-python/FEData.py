#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Global variables defining the FEM model
  Title: (str) Title of the problem to be solved.
  nsd  : (int) Number of space dimensions.
  ndof : (int) Number of degrees-of-freedom per node.
  nnp  : (int) Number of nodal points.
  nel  : (int) Number of elements.
  nen  : (int) Number of element nodes.
  neq  : (int) Number of equations (D.O.F)
  ngp  : (int) Number of gauss points.
         2 - Full integration
  
  flags: (numpy.array(nnp))  Nodal boundary condition flag:
         2 - located on the essential boundary;
         1 - located on the natural boundary.
  nd   : (int) Number of nodes on the essential boundary.
  e_bc : (numpy.array(nnp)) Value of essential B.C.

  nbe  : (int) Number of prescribed traction edges.
  n_bc : (numpy.array(nnp)) Value of natural B.C.
  
  P    : (numpy.array(neq)) Array of nodal external forces.
  b    : (numpy.array(nen*ndof, nel)) Element nodal forces.
  D	   : (numpy.array(3, 3)) elasticity matrix
  G	   : (float) shear modulus
  
  ANS   : (int) Assumed natural strain method ?
  EAS   : (int) Enhanced assumed strain method ?
  
  IEN  : (numpy.array(nen,nel)) Element connectivity array.
  ID   : (numpy.array(neq) Identification matrix.
  LM   : (numpy.array(nen,nel)) Location matrix.
  x    : (numpy.array(nnp))x coordinate.
  y    : (numpy.array(nnp))y coordinates.
  z    : (numpy.array(nnp))z coordinates.

  K    : (numpy.array(neq,neq)) Global stiffness matrix
  f    : (numpy.array(neq,1)) Global nodal force vector
  d    : (numpy.array(neq,1)) Solution vector

  print_disp      : print nodal displacement ?
"""

Title = None
nsd = 0
ndof = 0
nnp = 0
nel = 0
nen = 0
neq = 0
ngp = 2
nd = 0
nbe = 0

f = None
d = None
K = None

# boundary conditions
flags = None
e_bc = None
n_bc = None

# force conditions
P = None
b = None

# material
D = None
G = 0.0

# The methond for locking phenomena
ANS = 0
EAS = 0

# define the mesh
x = None
y = None
z = None
IEN = None

ID = None
LM = None

# parameter for postprocess
print_disp = None
