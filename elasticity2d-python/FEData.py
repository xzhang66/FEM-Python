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
  nd   : (int) Number of nodes on the essential boundary.

  flags: (numpy.array(nnp))  Nodal boundary condition flag:
         2 - located on the essential boundary;
         1 - located on the natural boundary.

  e_bc : (numpy.array(nnp)) Value of essential B.C.
  n_bc : (numpy.array(nnp)) Value of natural B.C.
  P    : (numpy.array(nnp)) Array of point forcess.
  b    : (numpy.array(nnp)) Array of body forcess.
  D	   : (numpy.array(3, 3)) elasticity matrix
  IEN  : (numpy.array(nen,nel)) Element connectivity array.
  ID   : (numpy.array(neq) Identification matrix.
  LM   : (numpy.array(nen,nel)) Location matrix.
  nbe  : (int) number of edges on the boundary
  x    : (numpy.array(nnp))x coordinate.
  y    : (numpy.array(nnp))y coordinates, used only for the bar plot.
  K    : (numpy.array(neq,neq)) Global stiffness matrix
  f    : (numpy.array(neq,1)) Global nodal force vector            
  d    : (numpy.array(neq,1)) Solution vector

  counter         : counter of nodes for stress plots
  nodestress      : stresses at nodes for the stress plots [sxx syy sxy]
  compute_stress  : computer stresses ?

  plot_mesh       : plot mesh ?
  plot_disp       : plot displaced mesh ?
  plot_nod        : plot node number ?
  plot_stress_xx  : plot stress xx ?
  plot_mises      : plot mises stress ?
  fact            : factor for scaled displacements plot
"""

Title = None
nsd = 0
ndof = 0
nnp = 0
nel = 0
nen = 0
neq = 0
ngp = 0
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

# define the mesh
x = None
y = None
IEN = None

ID = None
LM = None

# parameter for postprocess
counter = None
nodestress = None
compute_stress = None
plot_mesh = None
plot_disp = None
plot_nod = None
plot_stress_xx = None
plot_mises = None
fact = 1
