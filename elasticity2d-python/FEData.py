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
  E    : (numpy.array(nnp)) Nodal values Young's modulus.
  body : (numpy.arraay(nnp)) Nodal values body forces.
  CArea: (numpy.array(nnp)) Nodal values of cross-sectional area.
  ngp  : (int) Number of gauss points.
  flags: (numpy.array(nnp))  Nodal boundary condition flag:
         2 - located on the essential boundary;
         1 - located on the natural boundary.
  nd   : (int) Number of nodes on the essential boundary.
  e_bc : (numpy.array(nnp)) Value of essential B.C.
  n_bc : (numpy.array(nnp)) Value of natural B.C.
  np   : (int) Number of point forces.
  xp   : ((numpy.array(np))) Array of coordinates where point forces are applied.
  P    : (numpy.array(np)) Array of point forcess.
  plot_bar: (bool) Plot bar ?
  plot_nod: plot node number ?
  plot_tex : Convert figures into PGFPlots figures in LaTex file ?
  nplot: (int) Number of points in the element to plot displacements 
         and stresses (10*nnp).
  x    : (numpy.array(nnp))x coordinate.
  y    : (numpy.array(nnp))y coordinates, used only for the bar plot.
  IEN  : (numpy.array(nen,nel)) Element connectivity array.
  ID   : (numpy.array(neq) Identification matrix.
  LM   : (numpy.array(nen,nel)) Location matrix.
  K    : (numpy.array(neq,neq)) Global stiffness matrix
  f    : (numpy.array(neq,1)) Global nodal force vector            
  d    : (numpy.array(neq,1)) Solution vector
"""