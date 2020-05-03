#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides methods to create FE mode for a bar from a json file, to calculate
the element shape function matrix and its derivative, element stiffness matrix
and nodal force vector, to setup ID and LM matrix, to impose natural boundary 
conditions.

Created on Sun Apr 24 18:56:57 2020

@author: xzhang@tsinghua.edu.cn
"""

import numpy as np
import json

import FEData as model
from utitls import gauss


def create_model_json(DataFile):
    """ 
    Initialize the FEM model from file DataFile (in json format)
    """

    # input data from json file
    with open(DataFile) as f_obj:
        FEData = json.load(f_obj)

    model.Title= FEData['Title']
    model.nsd  = FEData['nsd']
    model.ndof = FEData['ndof']
    model.nnp  = FEData['nnp']
    model.nel  = FEData['nel']
    model.nen  = FEData['nen']    
    model.neq  = model.ndof*model.nnp

    # initialize K, d and f 
    model.f = np.zeros((model.neq,1))            
    model.d = np.zeros((model.neq,1))        
    model.K = np.zeros((model.neq,model.neq))    

    # element and material data (given at the element nodes)
    model.E     = np.array(FEData['E'])
    model.body  = np.array(FEData['body'])
    model.CArea = np.array(FEData['CArea'])

    # gauss integration
    model.ngp = FEData['ngp']

    # boundary conditions
    model.flags = np.array(FEData['flags'])
    model.e_bc = np.array(FEData['e_bc'])
    model.n_bc = np.array(FEData['n_bc'])
    model.nd = FEData['nd']

    # point forces
    model.np = FEData['np']
    if model.np > 0:
        model.xp = np.array(FEData['xp'])
        model.P  = np.array(FEData['P'])

    # output plots
    model.plot_bar = FEData['plot_bar']
    model.plot_nod = FEData['plot_nod']
    model.nplot = model.nnp*10
    model.plot_tex = FEData['plot_tex']

    # define the mesh
    model.x = np.array(FEData['x'])
    model.y = np.array(FEData['y'])  
    model.IEN = np.array(FEData['IEN'])

    model.ID  = np.zeros(model.neq,np.int)
    model.LM  = np.zeros((model.nen,model.nel),np.int)   

    # generate LM and ID arrays
    setup_ID_LM()


def setup_ID_LM():
    """ Setup ID and LM arrays """
    count  = 0
    count1 = 0
    
    # Reorder the D.O.F. to make the essential B.C. numbered first
    for i in range(model.neq):
        if model.flags[i] == 2:   # Essential boundary node
            count = count + 1      
            model.ID[i] = count   # The reordered number of essential B.C
            model.d[count] = model.e_bc[i] 
        else:
            count1 = count1 + 1
            model.ID[i] = model.nd + count1

    for i in range(model.nel):
        for j in range(model.nen):
            model.LM[j,i] = model.ID[model.IEN[j,i]-1]


def naturalBC():
    """ Compute and assemble nodal boundary force vector """
    for i in range(model.nnp):
       if model.flags[i] == 1:
          node = model.ID[i]-1
          model.f[node] += model.CArea[node]*model.n_bc[node]


def Nmatrix1D(xt, xe):
    """ 
    Calculate element shape function matrix N at coordinate xt
    
    Args:
        xt : The physical coordinate where N is calculated
        xe : (numnp.array(nen)) Element nodal coordinates
        
    Returns:
        Element shape function matrix N
    """
    if model.nen == 2:  # linear shape functions
        N = np.array([[
            (xt-xe[1])/(xe[0]-xe[1]),
            (xt-xe[0])/(xe[1]-xe[0]) ]])
    elif model.nen == 3:  # quadratic shape functions
        N = np.array([[
            (xt-xe[1])*(xt-xe[2])/((xe[0]-xe[1])*(xe[0]-xe[2])),
            (xt-xe[0])*(xt-xe[2])/((xe[1]-xe[0])*(xe[1]-xe[2])),
            (xt-xe[0])*(xt-xe[1])/((xe[2]-xe[0])*(xe[2]-xe[1])) ]])
    return N


def Bmatrix1D(xt, xe):
    """
    Calcualte derivative of element shape function matrix B at coordinate xt
    
    Args:
        xt : Physical coordinate where B is calculated
        xe : (numnp.array(nen)) Element nodal coordinates
        
    Returns:
        Derivative of element shape function matrix B
    """
    if model.nen == 2:   # derivative of linear shape functions (constant)
        B = 1/(xe[1]-xe[0])*np.array([[-1,1]])
    elif model.nen == 3: # derivative of quadratic shape functions
        B = np.array([[
            (2*xt-xe[1]-xe[2])/((xe[0]-xe[1])*(xe[0]-xe[2])),
            (2*xt-xe[0]-xe[2])/((xe[1]-xe[0])*(xe[1]-xe[2])),
            (2*xt-xe[0]-xe[1])/((xe[2]-xe[0])*(xe[2]-xe[1])) ]])
    return B


def BarElem(e):
    """ 
    Calculate element stiffness matrix and element nodal body force vector
    
    Args:
        e : (int) element number
        
    Returns: ke, fe
        ke : (numpy(nen,nen)) element stiffness matrix
        fe : (numpy(nen,1)) element nodal force vector
    """

    IENe = model.IEN[:,e]-1  # extract local connectivity information
    xe = model.x[IENe]       # extract element x coordinates
    J = (xe[model.nen-1] - xe[0])/2   # compute Jacobian
    [w , gp] = gauss(model.ngp)       # extract Gauss points and weights

    ke = np.zeros((model.nen,model.nen))  # initialize element stiffness matrix
    fe = np.zeros((model.nen,1))    # initialize element nodal force vector

    for i in range(model.ngp):
        # Compute Gauss points in physical coordinates
        xt = 0.5*(xe[0]+xe[-1])+J*gp[i]
        
        # Calculate the element shape function matrix and its derivative
        N = Nmatrix1D(xt,xe)
        B = Bmatrix1D(xt,xe)

        # cross-sectional area,Young's modulus and body force at gauss points
        Ae = N@model.CArea[IENe]
        Ee = N@model.E[IENe]
        be = N@model.body[IENe]

        # compute element stiffness matrix and nodal body force vector
        ke = ke + w[i]*Ae*Ee*(B.T@B)
        fe =  fe + w[i]*N.T*be 

    ke = J*ke
    fe = J*fe

    # check for point forces in this element
    for i in range(model.np):   # loop over all point forces
        Pi = model.P[i]    # extract point force
        xpi = model.xp[i]   # extract the location of point force within an element
        if xe[0] <= xpi & xpi < xe[-1]:
            # add to the nodal force vector
            fe = fe + Pi*np.transpose(Nmatrix1D(xpi,xe))   

    return ke,fe