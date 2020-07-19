#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exact solutions for some problems
  ExactSolution_TaperedBar: Plot the exact solution of the tapered elastic bar
  given in Example 5.2 in Fish's textbook
  ExactSolution_CompressionmBar: Plot the exact solution of the bar under
  compression given in Figure 5.13 in Fish's textbook

Created on Wed Apr 29 10:26:04 2020

@author: xzhang@tsinghua.edu.cn
"""

import numpy as np
from math import log, sqrt

from utitls import gauss
from Bar1DElem import Nmatrix1D, Bmatrix1D
import FEData as model

    
def ExactSolution_TaperedBar(ax1, ax2):
    """ 
    Plots the exact displacement and stress of the tapered elastic bar
    in axes ax1 and ax2, respectively.
    
    Args:
        ax1 : axis to draw displacement distribution
        ax2 : axis to draw stress distribution
    """
 
    # divide the problem domain into two regions 
    xa = np.arange(2,5,0.01)
    xb = np.arange(5,6,0.01)
         
    # exact displacement for xa 
    c1 = 72;  c2 = 1 - (c1/16)*log(2)
    u1 = -.5*xa + (c1/16)*np.log(xa) + c2
    
    # exact displacement for xb 
    c3 = 48;  c4 = log(5)/16*(c1-c3) + c2
    u2 = -.5*xb + (c3/16)*np.log(xb) + c4
    
    # plot displacement 
    ax1.plot(np.append(xa,xb),np.append(u1,u2), '--r', label='Exact')
     
    # exact stress for xa 
    ya = (36-4*xa)/xa
    
    # exact stress for xb 
    yb = (24-4*xb)/xb
    
    # plot stress 
    ax2.plot(np.append(xa,xb),np.append(ya,yb), '--r', label='Exact')
    
    
def ExactSolution_CompressionBar(ax1, ax2):
    """ 
    Plots the exact displacement and stress of a elastic bar under compression
    in axes ax1 and ax2, respectively.
    
    Args:
        ax1 : axis to draw displacement distribution
        ax2 : axis to draw stress distribution
    """
    xx = np.arange(0, 2, 0.01)

    # exact displacement for a bar under compression
    Ee = 10000
    ue = (-xx**3/6 + xx)/Ee 
    
    # plot displacement 
    ax1.plot(xx, ue, '--r',  label='Exact')
    
    # exact stress
    stre = (-xx**2/2 + 1)
    
    # plot stress 
    ax2.plot(xx,stre, '--r', label='Exact')


def ErrorNorm_CompressionBar():
    """ 
    Calculate and print the error norm (L2 and energy norm) of the elastic 
    bar under compression for convergence study
    """
    
    ngp = 3
    [w, gp] = gauss(ngp)    # extract Gauss points and weights
    
    L2Norm = 0
    EnNorm = 0
    
    L2NormEx = 0
    EnNormEx = 0
    
    for e in range(model.nel):
        
        de = model.d[model.LM[:,e]-1] # extract element nodal displacements
        IENe = model.IEN[:,e]-1       # extract local connectivity information
        xe = model.x[IENe]            # extract element x coordinates
        J = (xe[-1] - xe[0])/2        # compute Jacobian
        
        for i in range(ngp):
            xt = 0.5*(xe[0]+xe[-1])+J*gp[i]  # Gauss points in physical coordinates
            
            N = Nmatrix1D(xt,xe)     # shape functions matrix
            B = Bmatrix1D(xt,xe)     # derivative of shape functions matrix
            
            Ee = N@model.E[IENe]     # Young's modulus at element gauss points
            
            uh  = N@de               # displacement at gauss point
            uex = (-xt**3/6 + xt)/Ee # Exact displacement
            L2Norm += J*w[i]*(uex - uh)**2
            L2NormEx += J*w[i]*(uex)**2
            
            sh  = B@de               # strain at Gauss points
            sex = (-xt**2/2 + 1)/Ee  # Exact strain
            EnNorm += 0.5*J*w[i]*Ee*(sex-sh)**2
            EnNormEx += 0.5*J*w[i]*Ee*(sex)**2
    
    L2Norm = sqrt(L2Norm)
    L2NormEx = sqrt(L2NormEx)
    
    EnNorm = sqrt(EnNorm)
    EnNormEx = sqrt(EnNormEx)
    
    # print stresses at element gauss points
    print('\nError norms')
    print('%13s %13s %13s %13s %13s'
          %('h','L2Norm','L2NormRel','EnNorm','EnNormRel'))
    print('%13.6E %13.6E %13.6E %13.6E %13.6E\n'
          %(2/model.nel, L2Norm, L2Norm/L2NormEx, EnNorm, EnNorm/EnNormEx))
    
    return 2/model.nel, L2Norm, EnNorm
