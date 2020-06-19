#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
provides colsol function to solve finite element static equilibrium equations 
  in core using the active column solver

Created on Fri Jun 19 13:12:36 2020

@author: xzhang@tsinghua.edu.cn
"""

def colsol(n,m,K,R):
    """
    colsol - Solve finite element static equilibrium equations in core using
      the active column solver
    
    Usage:
        IERR, K, R = colsol(n, m, K, R)
    Input parameters
       n      - Number of equations
       m[n]   - Define the skyline of the stiffness matrix K
                The row number (1-based indexing) of the first nonzero 
                element in each column
       K[n,n] - The stiffness matrix
       R[n]   - Right-hand-side load vector
    
    Output parameters
       IERR   - Error indicator. If IERR > 0, K is not positive definite.
       K      - D and L (Factors of stiffness matrix)
                The elements of D are stored on its diagonal, 
                and L replaces its upper triangular part
       R      - Displacement vector
    """
    
    IERR = 0

    # Perform L*D*L(T) factorization of stiffness matrix
    # LDLT is an active column solver to obtain the LDLT factorization
    # of a stiffness matrix K
    # Note tath all indices are zero-based in Python
    for j in range(1, n):  
        for i in range(m[j]+1, j):
            c = 0.0
            for r in range(max(m[i],m[j]), i):
                c += K[r,i]*K[r,j]

            K[i,j] -= c
        
        for i in range(m[j], j):
            Lij = K[i,j]/K[i,i]
            K[j,j] = K[j,j] - Lij*K[i,j]
            K[i,j] = Lij
        
            if K[j,j] <= 0:
                print('Error - stiffness matrix is not positive definite !')
                print('        Nonpositive pivot for equation ', j)
                print('        Pivot = ', K[j,j])
                
                IERR = n
                return IERR, K, R
            
    # Reduce right-hand-side load vector
    for i in range(1, n):
        for j in range(m[i], i):
            R[i] -= K[j,i] * R[j]
    
    # Back-substitute
    for i in range(0, n):
        R[i] /= K[i,i]

    for i in range(n-1, 0, -1):
        for j in range(m[i], i):
            R[j] -= K[j,i]*R[i]

    return IERR, K, R