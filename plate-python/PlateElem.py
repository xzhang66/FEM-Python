#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides methods to calculate element stiffness matrix

Created on Fri Jun 19 18:56:57 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""
import FEData as model
from utitls import gauss
import numpy as np


def PlateElem(e):
	"""
	Calculate element stiffness matrix and element nodal body force vector

	Args:
		e : (int) element number

	Returns: ke, fe
		ke : (numpy(nen*ndof,nen*ndof)) element stiffness matrix
		fe : (numpy(nen*ndof,1)) element nodal force vector
	"""
	ke = np.zeros((model.nen*model.ndof, model.nen*model.ndof))
	fe = np.zeros((model.nen*model.ndof, 1))

	# get coordinates of element nodes
	je = model.IEN[:, e] - 1
	C = np.array([model.x[je], model.y[je]]).T

	# get gauss points and weights
	w, gp = gauss(model.ngp)

	# compute element stiffness matrix and element nodal force vector
	for i in range(model.ngp):
		for j in range(model.ngp):
			eta = gp[i]
			psi = gp[j]
			# shape functions matrix
			N = NmatPlate(eta, psi)
			# derivative of the shape functions
			B, detJ = BmatPlate(eta, psi, C)

			# element stiffness matrix
			ke = ke + w[i]*w[j]*detJ*(B.T@model.D@B)
			be = N@(model.b[:, e].reshape((-1, 1)))
			fe = fe + w[i]*w[j]*detJ*(N.T@be)
	return ke, fe


def NmatPlate(eta, psi):
	"""
	Calculate element shape function matrix N at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate

	Returns:
		Element shape function matrix N
	"""
	# parent coordinates at nodes
	eta_val = np.array([-1, 1, 1, -1])
	psi_val = np.array([-1, -1, 1, 1])
	
	N = np.zeros((1, 12))
	
	for i in range(model.nen):
		N[0,3*i] = 0.125*( 1 +eta_val[i]*eta)*(1 + psi_val[i]*psi) * \
				(2 + eta_val[i]*eta + psi_val[i]*psi - eta**2 - psi**2)
				
		N[0,3*i+1] = 0.125*( 1 +eta_val[i]*eta)*(1 + psi_val[i]*psi) * \
				(-model.be * psi_val[i] * (1 - psi**2))
				
		N[0,3*i+2] = 0.125*( 1 +eta_val[i]*eta)*(1 + psi_val[i]*psi) * \
				(model.ae * eta_val[i] * (1 - eta**2))

	return N


def BmatPlate(eta, psi, C):
	"""
	Calcualte derivative of element shape function matrix B at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate
		C   : The physical coordinates

	Returns:
		Derivative of element shape function matrix B and Jacobian determination
	"""
	# global coordinates at nodes
	x_val = np.array([C[0,0], C[1,0], C[2,0], C[3,0]])
	y_val = np.array([C[0,1], C[1,1], C[2,1], C[3,1]])
	
	#Calculate the B_M matrix
	B_M = np.zeros((12, 12))
	for i in range(model.nen):
		B_M[3*i:3*i+3,:] = np.array([[ 1, x_val[i], y_val[i], \
								x_val[i]**2, x_val[i]*y_val[i], y_val[i]**2, \
								x_val[i]**3, x_val[i]**2*y_val[i], x_val[i]*y_val[i]**2, \
								y_val[i]**3, x_val[i]**3*y_val[i], x_val[i]*y_val[i]**3 ], \
								[ 0, 0, 1, \
								0, x_val[i], 2*y_val[i], \
								0, x_val[i]**2, 2*x_val[i]*y_val[i], \
								3*y_val[i]**2, x_val[i]**3, 3*x_val[i]*y_val[i]**2 ], \
								[ 0, -1, 0, \
								-2*x_val[i], -y_val[i], 0, \
								-3*x_val[i]**2, -2*x_val[i]*y_val[i], -y_val[i]**2, \
								0, -3*x_val[i]**2*y_val[i], -y_val[i]**3 ]])
	
	# global coordinates at (eta, psi)
	xt = eta * model.ae + (x_val[0] + x_val[1]) / 2.0
	yt = psi * model.be + (y_val[1] + y_val[2]) / 2.0
	
	#Calculate the B_Q matrix
	B_Q = np.array([[0, 0, 0, 2, 0, 0, 6*xt, 2*yt, 0, 0, 6*xt*yt, 0], \
					[0, 0, 0, 0, 0, 2, 0, 0, 2*xt, 6*yt, 0, 6*xt*yt], \
					[0, 0, 0, 0, 2, 0, 0, 4*xt, 4*yt, 0, 6*xt**2, 6*yt**2]])
	
	B = B_Q @ np.linalg.inv(B_M)
	
	# Compute Jacobian determination
	detJ = model.ae * model.be

	return B, detJ
