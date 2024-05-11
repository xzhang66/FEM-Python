#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides methods to calculate element stiffness matrix

Created on Fri Apr 19 18:56:57 2024

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""
import FEData as model
from utitls import gauss
import numpy as np


def SolidShellElem(e):
	"""
	Calculate element stiffness matrix and element nodal body force vector

	Args:
		e : (int) element number

	Returns: ke, fe
		ke : (numpy(nen,nen)) element stiffness matrix
		fe : (numpy(nen,1)) element nodal force vector
	"""
	ke = np.zeros((model.nen*model.ndof, model.nen*model.ndof))
	fe = np.zeros((model.nen*model.ndof, 1))
	
	kaa = np.zeros((1, 1))
	kad = np.zeros((1, model.nen*model.ndof))

	# get coordinates of element nodes
	je = model.IEN[:, e] - 1
	C = np.array([model.x[je], model.y[je], model.z[je]]).T

	ngpf = 2  # No. of gauss points for force calculation
	if model.ngp == 2:	# Full integration
		ngpk = 2   # No. of gauss points for stiffness calculation
	else:
		print("Error : Invalid value of ngp ({}) !".format(model.ngp))

	temp1, detJ0, J0 = DNDxmatSolidShell(0, 0, 0, C)
	T0 = Transmat(J0)
	
	B23a,B23b,B23c,B23d,B13e,B13f,B13g,B13h,B33j,B33k,B33m,B33l = BximatANSSolidShell(C)

	# get gauss points and weights
	wk, gpk = gauss(abs(ngpk))
	wf, gpf = gauss(abs(ngpf))

	# compute element stiffness matrix
	for i in range(ngpk):
		for j in range(ngpk):
			for k in range(ngpk):
				xi = gpk[i]
				eta = gpk[j]
				zeta = gpk[k]
				
				Btmp = BximatSolidShell(xi, eta, zeta, C)

				if model.ANS == 1:
					B13tilde = 0.25*(1+eta)*(1+zeta)*B13e + 0.25*(1-eta)*(1+zeta)*B13f + 0.25*(1-eta)*(1-zeta)*B13g + 0.25*(1+eta)*(1-zeta)*B13h
					Btmp[4,:] = np.copy(B13tilde)
					
					B23tilde = 0.25*(1+xi)*(1+zeta)*B23a + 0.25*(1-xi)*(1+zeta)*B23b + 0.25*(1-xi)*(1-zeta)*B23c + 0.25*(1+xi)*(1-zeta)*B23d
					Btmp[5,:] = np.copy(B23tilde)
					
					B33tilde = 0.25*(1-xi)*(1-eta)*B33j + 0.25*(1+xi)*(1-eta)*B33k + 0.25*(1+xi)*(1+eta)*B33m + 0.25*(1-xi)*(1+eta)*B33l
					Btmp[2,:] = np.copy(B33tilde)
					
				temp2, detJ, JG = DNDxmatSolidShell(xi, eta, zeta, C)
				TG = Transmat(JG)
					
				B =TG@Btmp
				
				# element stiffness matrix
				ke = ke + wk[i]*wk[j]*wk[k]*detJ*(B.T@model.D@B)
				
				if model.EAS == 1:
					Ba_tmp = np.array([[0],[0],[zeta],[0],[0],[0]])
					Ba = detJ0/detJ * (T0 @ Ba_tmp);
					
					kaa = kaa + wk[i]*wk[j]*wk[k]*detJ*(Ba.T@model.D@Ba)
					kad = kad + wk[i]*wk[j]*wk[k]*detJ*(Ba.T@model.D@B)

	if model.EAS == 1:
		ke = ke - 1.0/kaa*(kad.T@kad)

	# compute element nodal force vector
	for i in range(ngpk):
		for j in range(ngpk):
			for k in range(ngpk):
				xi = gpk[i]
				eta = gpk[j]
				zeta = gpk[k]
			
				# shape functions matrix
				N = NmatSolidShell(xi, eta, zeta)
				# derivative of the shape functions
				temp3, detJ, temp4 = DNDxmatSolidShell(xi, eta, zeta, C)

				# element nodal force vector
				be = N@(model.b[:, e].reshape((-1, 1)))
				fe = fe + wf[i]*wf[j]*wf[k]*detJ*(N.T@be)

	return ke, fe


def NmatSolidShell(xi, eta, zeta):
	"""
	Calculate element shape function matrix N at coordinate xt

	Args:
		xi   : The first parent coordinate
		eta  : The second parent coordinate
		zeta : The third parent coordinate

	Returns:
		Element shape function matrix N
	"""

	N1 = (1.0 - xi) * (1.0 - eta) * (1.0 - zeta) / 8.0
	N2 = (1.0 + xi) * (1.0 - eta) * (1.0 - zeta) / 8.0
	N3 = (1.0 + xi) * (1.0 + eta) * (1.0 - zeta) / 8.0
	N4 = (1.0 - xi) * (1.0 + eta) * (1.0 - zeta) / 8.0
	N5 = (1.0 - xi) * (1.0 - eta) * (1.0 + zeta) / 8.0
	N6 = (1.0 + xi) * (1.0 - eta) * (1.0 + zeta) / 8.0
	N7 = (1.0 + xi) * (1.0 + eta) * (1.0 + zeta) / 8.0
	N8 = (1.0 - xi) * (1.0 + eta) * (1.0 + zeta) / 8.0

	return np.array([[N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0, N5, 0, 0, N6, 0, 0, N7, 0, 0, N8, 0, 0],
						 [0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0, N5, 0, 0, N6, 0, 0, N7, 0, 0, N8, 0],
						 [0, 0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0, N5, 0, 0, N6, 0, 0, N7, 0, 0, N8]])


def DNDxmatSolidShell(xi, eta, zeta, C):
	"""
	Calcualte derivative of element shape function matrix DN at coordinate xt

	Args:
		xi   : The first parent coordinate
		eta  : The second parent coordinate
		zeta : The third parent coordinate
		C    : The physical coordinates

	Return:
		Derivative of element shape function matrix
	"""
	dNdxi = 	np.array([[-((eta/2 - 1/2)*(zeta/2 - 1/2))/2, ((eta/2 - 1/2)*(zeta/2 - 1/2))/2, ((eta/2 + 1/2)*(zeta/2 - 1/2))/2, -((eta/2 + 1/2)*(zeta/2 - 1/2))/2, ((eta/2 - 1/2)*(zeta/2 + 1/2))/2, -((eta/2 - 1/2)*(zeta/2 + 1/2))/2, -((eta/2 + 1/2)*(zeta/2 + 1/2))/2, ((eta/2 + 1/2)*(zeta/2 + 1/2))/2],
				[-((xi/2 - 1/2)*(zeta/2 - 1/2))/2, ((xi/2 + 1/2)*(zeta/2 - 1/2))/2, ((xi/2 - 1/2)*(zeta/2 - 1/2))/2, -((xi/2 + 1/2)*(zeta/2 - 1/2))/2, ((xi/2 - 1/2)*(zeta/2 + 1/2))/2, -((xi/2 + 1/2)*(zeta/2 + 1/2))/2, -((xi/2 - 1/2)*(zeta/2 + 1/2))/2, ((xi/2 + 1/2)*(zeta/2 + 1/2))/2],
				[-((xi/2 - 1/2)*(eta/2 - 1/2))/2, ((xi/2 + 1/2)*(eta/2 - 1/2))/2, ((xi/2 - 1/2)*(eta/2 + 1/2))/2, -((xi/2 + 1/2)*(eta/2 + 1/2))/2, ((xi/2 - 1/2)*(eta/2 - 1/2))/2, -((xi/2 + 1/2)*(eta/2 - 1/2))/2, -((xi/2 - 1/2)*(eta/2 + 1/2))/2, ((xi/2 + 1/2)*(eta/2 + 1/2))/2]])
	
	dNdxi[:,[2,3]] = dNdxi[:,[3,2]]
	dNdxi[:,[6,7]] = dNdxi[:,[7,6]]

	# Compute Jacobian matrix
	J = dNdxi@C
	detJ = np.linalg.det(J)

	DNDx = np.linalg.solve(J, dNdxi)

	return DNDx, detJ, J

def BximatSolidShell(xi, eta, zeta, C):
	"""
	Calcualte derivative of element shape function matrix Bxi at coordinate xt

	Args:
		xi   : The first parent coordinate
		eta  : The second parent coordinate
		zeta : The third parent coordinate
		C    : The physical coordinates

	Returns:
		Derivative of element shape function matrix B and Jacobian determination
	"""
	dNdxi = 	np.array([[-((eta/2 - 1/2)*(zeta/2 - 1/2))/2, ((eta/2 - 1/2)*(zeta/2 - 1/2))/2, ((eta/2 + 1/2)*(zeta/2 - 1/2))/2, -((eta/2 + 1/2)*(zeta/2 - 1/2))/2, ((eta/2 - 1/2)*(zeta/2 + 1/2))/2, -((eta/2 - 1/2)*(zeta/2 + 1/2))/2, -((eta/2 + 1/2)*(zeta/2 + 1/2))/2, ((eta/2 + 1/2)*(zeta/2 + 1/2))/2],
				[-((xi/2 - 1/2)*(zeta/2 - 1/2))/2, ((xi/2 + 1/2)*(zeta/2 - 1/2))/2, ((xi/2 - 1/2)*(zeta/2 - 1/2))/2, -((xi/2 + 1/2)*(zeta/2 - 1/2))/2, ((xi/2 - 1/2)*(zeta/2 + 1/2))/2, -((xi/2 + 1/2)*(zeta/2 + 1/2))/2, -((xi/2 - 1/2)*(zeta/2 + 1/2))/2, ((xi/2 + 1/2)*(zeta/2 + 1/2))/2],
				[-((xi/2 - 1/2)*(eta/2 - 1/2))/2, ((xi/2 + 1/2)*(eta/2 - 1/2))/2, ((xi/2 - 1/2)*(eta/2 + 1/2))/2, -((xi/2 + 1/2)*(eta/2 + 1/2))/2, ((xi/2 - 1/2)*(eta/2 - 1/2))/2, -((xi/2 + 1/2)*(eta/2 - 1/2))/2, -((xi/2 - 1/2)*(eta/2 + 1/2))/2, ((xi/2 + 1/2)*(eta/2 + 1/2))/2]])
	
	dNdxi[:,[2,3]] = dNdxi[:,[3,2]]
	dNdxi[:,[6,7]] = dNdxi[:,[7,6]]

	# Compute Jacobian matrix
	J = dNdxi@C
	detJ = np.linalg.det(J)
	
	# [σxx σyy σzz τxy τxz τyz]
	Bxi = np.zeros((6, 24))
	for i in range(3):
		for j in range(8):
			for k in range(3):
				Bxi[i,3*j+k] = dNdxi[i,j]*J[i,k]
	
	for j in range(8):
		for k in range(3):
			Bxi[3,3*j+k] = dNdxi[1,j]*J[0,k] + dNdxi[0,j]*J[1,k]
	
	for j in range(8):
		for k in range(3):
			Bxi[4,3*j+k] = dNdxi[2,j]*J[0,k] + dNdxi[0,j]*J[2,k]
	
	for j in range(8):
		for k in range(3):
			Bxi[5,3*j+k] = dNdxi[1,j]*J[2,k] + dNdxi[2,j]*J[1,k]

	return Bxi

def BximatANSSolidShell(C):
	"""
	Calcualte derivative of element shape function matrix Bxi at coordinate xt

	Args:
		C    : The physical coordinates

	Returns:
		Derivative of element shape function matrix B and Jacobian determination
	"""
	
	Bmat = BximatSolidShell(1,0,1,C)
	B23a = np.copy(Bmat[5,:])
	
	Bmat = BximatSolidShell(-1,0,1,C)
	B23b = np.copy(Bmat[5,:])
	
	Bmat = BximatSolidShell(-1,0,-1,C)
	B23c = np.copy(Bmat[5,:])
	
	Bmat = BximatSolidShell(1,0,-1,C)
	B23d = np.copy(Bmat[5,:])
	
	
	Bmat = BximatSolidShell(0,1,1,C)
	B13e = np.copy(Bmat[4,:])
	
	Bmat = BximatSolidShell(0,-1,1,C)
	B13f = np.copy(Bmat[4,:])
	
	Bmat = BximatSolidShell(0,-1,-1,C)
	B13g = np.copy(Bmat[4,:])
	
	Bmat = BximatSolidShell(0,1,-1,C)
	B13h = np.copy(Bmat[4,:])
	
	
	Bmat = BximatSolidShell(-1,-1,0,C)
	B33j = np.copy(Bmat[2,:])
	
	Bmat = BximatSolidShell(1,-1,0,C)
	B33k = np.copy(Bmat[2,:])
	
	Bmat = BximatSolidShell(1,1,0,C)
	B33m = np.copy(Bmat[2,:])
	
	Bmat = BximatSolidShell(-1,1,0,C)
	B33l = np.copy(Bmat[2,:])
	
	return B23a,B23b,B23c,B23d,B13e,B13f,B13g,B13h,B33j,B33k,B33m,B33l

def Transmat(J):
	"""
	Calcualte coordinate transformation matrix

	Args:
		J   : Jacobian matrix

	Return:
		Coordinate transformation matrix
	"""

	TransMat = np.zeros((6, 6))
	
	I = np.array([[1, 0, 0],
					[0, 1, 0],
					[0, 0, 1]])
	
	invJ = np.linalg.solve(J, I)
	
	TransMat[0,0]=invJ[0,0]*invJ[0,0]
	TransMat[0,1]=invJ[0,1]*invJ[0,1]
	TransMat[0,2]=invJ[0,2]*invJ[0,2]
	TransMat[0,3]=invJ[0,0]*invJ[0,1]
	TransMat[0,4]=invJ[0,1]*invJ[0,2]
	TransMat[0,5]=invJ[0,0]*invJ[0,2]
	TransMat[1,0]=invJ[1,0]*invJ[1,0]
	TransMat[1,1]=invJ[1,1]*invJ[1,1]
	TransMat[1,2]=invJ[1,2]*invJ[1,2]
	TransMat[1,3]=invJ[1,0]*invJ[1,1]
	TransMat[1,4]=invJ[1,1]*invJ[1,2]
	TransMat[1,5]=invJ[1,0]*invJ[1,2]
	TransMat[2,0]=invJ[2,0]*invJ[2,0]
	TransMat[2,1]=invJ[2,1]*invJ[2,1]
	TransMat[2,2]=invJ[2,2]*invJ[2,2]
	TransMat[2,3]=invJ[2,0]*invJ[2,1]
	TransMat[2,4]=invJ[2,1]*invJ[2,2]
	TransMat[2,5]=invJ[2,0]*invJ[2,2]
	TransMat[3,0]=2.0*invJ[0,0]*invJ[1,0]
	TransMat[3,1]=2.0*invJ[0,1]*invJ[1,1]
	TransMat[3,2]=2.0*invJ[0,2]*invJ[1,2]
	TransMat[3,3]=invJ[0,0]*invJ[1,1]+invJ[1,0]*invJ[0,1]
	TransMat[3,4]=invJ[0,1]*invJ[1,2]+invJ[1,1]*invJ[0,2]
	TransMat[3,5]=invJ[0,0]*invJ[1,2]+invJ[1,0]*invJ[0,2]
	TransMat[4,0]=2.0*invJ[1,0]*invJ[2,0]
	TransMat[4,1]=2.0*invJ[1,1]*invJ[2,1]
	TransMat[4,2]=2.0*invJ[1,2]*invJ[2,2]
	TransMat[4,3]=invJ[1,0]*invJ[2,1]+invJ[2,0]*invJ[1,1]
	TransMat[4,4]=invJ[1,1]*invJ[2,2]+invJ[2,1]*invJ[1,2]
	TransMat[4,5]=invJ[1,0]*invJ[2,2]+invJ[2,0]*invJ[1,2]
	TransMat[5,0]=2.0*invJ[2,0]*invJ[0,0]
	TransMat[5,1]=2.0*invJ[2,1]*invJ[0,1]
	TransMat[5,2]=2.0*invJ[2,2]*invJ[0,2]
	TransMat[5,3]=invJ[2,0]*invJ[0,1]+invJ[0,0]*invJ[2,1]
	TransMat[5,4]=invJ[2,1]*invJ[0,2]+invJ[0,1]*invJ[2,2]
	TransMat[5,5]=invJ[2,0]*invJ[0,2]+invJ[0,0]*invJ[2,2]
	
	# [σxx σyy σzz τxy τxz τyz]
	TransMat[:,[4,5]] = TransMat[:,[5,4]]
	TransMat[[4,5],:] = TransMat[[5,4],:]

	return TransMat