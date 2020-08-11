#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides methods to setup ID and LM matrices, to calculate the element shape
function matrix and its derivative, element stiffness matrix and nodal force
vector, to setup ID and LM matrix, to impose natural boundary conditions.

Created on Aug. 11 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

import numpy as np

import FEData as model


def setup_ID_LM():
	""" Setup ID and LM arrays """
	count = 0
	count1 = 0

	# Reorder the D.O.F. to make the essential B.C. numbered first
	for i in range(model.neq):
		if model.flags[i] == 2:		# Essential boundary node
			count += 1
			model.ID[i] = count		# The reordered number of essential B.C
			model.d[count] = model.e_bc[i]
		else:
			count1 += 1
			model.ID[i] = model.nd + count1

	for i in range(model.nel):
		for j in range(model.nen):
			for k in range(model.ndof):
				ind = (j - 1)*model.ndof + k
				model.LM[ind, i] = model.ID[model.ndof*(model.IEN[j,i] - 1) + k]

def Nmatrix1D(s, xe):
	"""
	Calculate element shape function matrix N at parent coordinate s

	Args:
		s : The parent coordinate where N is calculated
		xe : (numnp.array(nen)) Element nodal coordinates

	Returns:
		Element shape function matrix N
	"""
	L = xe[1] - xe[0]
	N = np.array([[
		0.25 * (1 - s)^2 * (2 + s),
		0.125 * L * (1 - s)^2 * (1 + s),
		0.25 * (1 + s)^2 * (2 - s),
		0.125 * L * (1 + s)^2 * (s - 1)
	]])
	return N

def Bmatrix1D(s, xe):
	"""
	Calculate derivative of element shape function matrix B at parent coordinate s

	Args:
		s : parent coordinate where B is calculated
		xe : (numnp.array(nen)) Element nodal coordinates

	Returns:
		Derivative of element shape function matrix B
	"""
	L = xe[1] - xe[0]
	B = np.array([[
		1.5*s,
		L*(0.75*s-0.25),
		-1.5*s,
		L*(0.75*s+0.25)
	]])
	return B

def BeamElem(e):
	pass