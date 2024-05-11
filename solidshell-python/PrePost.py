#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Provide methods to create FE model from a json file, to plot the mesh,
to print stresses at nodes, to plot displacement and stress distributions
obtained by FE analysis.

Created on Fri Apr 19 18:56:57 2024

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""
import json
import FEData as model
import numpy as np


def create_model_json(DataFile):
	"""
	Initialize the FEM model from file DataFile (in json format)
	"""

	with open(DataFile) as f_obj:
		FEData = json.load(f_obj)

	model.Title = FEData['Title']
	model.nsd = FEData['nsd']
	model.ndof = FEData['ndof']
	model.nnp = FEData['nnp']
	model.nel = FEData['nel']
	model.nen = FEData['nen']
	model.neq = model.ndof * model.nnp

	# initialize K, d and f
	model.f = np.zeros((model.neq, 1))
	model.d = np.zeros((model.neq, 1))
	model.K = np.zeros((model.neq, model.neq))

	# material properties
	E = FEData['E']
	ne = FEData['nu']

	model.D = np.array([[1-ne, ne, ne, 0, 0, 0],
						[ne, 1-ne, ne, 0, 0, 0],
						[ne, ne, 1-ne, 0, 0, 0],
						[0, 0, 0, 0.5-ne, 0, 0],
						[0, 0, 0, 0, 0.5-ne, 0],
						[0, 0, 0, 0, 0, 0.5-ne]])*E/(1 + ne)/(1 - 2*ne)
	model.G = E / (2.0 * (1.0 + ne))
	
	model.ANS = FEData['ANS']
	model.EAS = FEData['EAS']

	# gauss integration
	model.ngp = FEData['ngp']

	# boundary conditions
	model.flags = np.array(FEData['flags'])
	model.nd = FEData['nd']
	if model.nbe > 0:
		model.n_bc = np.array(FEData['n_bc'])

	# The Essential B.C. is set to zero by default
	try:
		model.e_bc = np.array(FEData['e_bc'])
	except KeyError:
		model.e_bc = np.zeros((model.neq, 1))

	# force conditions
	# The F.C. is set to zero by default
	try:
		model.P = np.array(FEData['P'])
	except KeyError:
		model.P = np.zeros((model.neq, 1))

	try:
		model.b = np.array(FEData['b'])
	except KeyError:
		model.b = np.zeros((model.nen*model.ndof, model.nel))

	# define the mesh
	model.x = np.array(FEData['x'])
	model.y = np.array(FEData['y'])
	model.z = np.array(FEData['z'])
	model.IEN = np.array(FEData['IEN'], dtype=int)

	model.print_disp = FEData['print_disp']

	model.ID = np.zeros(model.neq, dtype=int)
	model.LM = np.zeros((model.nen*model.ndof, model.nel), dtype=int)
	setup_ID_LM()


def point_and_trac():
	"""
	Add the nodal forces to the global force vector.
	"""
	# Assemble point forces
	model.f[model.ID - 1] = model.f[model.ID - 1] + model.P


def setup_ID_LM():
	"""
	Calculate the ID and LM matrix according to model.flags and model.IEN matrix.
	"""
	count = 0
	count1 = 0
	for i in range(model.neq):
		if model.flags[i] == 2:
			# check if a node on essential boundary
			count += 1
			model.ID[i] = count
			model.d[count-1] = model.e_bc[i]
		else:
			count1 += 1
			model.ID[i] = model.nd + count1

	for i in range(model.nel):
		n = 0
		for j in range(model.nen):
			blk = model.ndof * (model.IEN[j, i] - 1)
			for k in range(model.ndof):
				model.LM[n, i] = model.ID[blk + k]
				n += 1

def postprocess():
	"""
	1. Calculate the coordinates of deformed configuration.
	2. Plot the initial and deformed configuration in one figure.
	3. Print the element stress on Gauss Point.
	4. Calculate the nodal stress and plot the stress contours.
	"""
	
	print('No. of Elements  {}'.format(model.nel))
	print('No. of Nodes     {}'.format(model.nnp))
	print('No. of Equations {}'.format(model.neq))
	
	# print the nodal displacement
	print_displacement()

def print_displacement():
	"""
	Print the displacement of all nodes.
	"""
	if model.print_disp == 'yes':
		dis = model.d[model.ID - 1]

		print("\n                           Nodal displacement")
		print("-------------------------------------------------------------------------------")
		print("\tnode\tx\ty\tz\t\t\tu_x\t\t\tu_y\t\t\tu_z")

		for i in range(model.nnp):
			print("\t{}\t{}\t{}\t{}\t{:.15e}\t{:.15e}\t{:.15e}".format(i+1, model.x[i], model.y[i], model.z[i], dis[3*i], dis[3*i+1], dis[3*i+2]))


