#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Provide methods to create FE model from a json file, to plot the mesh.

Created on Fri Jun 19 18:56:57 2020

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

import json
import FEData as model
import numpy as np
import matplotlib.pyplot as plt
from utitls import gauss
import tikzplotlib


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
	model.nbe = FEData['nbe']
	model.neq = model.ndof * model.nnp

	# initialize K, d and f
	model.f = np.zeros((model.neq, 1))
	model.d = np.zeros((model.neq, 1))
	model.K = np.zeros((model.neq, model.neq))

	# geometric data
	model.lx = FEData['lx']
	model.ly = FEData['ly']
	model.nelx = FEData['nelx']
	model.nely = FEData['nely']
	model.nenx = model.nelx + 1
	model.neny = model.nely + 1
	model.ae = model.lx / (2 * model.nelx)
	model.be = model.ly / (2 * model.nely)
	if model.nelx % 2 != 0:
		print('No. of Elements  {}  is not even, can not get the center deflection'.format(model.nelx))
	if model.nely % 2 != 0:
		print('No. of Elements  {}  is not even, can not get the center deflection'.format(model.nely))

	# material properties
	model.E = FEData['E']
	model.ne = FEData['nu']
	model.G = model.E / (2.0 * (1.0 + model.ne))
	model.r = FEData['r']

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

	try:
		model.q = FEData['q']
	except KeyError:
		model.q = 0.0
	
	# Add the uniform load
	for e in range(model.nel):
		for i in range(model.nen):
			model.b[i*model.ndof:i*model.ndof+model.ndof, e] +=np.array([0, 0, model.q])

	# define the mesh
	model.x = np.array(FEData['x'])
	model.y = np.array(FEData['y'])
	model.IEN = np.array(FEData['IEN'], dtype=np.int)

	# parameter for postprocess
	model.plot_mesh = FEData['plot_mesh']
	model.plot_nod = FEData['plot_nod']
	model.plot_tex = FEData['plot_tex']

	plot_mesh()

	model.ID = np.zeros(model.neq, dtype=np.int)
	model.LM = np.zeros((model.nen*model.ndof, model.nel), dtype=np.int)
	setup_ID_LM()


def point_and_trac():
	"""
	Add the nodal forces and natural B.C. to the global force vector.
	"""
	# Assemble point forces
	model.f[model.ID - 1] = model.f[model.ID - 1] + model.P

	# Compute nodal boundary force vector
	for i in range(model.nbe):
		ft = np.zeros((4, 1))							# initialize nodal boundary force vector
		node1 = int(model.n_bc[0, i])					# first node
		node2 = int(model.n_bc[1, i])					# second node
		n_bce = model.n_bc[2:, i].reshape((-1, 1))		# traction value at node1

		# coordinates
		x1 = model.x[node1 - 1]
		y1 = model.y[node1 - 1]
		x2 = model.x[node2 - 1]
		y2 = model.y[node2 - 1]

		# edge length
		leng = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
		J = leng/2.0

		ngp = abs(model.ngp)
		w, gp = gauss(ngp)

		for j in range(ngp):
			psi = gp[j]
			N = 0.5*np.array([[1-psi, 0, 1+psi, 0],
							  [0, 1-psi, 0, 1+psi]])

			traction = N@n_bce
			ft = ft + w[j]*J*(N.T@traction)

		# Assemble nodal boundary force vector
		ind1 = model.ndof*(node1 - 1)
		ind2 = model.ndof*(node2 - 1)

		model.f[model.ID[ind1] - 1, 0] += ft[0]
		model.f[model.ID[ind1 + 1] - 1, 0] += ft[1]
		model.f[model.ID[ind2] - 1, 0] += ft[2]
		model.f[model.ID[ind2 + 1] - 1, 0] += ft[3]


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


def plot_mesh():
	"""
	Plot the initial mesh and print the mesh parameters.
	"""
	if model.plot_mesh == 'yes':
		for i in range(model.nbe):
			# plot Natural B.C. in red lines
			node1 = int(model.n_bc[0, i])
			node2 = int(model.n_bc[1, i])

			# coordinates
			x1 = model.x[node1 - 1]
			y1 = model.y[node1 - 1]
			x2 = model.x[node2 - 1]
			y2 = model.y[node2 - 1]

			plt.plot([x1, x2], [y1, y2], color='r', linewidth=4)

		for i in range(model.nel):
			XX = [model.x[model.IEN[0, i] - 1], model.x[model.IEN[1, i] - 1], model.x[model.IEN[2, i] - 1],
				  model.x[model.IEN[3, i] - 1], model.x[model.IEN[0, i] - 1]]
			YY = [model.y[model.IEN[0, i] - 1], model.y[model.IEN[1, i] - 1], model.y[model.IEN[2, i] - 1],
				  model.y[model.IEN[3, i] - 1], model.y[model.IEN[0, i] - 1]]
			plt.plot(XX, YY, color='b')

			if model.plot_nod == 'yes':
				plt.text(XX[0], YY[0], str(model.IEN[0, i]))
				plt.text(XX[1], YY[1], str(model.IEN[1, i]))
				plt.text(XX[2], YY[2], str(model.IEN[2, i]))
				plt.text(XX[3], YY[3], str(model.IEN[3, i]))

		plt.title('Meshed plate')
		plt.xlabel(r'$X$')
		plt.ylabel(r'$Y$')

	print('  Mesh Params ')
	print('No. of Elements  {}'.format(model.nel))
	print('No. of Nodes     {}'.format(model.nnp))
	print('No. of Equations {}'.format(model.neq))


def postprocess(ratio, wc):
	"""
	Plot the curve of w_c*D/q/L^4 vs. L/h.
	
	Args:
		ratio : L/h, ratio of length to thickness of the square plate.
		wc    : w_c*D/q/L^4, where w_c is the center deflection.
	"""
	if model.plot_mesh:
		# Convert matplotlib figures into PGFPlots figures
		if model.plot_tex == "yes":
			tikzplotlib.save("plate-mesh.pgf")
            
		plt.savefig("plate-mesh.pdf")
		plt.show()
	
	nh = len(ratio)
	s = np.ones(nh) * 0.00126
	
	# plot the curve of w_c*D/q/L^4 vs. L/h
	line1, = plt.plot(ratio, wc, color='b', label='Selective reduced integration')
	line2, = plt.plot(ratio, s, color='k', label='Exact thin plate solution')
	
	plt.xlim(model.r[0],model.r[1])
	plt.xticks([10,50,100,1000])
	
	plt.xlabel(r'$L/h$')
	plt.ylabel(r'$w_cD/qL^4$')
	
	plt.legend()
	plt.grid()
	plt.savefig("wc vs h.pdf")
	plt.show()
