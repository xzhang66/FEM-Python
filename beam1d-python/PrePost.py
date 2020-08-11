#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provide methods to create FE model for a beam from a json file, to plot the
beam, to print stresses at Gauss points, to plot displacement and stress
distributions obtained by FE analysis and exact solution.

Created on Aug. 11 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

import numpy as np
import json
import matplotlib.pyplot as plt
import tikzplotlib

import FEData as model
from Beam1DElem import setup_ID_LM


def create_model_json(DataFile):
	"""
	Initialize the FEM model from file DataFile (in json format)
	"""

	# input data from json file
	with open(DataFile) as f_obj:
		FEData = json.load(f_obj)

	model.Title = FEData['Title']
	model.nsd = FEData['nsd']
	model.ndof = FEData['ndof']
	model.nnp = FEData['nnp']
	model.nel = FEData['nel']
	model.nen = FEData['nen']
	model.neq = model.ndof*model.nnp
	model.neqe = model.ndof*model.nen

	# initialize K, d and f
	model.f = np.zeros((model.neq, 1))
	model.d = np.zeros((model.neq, 1))
	model.K = np.zeros((model.neq, model.neq))

	# element and material data (given at the element nodes)
	model.E = np.array(FEData['E'])
	model.body = np.array(FEData['body'])
	model.CArea = np.array(FEData['CArea'])
	model.leng = np.array(FEData['leng'])

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
		model.P = np.array(FEData['P'])

	# output plots
	model.plot_beam = FEData['plot_beam']
	model.plot_nod = FEData['plot_nod']
	model.nplot = model.nen * 10
	model.plot_tex = FEData['plot_tex']

	# define the mesh
	model.x = np.array(FEData['x'])
	model.y = np.array(FEData['y'])
	model.IEN = np.array(FEData['IEN'], np.int)

	model.ID = np.zeros(model.neq, np.int)
	model.LM = np.zeros((model.nen, model.nel), np.int)

	# generate LM and ID arrays
	setup_ID_LM()

def plotbeam():
	""" Plot the beam """
	if model.plot_beam == 'yes':
		for i in range(model.nel):
			XX = np.array([model.x[model.IEN[0,i]-1], model.x[model.IEN[1,i]-1]])
			YY = np.array([model.y[model.IEN[0,i]-1], model.y[model.IEN[1,i]-1]])
			plt.figure(1)
			plt.plot(XX, YY)
			plt.plot(XX, -YY)
			plt.plot(XX, [0,0], '+r')

			if model.plot_nod == 'yes':
				plt.text(XX[0], 0, str(model.IEN[0, i]))
				plt.text(XX[1], 0, str(model.IEN[1, i]))

		plt.plot([ model.x[model.IEN[0,0] - 1], model.x[model.IEN[0,0] - 1]],
				 [-model.y[model.IEN[0,0] - 1], model.y[model.IEN[0,0] - 1]])
		plt.plot([ model.x[model.IEN[-1,-1] - 1], model.x[model.IEN[-1,-1] - 1]],
				 [-model.y[model.IEN[-1,-1] - 1], model.y[model.IEN[-1,-1] - 1]])
		plt.title('Beam Plot')
		plt.axis('equal')
		plt.show()

		# print some mesh parameters
		print('\n  Beam Params')
		print('No. of Elements  ' + str(model.nel))
		print('No. of Nodes     ' + str(model.nnp))
		print('No. of Equations ' + str(model.neq))