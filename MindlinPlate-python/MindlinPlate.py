#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MindlinPlate - Elastic plate bending FE analysis program.
  Element types  : 4-node quadratic Mindlin plate element.
  Problem solved : plate whose Young's modulus, Poisson's rate and thickness are known
    is deformed at surface pressure, boundary forces, nodal forces and body force.

Usage:
   >>> MindlinPlate file_name
or
   $ python MindlinPlate.py file_name

Command line arguments:
  file_name: File name in which the FE model is stored in json format

Created on Fri August 15 18:56:57 2021

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

from sys import argv, exit
from PrePost import create_model_json, point_and_trac, postprocess
from MindlinPlateElem import MindlinPlateElem
import FEData as model
import numpy as np
from utitls import assembly, solvedr

def FERun(DataFile):
	# create FE model from DataFile in json format
	create_model_json(DataFile)

	ratio = np.arange(model.r[0], model.r[1], 2)
	nh = len(ratio)
	wc = np.zeros(nh)
	for index, ri in enumerate(ratio):
		# initialize K, d and f
		model.K = np.zeros((model.neq, model.neq))
		model.d = np.zeros((model.neq, 1))
		model.f = np.zeros((model.neq, 1))
		
		# the plate gets thinner as L/h increases
		model.h = model.lx / ri
		model.Db = np.array([[1, model.ne, 0],
							[model.ne, 1, 0],
							[0, 0, (1-model.ne)/2]])*model.E*model.h**3/(12.0*(1-model.ne**2))
		shcof = 5/6.0 #shear correction factor
		model.Ds = np.array([[1, 0],
							[0, 1]])*shcof*model.G*model.h

		# Calculation and assembly of element matrices
		for e in range(model.nel):
			ke, fe = MindlinPlateElem(e)
			assembly(e, ke, fe)

		# Compute and assemble nodal boundary force vector and point forces
		point_and_trac()

		print('\nL/h =', ri)
		# Solution phase
		f_E = solvedr()

		wc[index] = model.wc

	# Post process
	postprocess(ratio, wc)


if __name__ == "__main__":
	nargs = len(argv)
	if nargs == 2:
		DataFile = argv[1]
	else:
		print("Usage ： MindlinPlate file_name")
		exit()

	# DataFile = "./plate_64.json"
	FERun(DataFile)
