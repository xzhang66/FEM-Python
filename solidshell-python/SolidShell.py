#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SolidShell - 3-dimension elasticity FE analysis program.
  Element types  : 3-dimension 8-node solid shell element.
  Problem solved : 3-dimension plane whose Young's modulus and Poisson's rate are known
    is deformed at nodal forces.

Usage:
   >>> SolidShell file_name
or
   $ python SolidShell.py file_name

Command line arguments:
  file_name: File name in which the FE model is stored in json format

Created on Fri Apr 19 18:56:57 2024

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""
from sys import argv, exit
from PrePost import create_model_json, point_and_trac, postprocess
from SolidShellElem import SolidShellElem
import FEData as model
from utitls import assembly, solvedr


def FERun(DataFile):
	# create FE model from DataFile in json format
	create_model_json(DataFile)

	# Calculation and assembly of element matrices
	for e in range(model.nel):
		ke, fe = SolidShellElem(e)
		assembly(e, ke, fe)

	# Compute and assemble nodal boundary force vector and point forces
	point_and_trac()

	# Solution phase
	f_E = solvedr()

	# Post process
	postprocess()


if __name__ == "__main__":
	nargs = len(argv)
	if nargs == 2:
		DataFile = argv[1]
	else:
		print("Usage ： SolidShell file_name")
		exit()

	# DataFile = "./elasticity_16.json"
	FERun(DataFile)
	
