#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exact solutions for some problems
	ExactSolution_Cantilever: Plot the exact solution of the cantilever given
	in Example 10.1 in Fish's textbook

Created on Aug. 15 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

import numpy as np

def ExactSolution_Fish_10_1(ax1, ax2, ax3):
	"""
	Plot the exact displacement, moment and shear force of the cantilever 
    beam (Example 10.1, Fish's book) in	ax1, ax2 and ax3, respectively.

	Args:
		ax1 : axis to draw displacement distribution
		ax2 : axis to draw moment distribution
		ax3 : axis to draw shear force distribution
	"""
	L = 12; a1 = 4; p1 = 10; a2 = 8; p2 = -5; a3 = 8; m = -20; p4 = 20
	p3 = 1; E = 1e4; I = 1.0

	x = np.arange(0, 12, 0.01)
	w = np.zeros(1200, np.float)
	M = np.zeros(1200, np.float)
	S = np.zeros(1200, np.float)

	for index, xi in enumerate(x):
		if xi < a1:
			w1 = -p1*xi**2/(6*E*I)*(3*a1-xi)
			M1 = p1 * a1 * (1 - xi / a1)
		else:
			w1 = -p1*a1**2/(6*E*I)*(3*xi-a1)
			M1 = 0

		if xi < a2:
			w2 = -p2*xi**2/(6*E*I)*(3*a2-xi)
			M2 = p2 * a2 * (1 - xi / a2)
		else:
			w2 = -p2*a2**2/(6*E*I)*(3*xi-a2)
			M2 = 0

		w3 = -m*(xi**2)/(2*E*I)
		M3 = m

		w4 = -p4*xi**2*(3*L-xi)/(6*E*I)
		M4 = p4*(L-xi)

		if xi < a3:
			w5 = -p3*xi**2*(6*a3**2-4*xi*a3+xi**2)/(24*E*I)
			M5 = 1/2*p3*(a3-xi)**2
		else:
			w5 = -p3*a3*(0.5*a3)**2/(6*E*I)*(3*xi-a3/2)
			M5 = 0

		if xi <= 12:
			S1 = 20
		else:
			S1 = 0

		if xi <= 8:
			S2 = -5
		else:
			S2 = 0

		if xi < 8:
			S3 = (8-xi)*1
		else:
			S3 = 0

		if xi < 4:
			S4 = 10
		else:
			S4 = 0

		w[index] = w1+w2+w3+w4+w5
		M[index] = -(M1+M2+M3+M4+M5)
		S[index] = S1+S2+S3+S4

	ax1.plot(x, w, '--r', label='Exact')
	ax2.plot(x, M, '--r', label='Exact')
	ax3.plot(x, S, '--r', label='Exact')
    

def ExactSolution_Ex_6_1(ax1, ax2, ax3):
	"""
	Plot the exact displacement, moment and shear force of the cantilever 
    beam (Example 6-1) in ax1, ax2 and ax3, respectively. 

	Args:
		ax1 : axis to draw displacement distribution
		ax2 : axis to draw moment distribution
		ax3 : axis to draw shear force distribution
	"""
    
	E = 1e4; I = 1.0

	x = np.arange(0, 8, 0.01)
	w = np.zeros(800, np.float)
	M = np.zeros(800, np.float)
	S = np.zeros(800, np.float)

	for index, xi in enumerate(x):
		if xi < 4:
			w[index] = (-xi**4/24 + 14*xi**3/3 - 71*xi**2)/(E*I)
			M[index] = -xi**2/2 + 28*xi - 142
			S[index] = 28 - xi
		else:
			w[index] = (-xi**4/24 + 3*xi**3 - 51*xi**2 - 80*xi + 320/3)/(E*I)
			M[index] = -xi**2/2 + 18*xi - 102
			S[index] = 18 - xi

	ax1.plot(x, w, '--r', label='Exact')
	ax2.plot(x, M, '--r', label='Exact')
	ax3.plot(x, S, '--r', label='Exact')    