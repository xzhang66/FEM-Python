#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
/*****************************************************************************/
/*  STAPpy : A python FEM code sharing the same input data file with STAP90  */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Created on Mon Jun 22, 2020                                           */
/*                                                                           */
/*     @author: thurcni@163.com, xzhang@tsinghua.edu.cn                      */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/
"""
import sys
sys.path.append('../')
from solver.Solver import CSolver
import numpy as np
import sys


class CLDLTSolver(CSolver):
	""" LDLT solver: A in core solver using skyline storage  and column reduction scheme """
	def __init__(self, K):
		self.K = K			# Global Stiffness matrix in Skyline storage

	def LDLT(self):
		""" LDLT facterization """
		N = self.K.dim()
		ColumnHeights = self.K.GetColumnHeights()

		for j in range(2, N+1): # Loop for column 2:n (Numbering starting from 1)
			# Row number of the first non-zero element in column j (Numbering starting from 1)
			mj = j - ColumnHeights[j - 1]

			for i in range(mj+1, j): # Loop for mj+1:j-1 (Numbering starting from 1)
				# Row number of the first nonzero element in column i (Numbering starting from 1)
				mi = i - ColumnHeights[i - 1]

				C = np.double(0.0)
				for r in range(max(mi, mj), i):
					C += (self.K[self.K.Index(r, i)]*self.K[self.K.Index(r, j)]) # C += L_ri * U_rj

				self.K[self.K.Index(i, j)] -= C		# U_ij = K_ij - C

			for r in range(mj, j):		# Loop for mj:j-1 (column j)
				Lrj = self.K[self.K.Index(r, j)]/self.K[self.K.Index(r, r)] # L_rj = U_rj / D_rr
				self.K[self.K.Index(j, j)] -= (Lrj*self.K[self.K.Index(r, j)]) # D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
				self.K[self.K.Index(r, j)] = Lrj

			if np.fabs(self.K[self.K.Index(j, j)] <= sys.float_info.min):
				error_info = "\n*** Error *** Stiffness matrix is not positive definite !" \
							 "\n    Euqation no = {}" \
							 "\n    Pivot = {}".format(j, self.K[self.K.Index(j, j)])
				raise ValueError(error_info)

	def BackSubstitution(self, Force):
		""" Solve displacement by back substitution """
		N = self.K.dim()
		ColumnHeights = self.K.GetColumnHeights()

		# Reduce right-hand-side load vector (LV = R)
		for i in range(2, N+1): # Loop for i=2:N (Numering starting from 1)
			mi = i - ColumnHeights[i - 1]

			for j in range(mi, i): # Loop for j=mi:i-1
				Force[i - 1] -= (self.K[self.K.Index(j, i)]*Force[j - 1]) # V_i = R_i - sum_j (L_ji V_j)

		# Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
		for i in range(1, N+1): # Loop for i=1:N
			Force[i - 1] /= self.K[self.K.Index(i, i)] # Vbar = D^(-1) V

		for j in range(N, 1, -1): # Loop for j=N:2
			mj = j - ColumnHeights[j - 1]

			for i in range(mj, j): # Loop for i=mj:j-1
				Force[i - 1] -= (self.K[self.K.Index(i, j)]*Force[j - 1]) # a_i = Vbar_i - sum_j(L_ij Vbar_j)
