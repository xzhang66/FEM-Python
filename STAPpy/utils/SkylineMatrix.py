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
import numpy as np
import sys


class CSkylineMatrix(object):
	"""
	CSkylineMatrix class is used to store the FEM stiffness matrix
	in skyline storage
	"""
	def __init__(self, N):
		super().__init__()
        
		# Dimension of the stiffness matrix
		self._NEQ = N

		# Maximum half bandwith
		self._MK = 0

		# Size of the storage used to store the stiffness matrkix in skyline
		self._NWK = 0

		# Store the stiffness matrkix in skyline storage
		self._data = None

		# Column hights
		self._ColumnHeights = np.zeros(N, dtype=np.int)

		# Diagonal address of all columns in data_
		self._DiagonalAddress = np.zeros(N+1, dtype=np.int)

	def Index(self, i, j):
		""" Return the index in self._data of (i, j) in K """
		if j >= i:
			return self._DiagonalAddress[j - 1] + (j - i) - 1
		else:
			return self._DiagonalAddress[i - 1] + (i - j) - 1

	def __getitem__(self, *item):
		(i, j), =item
		index = self.Index(i, j)
		return self._data[index]

	def __setitem__(self, *item):
		(i, j), value = item
		index = self.Index(i, j)
		self._data[index] = value

	def Allocate(self):
		""" Allocate storage for the matrix """
		self._NWK = self._DiagonalAddress[self._NEQ] - self._DiagonalAddress[0]
		self._data = np.zeros(self._NWK, dtype=np.double)

	def GetColumnHeights(self):
		""" Return pointer to the _ColumnHeights """
		return self._ColumnHeights

	def GetMaximumHalfBandwidth(self):
		""" Return the maximum half bandwidth """
		return self._MK

	def GetDiagonalAddress(self):
		""" Return pointer to the _DiagonalAddress """
		return self._DiagonalAddress

	def dim(self):
		""" Return the dimension of the stiffness matrix """
		return self._NEQ

	def size(self):
		"""
		Return the size of the storage used to store the stiffness matrkix
		in skyline
		"""
		return self._NWK

	def CalculateColumnHeight(self, LocationMatrix, ND):
		"""
		Calculate the column height, used with the skyline storage scheme
		"""
		# Look for the row number of the first non-zero element
		nfisrtrow = sys.maxsize
		for i in range(ND):
			if LocationMatrix[i] and LocationMatrix[i] < nfisrtrow:
				nfisrtrow = LocationMatrix[i]

		# Calculate the column height contributed by this element
		for i in range(ND):
			column = LocationMatrix[i]
			if not column:
				continue

			Height = column - nfisrtrow
			if self._ColumnHeights[column - 1] < Height:
				self._ColumnHeights[column - 1] = Height

	def CalculateMaximumHalfBandwidth(self):
		""" Maximum half bandwidth ( = max(ColumnHeights) + 1 ) """
		self._MK = self._ColumnHeights.max() + 1

	def Assembly(self, Matrix, LocationMatrix, ND):
		"""
		Assemble the banded global stiffness matrix (skyline storage scheme)
		"""
		# Assemble global stiffness matrix
		for j in range(ND):
			# Global equation number corresponding to jth DOF of the element
			Lj = LocationMatrix[j]
			if not Lj:
				continue

			# Address of diagonal element of column j
			# in the one dimensional element stiffness matrix
			DiagjElement = int((j + 1)*j/2 + 1)

			for i in range(j + 1):
				# Global equation number corresponding to ith DOF of the element
				Li = LocationMatrix[i]
				if not Li:
					continue

				self._data[self.Index(Li, Lj)] += Matrix[DiagjElement + j - i - 1]

	def CalculateDiagnoalAddress(self):
		"""
		Calculate address of diagonal elements in banded matrix
		Caution: Address is numbered from 1
		M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)

		:return: None
		"""
		self._DiagonalAddress[0] = 1
		for col in range(1, self._NEQ+1):
			self._DiagonalAddress[col] = self._DiagonalAddress[col - 1] \
										 + self._ColumnHeights[col - 1] + 1
