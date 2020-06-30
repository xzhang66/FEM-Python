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
import abc


class CElement(metaclass=abc.ABCMeta):
	"""
	Element base class
	All type of element classes should be derived from this base class
	"""
	def __init__(self):
		# Number of nodes per element
		self._NEN = 0

		# Node list of the element
		self._nodes = []

		# Material of the element
		self._ElementMaterial = None

		# Location Matrix of the element
		self._LocationMatrix = None

		# Dimension of the location matrix
		self._ND = 0

	@abc.abstractmethod
	def Read(self, input_file, Ele, MaterialSets, NodeList):
		""" Read element data from stream Input """
		pass

	@abc.abstractmethod
	def Write(self, output_file, Ele):
		""" Write element data to stream """
		pass

	@abc.abstractmethod
	def GenerateLocationMatrix(self):
		"""
		Generate location matrix: the global equation number that
		corresponding to each DOF of the element
		"""
		pass

	@abc.abstractmethod
	def ElementStiffness(self, stiffness):
		"""
		Calculate element stiffness matrix
		(Upper triangular matrix, stored as an array column by colum)
		"""
		pass

	@abc.abstractmethod
	def ElementStress(self, stress, displacement):
		""" Calculate element stress """
		pass

	def GetNodes(self):
		""" Return nodes of the element """
		return self._nodes

	def GetElementMaterial(self):
		""" Return material of the element """
		return self._ElementMaterial

	def GetLocationMatrix(self):
		""" Return the Location Matrix of the element """
		return self._LocationMatrix

	def GetND(self):
		""" Return the dimension of the location matrix """
		return self._ND

	@abc.abstractmethod
	def SizeOfStiffnessMatrix(self):
		"""
		Return the size of the element stiffness matrix
		(stored as an array column by column)
		"""
		pass