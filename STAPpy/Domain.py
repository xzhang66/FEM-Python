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
from utils.Singleton import Singleton


@Singleton
class Domain(object):
	"""
	Domain class : Define the problem domain
	Only a single instance of Domain class can be created
	"""
	def __init__(self):
		super().__init__()
		# Heading information for use in labeling the output
		self.Title = '0'

		# Solution MODEX
		# 		0 : Data check only
		# 		1 : Execution
		self.MODEX = 0

		# Total number of nodal points
		self.NUMNP = 0

		# List of all nodes in the domain
		self.NodeList = None

		# Total number of element groups
		self.NUMEG = 0

		# Element group list
		self.EleGrpList = None

		# Number of load cases
		self.NLCASE = 0

		# Number of concentrated loads applied in each load case
		self.NLOAD = None

		# List of all load cases
		self.LoadCases = None

		# Total number of equations in the system
		self.NEQ = 0

		# Global nodal force/displacement vector
		self.Force = None

		# Banded stiffness matrix
		# A one-dimensional array storing only the elements below the	skyline of the
		# global stiffness matrix.
		self.StiffnessMatrix = None

	def GetMODEX(self):
		return self.MODEX

	def GetTitle(self):
		return self.Title

	def GetNEQ(self):
		return self.NEQ

	def GetNUMNP(self):
		return self.NUMNP

	def GetNodeList(self):
		return self.NodeList

	def GetNUMEG(self):
		return self.NUMEG

	def GetEleGrpList(self):
		return self.EleGrpList

	def GetForce(self):
		return self.Force

	def GetDisplacement(self):
		return self.Force

	def GetNLCASE(self):
		return self.NLCASE

	def GetNLOAD(self):
		return self.NLOAD

	def GetLoadCases(self):
		return self.LoadCases

	def GetStiffnessMatrix(self):
		return self.StiffnessMatrix