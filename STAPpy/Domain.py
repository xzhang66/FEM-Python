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
from utils.Outputter import COutputter
from element.Node import CNode
from LoadCaseData import CLoadCaseData
from element.ElementGroup import CElementGroup
from utils.SkylineMatrix import CSkylineMatrix
import numpy as np
import sys


@Singleton
class Domain(object):
	"""
	Domain class : Define the problem domain
	Only a single instance of Domain class can be created
	"""
	def __init__(self):
		super().__init__()

		# Input file stream for reading data from input data file
		self.input_file = None

		# Heading information for use in labeling the output
		self.Title = '0'

		# Solution MODEX
		# 		0 : Data check only
		# 		1 : Execution
		self.MODEX = 0

		# Total number of nodal points
		self.NUMNP = 0

		# List of all nodes in the domain
		self.NodeList = []

		# Total number of element groups
		self.NUMEG = 0

		# Element group list
		self.EleGrpList = []

		# Number of load cases
		self.NLCASE = 0

		# Number of concentrated loads applied in each load case
		self.NLOAD = []

		# List of all load cases
		self.LoadCases = []

		# Total number of equations in the system
		self.NEQ = 0

		# Global nodal force/displacement vector
		self.Force = None

		# Banded stiffness matrix
		# A one-dimensional array storing only the elements below the
		# skyline of the global stiffness matrix.
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

	def ReadData(self, input_filename, output_filename):
		""" Read domain data from the input data file """
		try:
			self.input_file = open(input_filename)
		except FileNotFoundError as e:
			print(e)
			sys.exit(3)

		Output = COutputter(output_filename)

		# Read the heading line
		self.Title = self.input_file.readline()
		Output.OutputHeading()

		# Read the control line
		line = self.input_file.readline().split()
		self.NUMNP = int(line[0])
		self.NUMEG = int(line[1])
		self.NLCASE = int(line[2])
		self.MODEX = int(line[3])

		# Read nodal point data
		if self.ReadNodalPoints():
			Output.OutputNodeInfo()
		else:
			return False

		# Update equation number
		self.CalculateEquationNumber()
		Output.OutputEquationNumber()

		# Read load data
		if self.ReadLoadCases():
			Output.OutputLoadInfo()
		else:
			return False

		# Read element data
		if self.ReadElements():
			Output.OutputElementInfo()
		else:
			return False

		return True

	def ReadNodalPoints(self):
		""" Read nodal point data """
		self.NodeList = [CNode() for _ in range(self.NUMNP)]

		for np in range(self.NUMNP):
			try:
				self.NodeList[np].Read(self.input_file, np)
			except ValueError as e:
				print(e)
				return False

		return True

	def CalculateEquationNumber(self):
		"""
		Calculate global equation numbers corresponding to every
		degree of freedom of each node
		"""
		self.NEQ = 0

		for np in range(self.NUMNP):
			for dof in range(CNode.NDF):
				if self.NodeList[np].bcode[dof]:
					self.NodeList[np].bcode[dof] = 0
				else:
					self.NEQ += 1
					self.NodeList[np].bcode[dof] = self.NEQ

	def ReadLoadCases(self):
		""" Read load case data """
		self.LoadCases = [CLoadCaseData() for _ in range(self.NLCASE)]

		for lcase in range(self.NLCASE):
			try:
				self.LoadCases[lcase].Read(self.input_file, lcase)
			except ValueError as e:
				print(e)
				return False

		return True

	def ReadElements(self):
		""" Read element data """
		self.EleGrpList = [CElementGroup() for _ in range(self.NUMEG)]

		for EleGrp in range(self.NUMEG):
			if not self.EleGrpList[EleGrp].Read(self.input_file):
				return False

		return True

	def CalculateColumnHeights(self):
		""" Calculate column heights """
		for EleGrp in range(self.NUMEG):
			ElementGrp = self.EleGrpList[EleGrp]
			NUME = ElementGrp.GetNUME()

			for Ele in range(NUME):
				Element = ElementGrp[Ele]

				Element.GenerateLocationMatrix()

				self.StiffnessMatrix.CalculateColumnHeight(
					Element.GetLocationMatrix(), Element.GetND())

		self.StiffnessMatrix.CalculateMaximumHalfBandwidth()

	def AssembleStiffnessMatrix(self):
		""" Assemble the banded gloabl stiffness matrix """
		# Loop over for all element groups
		for EleGrp in range(self.NUMEG):
			ElementGrp = self.EleGrpList[EleGrp]
			NUME = ElementGrp.GetNUME()
			size = ElementGrp[0].SizeOfStiffnessMatrix()
			Matrix = np.zeros(size, dtype=np.double)

			# Loop over for all elements in group EleGrp
			for Ele in range(NUME):
				Element = ElementGrp[Ele]
				Element.ElementStiffness(Matrix)
				self.StiffnessMatrix.Assembly(Matrix,
					Element.GetLocationMatrix(), Element.GetND())

			del Matrix

	def AssembleForce(self, LoadCase):
		""" Assemble the global nodal force vector for load case LoadCase """
		if LoadCase > self.NLCASE:
			return False

		LoadData = self.LoadCases[LoadCase - 1]

		# Loop over for all concentrated loads in load case LoadCase
		for lnum in range(LoadData.nloads):
			dof = self.NodeList[LoadData.node[lnum]-1].bcode[LoadData.dof[lnum]-1]

			if dof:
				self.Force[dof - 1] += LoadData.load[lnum]

		return True

	def AllocateMatrices(self):
		"""
		Allocate storage for matrices Force, ColumnHeights, DiagonalAddress
		and StiffnessMatrix and calculate the column heights and address
		of diagonal elements
		"""
		# Allocate for global force/displacement vector
		self.Force = np.zeros(self.NEQ, dtype=np.double)

		# Create the banded stiffness matrix
		self.StiffnessMatrix = CSkylineMatrix(self.NEQ)

		# Calculate column heights
		self.CalculateColumnHeights()

		# Calculate address of diagonal elements in banded matrix
		self.StiffnessMatrix.CalculateDiagnoalAddress()

		# Allocate for banded global stiffness matrix
		self.StiffnessMatrix.Allocate()

		Output = COutputter()
		Output.OutputTotalSystemData()
