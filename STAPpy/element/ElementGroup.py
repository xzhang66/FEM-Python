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
from element.Bar import CBar
from element.Material import CBarMaterial

# dictionary: Define set of element types
ElementTypes = {0:'UNDEFINED',
				1:'Bar',
				2:'Q4',
				3:'T3',
				4:'H8',
				5:'Beam',
				6:'Plate',
				7:'Shell'}


class CElementGroup(object):
	""" Element group class """
	def __init__(self):
		from Domain import Domain
		FEMData = Domain()
		# List of all nodes in the domain, obtained from CDomain object
		self._NodeList = FEMData.GetNodeList()

		# Element type of this group
		self._ElementType = 0

		# Number of elements in this group
		self._NUME = 0

		# Element List in this group
		self._ElementList = []

		# Number of material/section property sets in this group
		self._NUMMAT = 0

		# Material list in this group
		self._MaterialList = []

	def __getitem__(self, item):
		""" operator [] """
		return self._ElementList[item]

	def GetMaterial(self, index):
		return self._MaterialList[index]

	def GetElementType(self):
		return self._ElementType

	def GetNUME(self):
		return self._NUME

	def GetNUMMAT(self):
		return self._NUMMAT

	def AllocateElements(self, amount):
		"""
		Allocate array of derived elements

		:param amount: (int) the amount of elements
		:return:
		"""
		element_type = ElementTypes.get(self._ElementType)
		if element_type == 'Bar':
			self._ElementList = [CBar() for _ in range(amount)]
		elif element_type == 'Q4':
			# implementation for other element types by yourself
			# ...
			pass # comment or delete this line after implementation
		else:
			error_info = "\nType {} not available. See CElementGroup." \
						 "AllocateElement.".format(self._ElementType)
			raise ValueError(error_info)

	def AllocateMaterials(self, amount):
		"""
		Allocate array of derived materials

		:param amount: (int) the amount of material sets
		:return: None
		"""
		element_type = ElementTypes.get(self._ElementType)
		if element_type == 'Bar':
			self._MaterialList = [CBarMaterial() for _ in range(amount)]
		elif element_type == 'Q4':
			# implementation for other element types by yourself
			# ...
			pass # comment or delete this line after implementation
		else:
			error_info = "\nType {} not available. See CElementGroup." \
						 "AllocateMaterials.".format(self._ElementType)
			raise ValueError(error_info)

	def Read(self, input_file):
		""" Read element group data from stream Input """
		line = input_file.readline().split()
		self._ElementType = int(line[0])
		self._NUME = int(line[1])
		self._NUMMAT = int(line[2])

		if not self.ReadElementData(input_file):
			return False

		return True

	def ReadElementData(self, input_file):
		""" Read bar element data from the input data file """
		# Read material/section property lines
		self.AllocateMaterials(self._NUMMAT)

		# Loop over for all material property sets in this element group
		for mset in range(self._NUMMAT):
			try:
				self.GetMaterial(mset).Read(input_file, mset)
			except ValueError as e:
				print(e)
				return False

		# Read element data lines
		self.AllocateElements(self._NUME)

		# Loop over for all elements in this element group
		for Ele in range(self._NUME):
			try:
				self[Ele].Read(input_file, Ele, self._MaterialList, self._NodeList)
			except ValueError as e:
				print(e)
				return False

		return True
