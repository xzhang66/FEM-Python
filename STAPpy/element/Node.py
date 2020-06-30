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


class CNode(object):
	# Maximum number of degrees of freedom per node
	# For 3D bar and solid elements, NDF = 3.
	# For 3D beam or shell elements, NDF = 5 or 6
	NDF = 3

	def __init__(self, x=0.0, y=0.0, z=0.0):
		super().__init__()
		# x, y and z coordinates of the node
		self.XYZ = np.zeros(CNode.NDF)
		self.XYZ[0] = x; self.XYZ[1] = y; self.XYZ[2] = z

		# Boundary code of each degree of freedom of the node
		#     0: The corresponding degree of freedom is active
		#     		(defined in the global system)
		#     1: The corresponding degree of freedom in nonactive
		#     		(not defined)
		# After call Domain.CalculateEquationNumber(),
		# bcode stores the global equation number
		# corresponding to each degree of freedom of the node
		self.bcode = np.zeros(CNode.NDF, dtype=np.int)

		# Node numer
		self.NodeNumber = 0

	def Read(self, input_file, check_np):
		"""
		Read element data from stream Input
		"""
		line = input_file.readline().split()

		N = int(line[0])
		if N != check_np + 1:
			error_info = "\n*** Error *** Nodes must be inputted in order !" \
						 "\n   Expected node number : {}" \
						 "\n   Provided node number : {}".format(check_np+1, N)
			raise ValueError(error_info)

		self.NodeNumber = N

		self.bcode[0] = np.int(line[1])
		self.bcode[1] = np.int(line[2])
		self.bcode[2] = np.int(line[3])
		self.XYZ[0] = np.double(line[4])
		self.XYZ[1] = np.double(line[5])
		self.XYZ[2] = np.double(line[6])

	def Write(self, output_file):
		"""
		Output nodal point data to stream
		"""
		node_info = "%9d%5d%5d%5d%18.6e%15.6e%15.6e\n"%(
			self.NodeNumber, self.bcode[0], self.bcode[1], self.bcode[2],
			self.XYZ[0], self.XYZ[1], self.XYZ[2])
		# print the nodal info on the screen
		print(node_info, end='')
		# write the nodal info to output file
		output_file.write(node_info)

	def WriteEquationNo(self, output_file):
		"""
		Output equation numbers of nodal point to stream
		"""
		equation_info = "%9d       "%self.NodeNumber

		for dof in range(CNode.NDF):
			equation_info += "%5d"%self.bcode[dof]

		equation_info += '\n'
		# print the nodal info on the screen
		print(equation_info, end='')
		# write the nodal info to output file
		output_file.write(equation_info)

	def WriteNodalDisplacement(self, output_file, displacement):
		"""
		Write nodal displacement
		"""
		displacement_info = "%5d        "%self.NodeNumber

		for dof in range(CNode.NDF):
			if self.bcode[dof] == 0:
				displacement_info += "%18.6e"%0.0
			else:
				displacement_info += "%18.6e"%displacement[self.bcode[dof] - 1]

		displacement_info += '\n'
		# print the nodal info on the screen
		print(displacement_info, end='')
		# write the nodal info to output file
		output_file.write(displacement_info)