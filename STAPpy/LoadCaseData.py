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


class CLoadCaseData(object):
	""" Class LoadData is used to store load data """
	def __init__(self):
		self.nloads = 0    #!< Number of concentrated loads in this load case
		self.node = None   #!< Node number to which this load is applied
		self.dof = None    #!< Degree of freedom number for this load component
		self.load = None   #!< Magnitude of load

	def Allocate(self, num):
		self.nloads = num
		self.node = np.zeros(num, dtype=np.int)
		self.dof = np.zeros(num, dtype=np.int)
		self.load = np.zeros(num, dtype=np.double)

	def Read(self, input_file, lcase):
		"""
		Read load case data from stream Input

		:param input_file: (_io.TextIOWrapper) the object of input file
		:param lcase: check index
		:return: None
		"""
		line = input_file.readline().split()

		LL = int(line[0])
		NL = int(line[1])

		if LL != lcase + 1:
			error_info = "\n*** Error *** Load case must be inputted in order !" \
						 "\n   Expected load case : {}" \
						 "\n   Provided load case : {}".format(lcase + 1, LL)
			raise ValueError(error_info)

		self.Allocate(NL)

		for i in range(NL):
			line = input_file.readline().split()
			self.node[i] = np.int(line[0])
			self.dof[i] = np.int(line[1])
			self.load[i] = np.double(line[2])

	def Write(self, output_file, lcase):
		"""
		Write load case data to stream

		:param output_file: (_io.TextIOWrapper) the object of output file
		:param lcase: the index of load case
		:return: None
		"""
		for i in range(self.nloads):
			load_info = "%7d%13d%19.6e\n"%(self.node[i], self.dof[i],
										   self.load[i])
			print(load_info, end="")
			output_file.write(load_info)