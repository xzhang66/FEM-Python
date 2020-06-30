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
import abc


class CMaterial(metaclass=abc.ABCMeta):
	"""
		Material base class which only define one data member
		All type of material classes should be derived from this base class
	"""
	def __init__(self):
		self.nset = 0			# Number of set
		self.E = 0				# Young's modulus

	@abc.abstractmethod
	def Read(self, input_file, mset):
		pass

	@abc.abstractmethod
	def Write(self, output_file):
		pass


class CBarMaterial(CMaterial):
	""" Material class for bar element """
	def __init__(self):
		super().__init__()
		self.Area = 0			# Sectional area of a bar element

	def Read(self, input_file, mset):
		"""
		Read material data from stream Input
		"""
		line = input_file.readline().split()

		self.nset = np.int(line[0])
		if self.nset != mset + 1:
			error_info = "\n*** Error *** Material sets must be inputted in order !" \
						 "\n   Expected set : {}" \
						 "\n   Provided set : {}".format(mset + 1, self.nset)
			raise ValueError(error_info)

		self.E = np.double(line[1])
		self.Area = np.double(line[2])

	def Write(self, output_file):
		"""
		Write material data to Stream
		"""
		material_info = "%5d%16.6e%16.6e\n"%(self.nset, self.E, self.Area)

		# print the material info on the screen
		print(material_info, end='')
		# write the material info to output file
		output_file.write(material_info)