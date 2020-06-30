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


class CSolver(metaclass=abc.ABCMeta):
	"""
	Interface class for all the solver. New solver based on LDLT factorization
	should be derived from this base class.
	This class defines the common methods in these solvers.
	"""
	@abc.abstractmethod
	def LDLT(self):
		""" Perform L*D*L(T) factorization of the stiffness matrix """
		pass

	@abc.abstractmethod
	def BackSubstitution(self, Force):
		""" Reduce right-hand-side load vector and back substitute """
		pass