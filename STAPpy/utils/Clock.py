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
import datetime


class Clock(object):
	""" Clock class for timing """
	def __init__(self):
		super().__init__()
		self._t0 = None
		self._t1 = None
		self._ct = datetime.timedelta(0)
		self._st0 = False					# Flag for Start method
		self._st1 = False					# Flag for Stop method

	def Start(self):
		""" Start the clock """
		self._t0 = datetime.datetime.now()
		self._st0 = True

	def Stop(self):
		""" Stop the clock """
		if not self._st0:
			error_info = "\n*** Error *** In Clock.Stop()" \
						 "\n : Method Start() must have been called before."
			raise RuntimeError(error_info)

		if not self._st1:
			self._t1 = datetime.datetime.now()
			self._ct += (self._t1 - self._t0)
			self._st1 = True

	def Resume(self):
		""" Resume the stopped clock """
		if not self._st0:
			error_info = "\n*** Error *** In Clock.Resume()" \
						 "\n : Method Start() must have been called before."
			raise RuntimeError(error_info)

		if not self._st1:
			error_info = "\n*** Error *** In Clock.Resume()" \
						 "\n : Method Stop() must have been called before."
			raise RuntimeError(error_info)
		else:
			self._t0 = datetime.datetime.now()
			self._st1 = False

	def Clear(self):
		""" Clear the clock """
		self._ct = datetime.timedelta(0)
		self._st0 = False
		self._st1 = False

	def ElapsedTime(self):
		""" Return the elapsed time since the clock started """
		if not self._st0:
			error_info = "\n*** Error *** In Clock.ElapsedTime()" \
						 "\n : Method Start() must have been called before."
			raise RuntimeError(error_info)

		if self._st1:
			elapsed = self._ct
		else:
			self._t1 = datetime.datetime.now()
			elapsed = self._ct + (self._t1 - self._t0)

		return elapsed
