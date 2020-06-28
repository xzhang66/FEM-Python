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


class Singleton(object):
	_Instance = {}

	def __init__(self, cls):
		self.cls = cls

	def __call__(self, *args, **kwargs):
		instance = self._Instance.get(self.cls, None)
		if not instance:
			instance = self.cls(*args, **kwargs)
			self._Instance[self.cls] = instance
		return instance

	def __getattr__(self, item):
		return getattr(self.cls, item, None)
