#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv,exit

from PrePost import create_model_json, plotbeam

def FERun(DataFile, BeamType):
	# create FE model from DataFile in json format
	create_model_json(DataFile)

	# plot the models
	plotbeam()


if __name__ == "__main__":
	# nargs = len(argv)
	# if nargs == 2:
	# 	DataFile = argv[1]
	# 	BeamType = None
	# elif nargs == 3:
	# 	DataFile = argv[1]
	# 	BeamType = argv[2]
	# else:
	# 	print("Usage ： Beam1D file_name [BeamType]")
	# 	exit()
	DataFile = "beam_10_1.json"
	BeamType = None

	FERun(DataFile, BeamType)
