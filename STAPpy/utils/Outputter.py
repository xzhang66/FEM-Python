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
from utils.Singleton import Singleton
from element.ElementGroup import ElementTypes
import datetime
import numpy as np

weekday = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", 
           "Saturday", "Sunday"]
month = ["January", "February", "March", "April", "May", "June",
		 "July", "August", "September", "October", "November", "December"]


@Singleton
class COutputter(object):
	""" Singleton: Outputer class is used to output results """
	def __init__(self, filename=""):
		try:
			self._output_file = open(filename, 'w')
		except FileNotFoundError as e:
			print(e)
			sys.exit(3)

	def GetOutputFile(self):
		return self._output_file

	def PrintTime(self):
		""" Output current time and date """
		t = datetime.datetime.now()
		time_info = t.strftime("        (%H:%M:%S on ")
		time_info += (month[t.month - 1] + " ")
		time_info += (str(t.day) + ", ")
		time_info += (str(t.year) + ", ")
		time_info += (weekday[t.weekday()] + ")\n\n")

		print(time_info, end="")
		self._output_file.write(time_info)

	def OutputHeading(self):
		""" Print program logo """
		from Domain import Domain
		FEMData = Domain()

		title_info = "TITLE : " + FEMData.GetTitle() + "\n"
		print(title_info, end="")
		self._output_file.write(title_info)

		self.PrintTime()

	def OutputNodeInfo(self):
		""" Print nodal data """
		from Domain import Domain
		FEMData = Domain()

		NodeList = FEMData.GetNodeList()

		pre_info = "C O N T R O L   I N F O R M A T I O N\n\n"
		print(pre_info, end="")
		self._output_file.write(pre_info)

		NUMNP = FEMData.GetNUMNP()
		NUMEG = FEMData.GetNUMEG()
		NLCASE = FEMData.GetNLCASE()
		MODEX = FEMData.GetMODEX()

		pre_info = "\t  NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =%6d\n" \
				   "\t  NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =%6d\n" \
				   "\t  NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =%6d\n" \
				   "\t  SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =%6d\n" \
				   "\t\t EQ.0, DATA CHECK\n" \
				   "\t\t EQ.1, EXECUTION\n\n"%(NUMNP, NUMEG, NLCASE, MODEX)
		print(pre_info, end="")
		self._output_file.write(pre_info)

		pre_info = " N O D A L   P O I N T   D A T A\n\n" \
				   "    NODE       BOUNDARY                         NODAL POINT\n" \
				   "   NUMBER  CONDITION  CODES                     COORDINATES\n"
		print(pre_info, end="")
		self._output_file.write(pre_info)

		for n in range(NUMNP):
			NodeList[n].Write(self._output_file)

		print("\n", end="")
		self._output_file.write("\n")

	def OutputEquationNumber(self):
		""" Output equation numbers """
		from Domain import Domain
		FEMData = Domain()

		NodeList = FEMData.GetNodeList()

		NUMNP = FEMData.GetNUMNP()

		pre_info = " EQUATION NUMBERS\n\n" \
				   "   NODE NUMBER   DEGREES OF FREEDOM\n" \
				   "        N           X    Y    Z\n"
		print(pre_info, end="")
		self._output_file.write(pre_info)

		for n in range(NUMNP):
			NodeList[n].WriteEquationNo(self._output_file)

		print("\n", end="")
		self._output_file.write("\n")

	def OutputElementInfo(self):
		""" Output element data """
		# Print element group control line
		from Domain import Domain
		FEMData = Domain()

		NUMEG = FEMData.GetNUMEG()

		pre_info = " E L E M E N T   G R O U P   D A T A\n\n\n"
		print(pre_info, end="")
		self._output_file.write(pre_info)

		for EleGrp in range(NUMEG):
			ElementType = FEMData.GetEleGrpList()[EleGrp].GetElementType()
			NUME = FEMData.GetEleGrpList()[EleGrp].GetNUME()

			pre_info = " E L E M E N T   D E F I N I T I O N\n\n" \
					   " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =%5d\n" \
					   "     EQ.1, TRUSS ELEMENTS\n" \
					   "     EQ.2, ELEMENTS CURRENTLY\n" \
					   "     EQ.3, NOT AVAILABLE\n\n" \
					   " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =%5d\n\n" \
					   %(ElementType, NUME)
			print(pre_info, end="")
			self._output_file.write(pre_info)

			element_type = ElementTypes.get(ElementType)
			if element_type == 'Bar':
				self.PrintBarElementData(EleGrp)
			elif element_type == 'Q4':
				# implementation for other element types by yourself
				# ...
				pass  # comment or delete this line after implementation
			else:
				error_info = "\n*** Error *** Elment type {} has not been " \
							 "implemented.\n\n".format(ElementType)
				raise ValueError(error_info)

	def PrintBarElementData(self, EleGrp):
		""" Output bar element data """
		from Domain import Domain
		FEMData = Domain()

		ElementGroup = FEMData.GetEleGrpList()[EleGrp]
		NUMMAT = ElementGroup.GetNUMMAT()

		pre_info = " M A T E R I A L   D E F I N I T I O N\n\n" \
				   " NUMBER OF DIFFERENT SETS OF MATERIAL\n" \
				   " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =%5d\n\n" \
				   "  SET       YOUNG'S     CROSS-SECTIONAL\n" \
				   " NUMBER     MODULUS          AREA\n" \
				   "               E              A\n"%NUMMAT
		print(pre_info, end="")
		self._output_file.write(pre_info)

		for mset in range(NUMMAT):
			ElementGroup.GetMaterial(mset).Write(self._output_file)

		pre_info = "\n\n E L E M E N T   I N F O R M A T I O N\n" \
				   " ELEMENT     NODE     NODE       MATERIAL\n" \
				   " NUMBER-N      I        J       SET NUMBER\n"
		print(pre_info, end="")
		self._output_file.write(pre_info)

		NUME = ElementGroup.GetNUME()
		for Ele in range(NUME):
			ElementGroup[Ele].Write(self._output_file, Ele)

		print("\n", end="")
		self._output_file.write("\n")

	def OutputLoadInfo(self):
		""" Print load data """
		from Domain import Domain
		FEMData = Domain()

		for lcase in range(FEMData.GetNLCASE()):
			LoadData = FEMData.GetLoadCases()[lcase]

			pre_info = " L O A D   C A S E   D A T A\n\n" \
					   "     LOAD CASE NUMBER . . . . . . . =%6d\n" \
					   "     NUMBER OF CONCENTRATED LOADS . =%6d\n\n" \
					   "    NODE       DIRECTION      LOAD\n" \
					   "   NUMBER                   MAGNITUDE\n"%(lcase + 1,
																  LoadData.nloads)
			print(pre_info, end="")
			self._output_file.write(pre_info)

			LoadData.Write(self._output_file, lcase+1)

			print("\n", end="")
			self._output_file.write("\n")

	def OutputNodalDisplacement(self, lcase):
		""" Print nodal displacement """
		from Domain import Domain
		FEMData = Domain()
		NodeList = FEMData.GetNodeList()
		displacement = FEMData.GetDisplacement()

		pre_info = " LOAD CASE%5d\n\n\n" \
				   " D I S P L A C E M E N T S\n\n" \
				   "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT\n" \
				   %(lcase+1)
		print(pre_info, end="")
		self._output_file.write(pre_info)

		for n in range(FEMData.GetNUMNP()):
			NodeList[n].WriteNodalDisplacement(self._output_file, displacement)

		print("\n", end="")
		self._output_file.write("\n")

	def OutputElementStress(self):
		""" Calculate stresses """
		from Domain import Domain
		FEMData = Domain()

		displacement = FEMData.GetDisplacement()

		NUMEG = FEMData.GetNUMEG()

		for ELeGrpIndex in range(NUMEG):
			pre_info = " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P%5d\n\n" \
					   %(ELeGrpIndex+1)
			print(pre_info, end="")
			self._output_file.write(pre_info)

			EleGrp = FEMData.GetEleGrpList()[ELeGrpIndex]
			NUME = EleGrp.GetNUME()
			ElementType = EleGrp.GetElementType()

			element_type = ElementTypes.get(ElementType)
			if element_type == 'Bar':
				pre_info = "  ELEMENT             FORCE            STRESS\n" \
						   "  NUMBER\n"
				print(pre_info, end="")
				self._output_file.write(pre_info)

				stress = np.zeros(1)

				for Ele in range(NUME):
					Element = EleGrp[Ele]
					Element.ElementStress(stress, displacement)

					material = Element.GetElementMaterial()
					stress_info = "%5d%22.6e%18.6e\n"%(Ele+1, stress[0]*material.Area, stress[0])
					print(stress_info, end="")
					self._output_file.write(stress_info)
			elif element_type == 'Q4':
				# implementation for other element types by yourself
				# ...
				pass  # comment or delete this line after implementation
			else:
				error_info = "\n*** Error *** Elment type {} has not been " \
							 "implemented.\n\n".format(ElementType)
				raise ValueError(error_info)

	def OutputTotalSystemData(self):
		""" Print total system data """
		from Domain import Domain
		FEMData = Domain()

		pre_info = "	TOTAL SYSTEM DATA\n\n" \
				   "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = {}\n" \
				   "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = {}\n" \
				   "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = {}\n" \
				   "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = {}\n\n\n".format(
			FEMData.GetNEQ(), FEMData.GetStiffnessMatrix().size(),
			FEMData.GetStiffnessMatrix().GetMaximumHalfBandwidth(),
			FEMData.GetStiffnessMatrix().size()/FEMData.GetNEQ()
		)
		print(pre_info, end="")
		self._output_file.write(pre_info)

	def OutputSolutionTime(self, time_info):
		""" Print CPU time used for solution """
		print(time_info, end="")
		self._output_file.write(time_info)
