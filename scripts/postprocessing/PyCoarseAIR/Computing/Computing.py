##==============================================================================================================
# 
# Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
# 
# Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
#
# Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
# 
# This program is free software; you can redistribute it and/or modify it under the terms of the 
# Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU Lesser General Public License for more details. 
# 
# You should have received a copy of the GNU Lesser General Public License along with this library; 
# if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
# 
#---------------------------------------------------------------------------------------------------------------
##==============================================================================================================
import numpy as np
import pandas
import csv
import sys

from matplotlib import rc 
import matplotlib.pyplot as plt

sys.path.insert(0, '../Plotting/')

from Plotting      import Plot_DissRates_Thermal


def Compute_Correction_To_DissRates(InputData, Syst, iT):

    Syst.T[iT-1].Proc[0].Rates[:,0] = Syst.T[iT-1].Proc[0].Rates[:,0] * InputData.Kin.CorrFactor

    return Syst



def Compute_Rates_Overall(Syst, iT):

    Syst.T[iT-1].ProcTot[0].Rates = np.sum(Syst.T[iT-1].Proc[0].Rates, axis=1)
    Syst.T[iT-1].ProcTot[1].Rates = np.sum(Syst.T[iT-1].Proc[1].Rates, axis=1)
    for jProc in range(2, Syst.NProcTypes):
        Syst.T[iT-1].ProcTot[jProc].Rates = np.sum(Syst.T[iT-1].ProcExch[jProc-2].Rates, axis=1)

    return Syst



def Compute_Rates_Thermal(Syst, iT):

    for jProc in range(Syst.NProcTypes):
        Syst.RatesTh[iT-1,jProc] = sum( Syst.Molecule[0].T[iT-1].LevelQExp * Syst.T[iT-1].ProcTot[jProc].Rates )

    return Syst



def Compute_Rates_Thermal_FromOverall(Syst, Temp, InputData):

    DissFile = InputData.Kin.ReadFldr   + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name + '.csv'
    print('  [Compute_Rates_Thermal_FromOverall]: Reading Dissociation Rates From File: ' + DissFile)
    for iT in Temp.iTVec:
        LevelKDiss = np.zeros(Syst.Molecule[0].NBins)
        with open(DissFile) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')
            next(readCSV)
            for row in readCSV:
                if (float(row[iT]) > 0.0):
                    LevelKDiss[int(row[0])-1] = float(row[iT])
        csvfile.close()

        Syst.RatesTh[iT-1,0] = sum( Syst.Molecule[0].T[iT-1].LevelQExp * LevelKDiss )
    
    print('  [Compute_Rates_Thermal_FromOverall]: Thermal Dissociation Rates = ', Syst.T[iT-1].ProcTot[0].Rates)
    Write_DissRates_Thermal(Syst, Temp, InputData)

    print('  [Compute_Rates_Thermal_FromOverall]: Plotting Dissociation Thermal Rates')
    Plot_DissRates_Thermal(Syst, Temp, InputData)

    return Syst



def Compute_QSS(Syst, ME, iT):

    KQSS_Eps   = 1.e-7

    KDer = np.gradient(Syst.T[iT-1].ProcTot[0].RatesAveraged, ME.Time)
    iQSS = 0
    while (KDer[iQSS] > KQSS_Eps):
        iQSS = iQSS + 1
    iQSS_Start = iQSS - 1
    while (KDer[iQSS] < KQSS_Eps):
        iQSS = iQSS + 1
    iQSS_End = iQSS - 1
    Syst.T[iT-1].QSS.iTime[0] = iQSS_Start
    Syst.T[iT-1].QSS.iTime[1] = iQSS_End
    Syst.T[iT-1].QSS.Time[0]  = ME.Time[iQSS_Start]
    Syst.T[iT-1].QSS.Time[1]  = ME.Time[iQSS_End]
    
    Syst.T[iT-1].QSS.Rate[0]  = ( Syst.T[iT-1].ProcTot[0].RatesAveraged[iQSS_Start] + Syst.T[iT-1].ProcTot[0].RatesAveraged[iQSS_End] ) / 2.0
    for jProc in range(2, Syst.NProcTypes):
        Syst.T[iT-1].QSS.Rate[jProc] = ( Syst.T[iT-1].ProcTot[jProc].RatesAveraged[iQSS_Start] + Syst.T[iT-1].ProcTot[jProc].RatesAveraged[iQSS_End] ) / 2.0

    return Syst
