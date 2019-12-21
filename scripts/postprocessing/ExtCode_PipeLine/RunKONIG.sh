#!/bin/bash
#===============================================================================================================
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
#===============================================================================================================


# Ex. of Call: bash RunMeCvode.sh O3 10000 $WORKSPACE_PATH/neqplasma_QCT/cvode_code/install $WORKSPACE_PATH/Mars_Database/Run_0D/database/kinetics/ $WORKSPACE_PATH/Mars_Database/Run_0D/ 1 1 0 0

echo '------------------------------------------------------------------------------------------'
echo ' CoarseAIR: Coarse-Grained Quasi-Classical Trajectories 								    '
echo '------------------------------------------------------------------------------------------'
echo ' '
echo '------------------------------------------------------------------------------------------'
echo '   PipeLine for Running External Codes				 								    '
echo '------------------------------------------------------------------------------------------'
echo ' '

export System=${1}
export TTran=${2}
export PathToMECVODEFldr={3}
export PathToDtbFldr=${4}
export PathToRunFldr=${5}
export DissFlg=${6}
export InelFlg=${7}
export ExchFlg1=${8}
export ExchFlg2=${9}

ExtCode_SH_DIR=${COARSEAIR_SOURCE_DIR}"/scripts/postprocessing/ExtCode_PipeLine/"

echo '------------------------------------------------------'
echo '  Paths:'
echo '------------------------------------------------------'
echo '  ExtCode .sh     directory = '${ExtCode_SH_DIR}
echo '  MeCvode install directory = '${PathToMECVODEFldr}
echo '  MeCvode Dtb     directory = '${PathToDtbFldr}
echo '  MeCvode running directory = '${PathToRunFldr}
echo '------------------------------------------------------'
echo ' '

echo '------------------------------------------------------'
echo '  Inputs:'
echo '------------------------------------------------------'
echo '  System                         = '${System}
echo '  Translational Temp.            = '${TTran}
echo '  Writing Dissociation?          = '${DissFlg}
echo '  Writing Inelastic?             = '${InelFlg}
echo '  Firts Exchanges to be Written? = '${ExchFlg1}
echo '  Last Exchanges  to be Written? = '${ExchFlg2}
echo '------------------------------------------------------'
echo ' '


function Load_Initialize_0D{
  source $ExtCode_SH_DIR/Initialize_0D_Database_Function.sh
  Initialize_0D_Database_Function
}


function Call_MeCvode{
  cd ${PathToRunFldr}
  export OutputFldr='output_'${System}'_T'${TTran}'K_'${DissFlg}'_'${InelFlg}'_'${ExchFlg1}'_'${ExchFlg2}
  echo "[RunMECVODE]: MeCvode will be executed in the Folder "${OutputFldr}
  mkdir -p ./${OutputFldr}
  cd ./{OutputFldr}
  scp ${PathToMECVODEFldr}/${System}/'Mars_T'${TTran}'K/exec/run_' ./
  ./run_ ${OMP_NUM_THREADS}
}



echo "[RunMECVODE]: Calling Load_Initialize_0D"
Load_Initialize_0D
echo " "


echo "[RunMECCVODE]: Calling Call_MeCvode"
Call_MeCvode
echo " "