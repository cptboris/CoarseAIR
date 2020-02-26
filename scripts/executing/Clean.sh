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

# ----------------------------------------------------------------------------------------------------------------- Clean #
function RmAll {
  
  export COARSEAIR_OUTPUT_DIR=$(pwd)/Test
  export TranFlg=${1}
  export Tran=${2}
  export Tint=${3}
  export NMolecules=${4}
  export MinLevel1=${5}
  export MaxLevel1=${6}
  export MinLevel1=${7}
  export MaxLevel1=${8}

  echo "  [RmAll]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
  echo "  [RmAll]: TranFlg               = "${TranFlg}
  echo "  [RmAll]: Tran                  = "${Tran}
  echo "  [RmAll]: Tint                  = "${Tint}
  echo "  [RmAll]: NMolecules            = "${NMolecules}
  echo "  [RmAll]: MinLevel1             = "${MinLevel1}
  echo "  [RmAll]: MaxLevel1             = "${MaxLevel1}
  echo "  [RmAll]: MinLevel2             = "${MinLevel2}
  echo "  [RmAll]: MaxLevel2             = "${MaxLevel2}
   

  MinLevel1Temp=${MinLevel1}
  MaxLevel1Temp=${MaxLevel1}
  MinLevel2Temp=${MinLevel2}
  MaxLevel2Temp=${MaxLevel2}
  for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
    if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
      if [ ${MinLevel2} -eq 0 ] && [ ${NMolecules} -gt 1 ]; then
        MinLevel2Temp=${iLevel1}
      fi
      for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
        if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
          echo "  [RmAll]: --- Molecule 1, Level/Bin " ${iLevels1} " ----------------------------- "
          echo "  [RmAll]: ----- Molecule 2, Level/Bin = " ${iLevels2} " --------------------- "
          echo "  [RmAll]"

      
          if [ ${TranFlg} -eq 0 ]; then 
            export COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
          else
            export COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
          fi


          #if [ -f ${COARSEAIR_BIN_OUTPUT_DIR} ]; then
          #  echo ${COARSEAIR_BIN_OUTPUT_DIR}
          rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}
          #fi
        fi
        echo "  [RmAll]: ---------------------------------------------------------- "
      done
    fi
    echo "  [RmAll]: ------------------------------------------------------------ "
    echo " "
  done         
      
}
#================================================================================================================================#


# ------------------------------------------------------------------------------------------------------------------------ Clean #
function Clean {
  # Ex: Clean 1 10000.0 10000.0 1 9390 0 0

  export COARSEAIR_OUTPUT_DIR=$(pwd)/Test
  export TranFlg=${1}
  export Tran=${2}
  export Tint=${3}
  export NMolecules=${4}
  export MinLevel1=${5}
  export MaxLevel1=${6}
  export MinLevel1=${7}
  export MaxLevel1=${8}
  
  echo "  [Clean]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
  echo "  [Clean]: TranFlg               = "${TranFlg}
  echo "  [Clean]: Tran                  = "${Tran}
  echo "  [Clean]: Tint                  = "${Tint}
  echo "  [Clean]: NMolecules            = "${NMolecules}
  echo "  [Clean]: MinLevel1             = "${MinLevel1}
  echo "  [Clean]: MaxLevel1             = "${MaxLevel1}
  echo "  [Clean]: MinLevel2             = "${MinLevel2}
  echo "  [Clean]: MaxLevel2             = "${MaxLevel2}


  MinLevel1Temp=${MinLevel1}
  MaxLevel1Temp=${MaxLevel1}
  MinLevel2Temp=${MinLevel2}
  MaxLevel2Temp=${MaxLevel2}
  for (( iLevel1=1; iLevel1<=${NLevels1}; iLevel1++ )); do
    if [ ${iLevel1} -ge ${MinLevel1Temp} ] && [ ${iLevel1} -le ${MaxLevel1Temp} ]; then
      if [ ${MinLevel2} -eq 0 ] && [ ${NMolecules} -gt 1 ]; then
        MinLevel2Temp=${iLevel1}
      fi
      for (( iLevel2=1; iLevel2<=${NLevels2}; iLevel2++ )); do
        if [ ${iLevel2} -ge ${MinLevel2Temp} ] && [ ${iLevel2} -le ${MaxLevel2Temp} ]; then
          echo "  [Clean]: --- Molecule 1, Level/Bin " ${iLevels1} " ----------------------------- "
          echo "  [Clean]: ----- Molecule 2, Level/Bin = " ${iLevels2} " --------------------- "
          echo "  [Clean]"

      
          if [ ${TranFlg} -eq 0 ]; then 
            export COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"E_"${Tran%.*}"_T_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
          else
            export COARSEAIR_BIN_OUTPUT_DIR=${COARSEAIR_OUTPUT_DIR}/"T_"${Tran%.*}"_"${Tint%.*}/"Bins_"${iLevels1}"_"${iLevels2}
          fi


          if [ -f ${COARSEAIR_BIN_OUTPUT_DIR}/NConvTraj.dat ]; then
            rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/Node*
            rm -rf ${COARSEAIR_BIN_OUTPUT_DIR}/statistics-*
          fi


        fi
        echo "  [Clean]: ---------------------------------------------------------- "
      done
    fi
    echo "  [Clean]: ------------------------------------------------------------ "
    echo " "
  done     
      
}
#================================================================================================================================#
