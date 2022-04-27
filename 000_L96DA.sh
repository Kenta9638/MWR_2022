#!/bin/sh -f

set -e
ulimit -s unlimited

#--- Setting 1 ---
F90=gfortran
#F90=ifort
MY_LAPACK1=-I/homes/metogra/kkurosaw/WORK/DA/lorenz_DA/KK_L96DA_FORTRAN/20220214/LAPACK/test_include
MY_LAPACK2=/homes/metogra/kkurosaw/WORK/DA/lorenz_DA/KK_L96DA_FORTRAN/20220214/LAPACK/test_lib
MY_LAPACK3="$MY_LAPACK2/liblapack95.a $MY_LAPACK2/liblapack.a $MY_LAPACK2/librefblas.a"

#--- Setting 2 ---
CWD=`pwd`
export OMP_NUM_THREADS=20

#--- Setting 3 ---
PARA_DAMODE=LPF
LINEAR_CASE=linear
LINEAR_CASE=nonlinear1
LINEAR_CASE=nonlinear2
#LINEAR_CASE=mix1 
#LINEAR_CASE=mix2 
#LINEAR_CASE=mix3 
PARA_PARADIR=parq1
PARA_MEM="20"
PARA_INF="`seq 0.3 0.1 0.3`"
PARA_LOC="`seq 5 1 5`"
PARA_OBS_num="`seq 20 10 20`"
PARA_Outper=0.05
PARA_dt=0.05
PARA_LPF_Neff_para=0.50
PARA_LPF_ini_minres="999"
#PARA_LPF_ini_minres="0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 999"

#EXP_NAME=20210622/$LINEAR_CASE/$PARA_LPF_adaptive_hybrid/
#EXP_NAME=test/$LINEAR_CASE/
#EXP_NAME=20220401_1st_variable/$LINEAR_CASE/
#EXP_NAME=20220401/obs_void/$LINEAR_CASE/
#EXP_NAME=20220401/model_bias/$LINEAR_CASE/
#EXP_NAME=20220210/model_bias/$LINEAR_CASE/
#EXP_NAME=20220210/obs_mix/$LINEAR_CASE/
EXP_NAME=20220210/obs_avoid/$LINEAR_CASE/

#--- Setting 4 ---
tmp1=/Outper-`printf "%3.2f" $PARA_Outper`
tmp2=/dt-`printf "%3.2f" $PARA_dt`
if [ $LINEAR_CASE == "linear" ]; then 
  PARA_linear_case=0
  PARA_R=1
elif [ $LINEAR_CASE == "nonlinear1" ]; then 
  PARA_linear_case=1
  PARA_R=1
elif [ $LINEAR_CASE == "nonlinear2" ]; then 
  PARA_linear_case=2
  PARA_R=0.1
elif [ $LINEAR_CASE == "mix1" ]; then 
  PARA_linear_case=91
  PARA_R=1
elif [ $LINEAR_CASE == "mix2" ]; then 
  PARA_linear_case=92
  PARA_R=999
elif [ $LINEAR_CASE == "mix3" ]; then 
  PARA_linear_case=93
  PARA_R=999
else
  echo ERROR!!
  exit
fi

EXP_NAME=$EXP_NAME/$PARA_DAMODE/$tmp1/
BASE_DIR=/Users/kurosawa/lorenz_DA/KK_L96DA_FORTRAN
WORK_DIR=$BASE_DIR/WORK/$EXP_NAME/$PARA_PARADIR

OUT_BASE_DIR=$BASE_DIR/OUT/$EXP_NAME
TRUTH_DIR=$BASE_DIR/TRUTH/$tmp2/
OBS_DIR=$BASE_DIR/OBS/$tmp2/
INI_DIR=$BASE_DIR/INI

if [ ! -d $WORK_DIR ]; then mkdir -p $WORK_DIR; fi
if [ ! -d $OUT_BASE_DIR  ]; then mkdir -p $OUT_BASE_DIR ; fi

#--- Compile ---
cd $WORK_DIR
# TK
cp -pr $CWD/../common/SFMT.f90         $WORK_DIR
cp -pr $CWD/../common/common.f90       $WORK_DIR
cp -pr $CWD/../common/netlib.f         $WORK_DIR
cp -pr $CWD/../common/common_mtx.f90   $WORK_DIR
# KK
cp -pr $CWD/../common/common_cdf.f90   $WORK_DIR
cp -pr $CWD/../common/prob.f90         $WORK_DIR
cp -pr $CWD/../common/SW_test.f90      $WORK_DIR
cp -pr $CWD/functions/common/common_L96_ver2.f90        $WORK_DIR
cp -pr $CWD/functions/DA/common/common_L96DA.f90        $WORK_DIR
#cp -pr $CWD/functions/DA/common/FDVAR_tools.f90         $WORK_DIR
cp -pr $CWD/functions/DA/common/Serial_EnSRF_tools.f90  $WORK_DIR
cp -pr $CWD/functions/DA/common/LPF_tools.f90           $WORK_DIR
#cp -pr $CWD/functions/DA/E4DVAR.f90                     $WORK_DIR
#cp -pr $CWD/functions/DA/FDEnVAR.f90                    $WORK_DIR
#cp -pr $CWD/functions/DA/ESMDA.f90                      $WORK_DIR
#cp -pr $CWD/functions/DA/ESMDA_fixedlag.f90             $WORK_DIR
cp -pr $CWD/functions/DA/LPF.f90                        $WORK_DIR
cp -pr $CWD/L96DA_main.f90                              $WORK_DIR
cp -pr $CWD/001_L96DA_para.sh                           $WORK_DIR

cd $WORK_DIR

echo  $WORK_DIR

MY_COMPILE1="L96DA SFMT.f90 common.f90 netlib.f common_mtx.f90 common_cdf.f90 prob.f90 SW_test.f90"
MY_COMPILE2="common_L96DA.f90 common_L96_ver2.f90"
MY_COMPILE3="Serial_EnSRF_tools.f90 LPF_tools.f90"
MY_COMPILE4="LPF.f90"
MY_COMPILE5="L96DA_main.f90"
$F90 -fopenmp -o $MY_COMPILE1 $MY_COMPILE2 $MY_COMPILE3 $MY_COMPILE4 $MY_COMPILE5 $MY_LAPACK1 $MY_LAPACK3

#--- Main Loop ---
for LOOP_MEM in $PARA_MEM
do
for LOOP_INF in $PARA_INF
do
for LOOP_LOC in $PARA_LOC
do
for LOOP_OBS in $PARA_OBS_num
do
for LOOP_MRS in $PARA_LPF_ini_minres
do
  echo $LOOP_LOC/$LOOP_INF/$LOOP_MEM
  tmp1=/MEM-`printf "%0*d"  3 $LOOP_MEM`
  tmp2=/INF-`printf "%3.2f"   $LOOP_INF`
  tmp3=/LOC-`printf "%0*d"  3 $LOOP_LOC`
#  tmp6=/OBSnum-`printf "%0*d"  3 $LOOP_OBS`
  tmp7=/MRS-`printf "%3.2f"   $LOOP_MRS`
  OUT_DIR=$OUT_BASE_DIR/$tmp1/$tmp2/$tmp3/$tmp7
  echo $OUT_DIR
  if [ ! -d $OUT_DIR  ]; then mkdir -p $OUT_DIR ; fi
  
  PARA_LPF_alpha=$LOOP_INF

  tmp_INP1="$LOOP_MEM $LOOP_INF $LOOP_LOC"                   # loop para
  tmp_INP2="$OUT_DIR $TRUTH_DIR $OBS_DIR $INI_DIR"           # directory
  tmp_INP4="$PARA_Outper $PARA_DAMODE $PARA_dt"              # L96DA para
  tmp_INP5="$LOOP_OBS $PARA_linear_case $PARA_R"             # obs
  tmp_INP6="$PARA_LPF_Neff_para $PARA_LPF_alpha"             # others1
  tmp_INP7="$LOOP_MRS"                                       # others2

  #- namelist -
  sh 001_L96DA_para.sh $tmp_INP1 $tmp_INP2 $tmp_INP3 $tmp_INP4 $tmp_INP5 $tmp_INP6 $tmp_INP7
  
  cp -p L96DA_para.cnf $OUT_DIR
 
  #- main -
  ./L96DA
  
done
done
done
done
done






exit

exit


