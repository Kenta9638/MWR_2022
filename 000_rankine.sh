#!/bin/sh -f

set -e
ulimit -s unlimited

#--- Setting 1 ---
CWD=`pwd`
F90=gfortran
F90=mpif90
MY_LAPACK1=-I/$CWD/../LAPACK/test_include
MY_LAPACK2=/$CWD/../LAPACK/test_lib
MY_LAPACK3="$MY_LAPACK2/liblapack95.a $MY_LAPACK2/liblapack.a $MY_LAPACK2/librefblas.a"
MY_LAPACK4="-L/usr/lib64/openmpi/lib -lnetcdff -lnetcdf"
MY_PASS1=/homes/metogra/kkurosaw/WORK/TOOLS/kk_tools/kk_tools/netcdf_common/

#--- Setting 2 ---
export OMP_NUM_THREADS=1

#--- Setting 3 ---
PARA_PARADIR=para1

PARA_MEM="40"
PARA_LOC="`seq 2000 600 2000`"
PARA_OBn="`seq 100 20 100`"
PARA_RMm=12     # Rmax mean
PARA_RMs="3" # "`seq 0 5 10`"      # Rmax sig
PARA_VMm=30     # Vmax mean
PARA_VMs="3"      # Vmax sig
PARA_POb=6      # Position bias
PARA_POs="12" # 12 18 24 36"      # Position sig
PARA_Mra=45     # Maximum radius
PARA_OBs="3"      # Observation sig

PARA_RAf=1      # Set flag to 1 to assimilate obs located in the vicinity of the radar location.
PARA_IOt=25     # Radar location
PARA_JOt=25     # Radar location

PARA_INF="`seq 0.3 0.1 0.3`"
PARA_Nef=0.50

PARA_LOOP_NUM=1
FLAG_READ_MINRES=.false.

#--- Setting 4 ---
#EXP_NAME=rankine_20211121/POb6_convert/
EXP_NAME=test/
EXP_NAME=$EXP_NAME/$PARA_DAMODE/
BASE_DIR=/Users/kurosawa/rankine/fortran/
WORK_DIR=$BASE_DIR/WORK/$EXP_NAME/$PARA_PARADIR
OUT_BASE_DIR=$BASE_DIR/OUT/$EXP_NAME
OBS_DIR=$OUT_BASE_DIR/
INI_DIR=$OUT_BASE_DIR/

EXP_NAME_minres=rankine_20211121/POb6/
EXP_NAME_minres=$EXP_NAME_minres/$PARA_DAMODE/
OUT_BASE_DIR_minres=$BASE_DIR/OUT/$EXP_NAME_minres

if [ ! -d $WORK_DIR ]; then mkdir -p $WORK_DIR; fi
if [ ! -d $OUT_BASE_DIR  ]; then mkdir -p $OUT_BASE_DIR ; fi

#--- Compile ---
cd $WORK_DIR
# TK
cp -pr $CWD/../DA_common/SFMT.f90         $WORK_DIR
cp -pr $CWD/../DA_common/common.f90       $WORK_DIR
cp -pr $CWD/../DA_common/common_mpi.f90   $WORK_DIR
cp -pr $CWD/../DA_common/netlib.f         $WORK_DIR
cp -pr $CWD/../DA_common/common_mtx.f90   $WORK_DIR
# KK
cp -pr $CWD/../SW_common/common_cdf.f90   $WORK_DIR
cp -pr $CWD/../SW_common/prob.f90         $WORK_DIR
cp -pr $CWD/../SW_common/SW_test.f90      $WORK_DIR
cp -pr $CWD/functions/common_rankine.f90        $WORK_DIR
cp -pr $CWD/functions/Serial_EnSRF_tools.f90  $WORK_DIR
cp -pr $CWD/functions/LPF_tools.f90           $WORK_DIR
cp -pr $CWD/rankine_main.f90                              $WORK_DIR
cp -pr $CWD/001_rankine_para.sh                           $WORK_DIR
# netcdf
cp -pr $MY_PASS1/kk_netcdf_common.f90        $WORK_DIR
cp -pr $MY_PASS1/kk_netcdf_tools.f90         $WORK_DIR


cd $WORK_DIR

echo  $WORK_DIR

MY_COMPILE1="rankine SFMT.f90 common.f90 common_mpi.f90 netlib.f common_mtx.f90 common_cdf.f90 prob.f90 SW_test.f90"
MY_COMPILE2="common_rankine.f90"
MY_COMPILE3="Serial_EnSRF_tools.f90 LPF_tools.f90"
MY_COMPILE4="kk_netcdf_common.f90 kk_netcdf_tools.f90"
MY_COMPILE5="rankine_main.f90"
$F90 -fopenmp -o $MY_COMPILE1 $MY_COMPILE2 $MY_COMPILE3 $MY_COMPILE4 $MY_COMPILE5 $MY_LAPACK1 $MY_LAPACK3 $MY_LAPACK4

#--- Main Loop ---
for LOOP_MEM in $PARA_MEM
do
for LOOP_INF in $PARA_INF
do
for LOOP_LOC in $PARA_LOC
do
for LOOP_OBS in $PARA_OBn
do
for LOOP_RMs in $PARA_RMs
do
for LOOP_VMs in $PARA_VMs
do
for LOOP_POs in $PARA_POs
do
for LOOP_OBs in $PARA_OBs
do
  echo $LOOP_LOC/$LOOP_INF/$LOOP_MEM
  echo `date`
  tmp1=/MEM-`printf "%0*d"  3 $LOOP_MEM`
  tmp2=/INF-`printf "%3.2f"   $LOOP_INF`
  tmp3=/LOC-`printf "%0*d"  4 $LOOP_LOC`
  tmp4=/OBn-`printf "%0*d"  3 $LOOP_OBS`
  tmp5=/RMs-`printf "%0*d"  3 $LOOP_RMs`
  tmp6=/VMs-`printf "%0*d"  3 $LOOP_VMs`
  tmp7=/POs-`printf "%0*d"  3 $LOOP_POs`
  tmp8=/OBs-`printf "%0*d"  3 $LOOP_OBs`
  OUT_DIR=$OUT_BASE_DIR/$tmp1/$tmp2/$tmp3/$tmp4/$tmp5/$tmp6/$tmp7/$tmp8/
  echo $OUT_DIR
  if [ ! -d $OUT_DIR  ]; then mkdir -p $OUT_DIR ; fi
  OUT_DIR_minres=$OUT_BASE_DIR_minres/$tmp1/$tmp2/$tmp3/$tmp4/$tmp5/$tmp6/$tmp7/$tmp8/
  
  tmp_INP1="$OUT_DIR $OBS_DIR $INI_DIR"                                     # directory
  tmp_INP2="$LOOP_MEM $LOOP_INF $LOOP_LOC $LOOP_OBS"                        # loop para
  tmp_INP3="$PARA_Nef $PARA_RMm $LOOP_RMs $PARA_VMm $LOOP_VMs $PARA_POb $LOOP_POs "    # 1
  tmp_INP4="$PARA_Mra $LOOP_OBs $PARA_RAf $PARA_IOt $PARA_JOt"                         # 2
  tmp_INP5="$PARA_LOOP_NUM $FLAG_READ_MINRES $OUT_DIR_minres"    # 3

  #- namelist -
  sh 001_rankine_para.sh $tmp_INP1 $tmp_INP2 $tmp_INP3 $tmp_INP4 $tmp_INP5 
  cp -p rankine_para.cnf $OUT_DIR
 
  # main -
#  ./rankine
  /usr/lib64/openmpi/bin/mpirun -np 25 ./rankine >& $OUT_DIR/output.log
#  exit

done
done
done
done
done
done
done
done






exit

exit


