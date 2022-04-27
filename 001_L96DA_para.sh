#!/bin/sh -f

set -e

INP1=$1       # Member
INP2=$2       # Inflation
INP3=$3       # Localization Scale

INP11=$4   # OUT Directory
INP12=$5   # TRUTH Directory
INP13=$6   # OBS Directory
INP14=$7   # INITIAL Directory

INP31=${8} # Outper
INP32=${9} # DA mode
INP33=${10} # dt

INP41=${11} # OBS_num
INP42=${12} # Linear_CASE
INP43=${13} # R_std

INP51=${14} # LPF Neff para
INP52=${15} # LPF alpha

INP61=${16} # LPF minres


if [ -e L96DA_para.cnf ]; then rm L96DA_para.cnf; fi
cat >> L96DA_para.cnf << EOF
!
!===== DIRECTORY & FILE =====-
&DIR_FILE_list
OUT_DIR   = '$INP11'
TRUTH_DIR = '$INP12'
OBS_DIR   = '$INP13'
INI_DIR   = '$INP14'
/
!
!===== Lorenz96 =====-
&L96DA_para_list_L96
PARA_L96_J       = 40,
PARA_L96_F       = 9.00,
PARA_L96_P       = 0.008,
PARA_L96_dt      = $INP33,     ! unit:time_step, (dt in Model cal.)
PARA_L96_calyear = 101,        ! unit:year,      (include spin-up)
PARA_L96_spinup  = 1,          ! unit:year
PARA_L96_oneday  = 0.2,        ! unit:time_step, (time_step/1day in Lorenz96)
PARA_L96_oneyear = 365,        ! unit:day
PARA_L96_outper  = $INP31,     ! unit:time_step, (default:out every 6hr (0.05 timestep, oneday/4))
PARA_L96_bias    = 0.0,        ! (default:0.0)
PARA_L96_NL_M    = 'JP',       ! Kenta Kurosawa or Jonathan Poterjoy
PARA_L96_TL_M    = 'JP',       ! Kenta Kurosawa or Jonathan Poterjoy
/
!
!===== Truth =====
&L96DA_para_list_truth
PARA_truth_create = .false.
/
!
!===== Observation =====
&L96DA_para_list_obs
PARA_OBS_num     = $INP41,
PARA_OBS_bias    = 0.0,
PARA_obs_stable  = .TRUE.
PARA_obs_create  = .false.
PARA_obs_nonlinear = $INP42,   ! 1: Y=H*(x^2), 2: Y=log(abs(x))
!!! PARA_OBS_avail_hr
/
!
!===== Data Assimilation =====
&L96DA_para_list_DA
PARA_DA_MODE       = '$INP32',
PARA_DA_calyear    = 10,     ! unit:year
PARA_DA_spinupday  = 300,   ! unit:day
PARA_DA_MEM        = $INP1,
PARA_DA_VAR_B      = 999,   ! Only for 3DVAR, 4DVAR, or En4DVAR
PARA_DA_tmp_Q      = 0,      ! (default:zero)
PARA_DA_R_std      = $INP43,      ! (default:one)
!
!--- EnKF ---
PARA_DA_infl         = $INP2,
PARA_DA_adap_infl    = .false.,
PARA_DA_localize     = .true.,
PARA_DA_local_scale  = $INP3,
!
!--- LPF ---
PARA_DA_LPF_Neff_para       = $INP51,
PARA_DA_LPF_alpha           = $INP52,
PARA_DA_LPF_ini_minres      = $INP61,
PARA_DA_LPF_SWT_rate        = 1.0d0,
/
EOF


exit


