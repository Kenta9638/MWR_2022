#!/bin/sh -f

set -e

INP01=$1   # OUT Directory
INP02=$2   # OBS Directory
INP03=$3   # INITIAL Directory

INP04=$4       # Member
INP05=$5       # Inflation
INP06=$6       # Localization Scale
INP07=$7       # OBS num

INP08=$8       # Neff
INP09=$9       # RMm
INP10=${10}    # RMs 
INP11=${11}    # VMm
INP12=${12}    # VMs
INP13=${13}    # POb
INP14=${14}    # POs
INP15=${15}    # Mra
INP16=${16}    # OBs
INP17=${17}    # RAf
INP18=${18}    # IOt
INP19=${19}    # JOt

INP20=${20}    # PARA_LOOP_NUM
INP21=${21}    # FLAG_READ_MINRES
INP22=${22}    # minres_file


if [ -e rankine_para.cnf ]; then rm rankine_para.cnf; fi
cat >> rankine_para.cnf << EOF
!
!===== DIRECTORY & FILE =====-
&DIR_FILE_list
OUT_DIR   = '$INP01'
OBS_DIR   = '$INP02'
INI_DIR   = '$INP03'
/
!===== rankine vortex =====
&rankine_para_list
PARA_MEM  = $INP04,
PARA_INF  = $INP05,
PARA_LOC  = $INP06,
PARA_OBn  = $INP07,
!
PARA_Nef  = $INP08,
PARA_RMm  = $INP09,
PARA_RMs  = $INP10,
PARA_VMm  = $INP11,
PARA_VMs  = $INP12,
PARA_POb  = $INP13,
PARA_POs  = $INP14,
PARA_Mra  = $INP15,
PARA_OBs  = $INP16,
PARA_RAf  = $INP17,
PARA_IOt  = $INP18,
PARA_JOt  = $INP19,
!
PARA_LOOP_NUM    = $INP20,
FLAG_READ_MINRES = $INP21,
MINRES_NCDIR     = '$INP22',
/
EOF


exit

