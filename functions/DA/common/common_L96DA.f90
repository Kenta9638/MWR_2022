MODULE common_L96DA

!======================================================================================
  USE common
  IMPLICIT NONE
  PUBLIC

!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  !--- Parameter (DIRECTORY & FILE) ---
  CHARACTER(256),SAVE :: OUT_DIR
  CHARACTER(256),SAVE :: TRUTH_DIR
  CHARACTER(256),SAVE :: OBS_DIR
  CHARACTER(256),SAVE :: INI_DIR

  !--- Parameter (L96) ---
  INTEGER,SAVE        :: PARA_L96_J
  INTEGER,SAVE        :: PARA_L96_calyear,PARA_L96_spinup  ! unit: year
  INTEGER,SAVE        :: PARA_L96_oneyear                  ! unit: day
  INTEGER,SAVE        :: PARA_L96_NWRT,PARA_L96_NTIMES     ! unit: number of timestep
  REAL(r_size),SAVE   :: PARA_L96_F,PARA_L96_P
  REAL(r_size),SAVE   :: PARA_L96_dt,PARA_L96_oneday,PARA_L96_outper ! unit: timestep
  REAL(r_size),SAVE   :: PARA_L96_bias
  CHARACTER(256),SAVE :: PARA_L96_NL_M,PARA_L96_TL_M,PARA_L96_AD_M
  
  !--- Parameter (Truth) ---
  LOGICAL,SAVE :: PARA_truth_create
  
  !--- Parameter (OBS) ---
  INTEGER,SAVE      :: PARA_obs_num
  REAL(r_size),SAVE :: PARA_obs_bias
  REAL(r_size),SAVE :: PARA_obs_avail_hr ! unit: timestep
  LOGICAL,SAVE      :: PARA_obs_stable,PARA_obs_create
  INTEGER,SAVE      :: PARA_obs_nonlinear
  
  !--- Parameter (DA) ---
  INTEGER,SAVE        :: PARA_DA_MEM
  INTEGER,SAVE        :: PARA_DA_calyear                        ! unit: year
  INTEGER,SAVE        :: PARA_DA_spinupday                      ! unit: day
  INTEGER,SAVE        :: PARA_DA_TIMES,PARA_DA_DAwindow_times   ! unit: number of timestep
  INTEGER,SAVE        :: PARA_DA_DAwindow_obstimes
  INTEGER,SAVE        :: PARA_DA_LAG,PARA_DA_Na,PARA_DA_outer_loop 
  INTEGER,SAVE        :: PARA_DA_spinupstep 
  REAL(r_size),SAVE   :: PARA_DA_DAwindow_hr               ! unit: timestep
  REAL(r_size),SAVE   :: PARA_DA_VAR_B,PARA_DA_tmp_Q,PARA_DA_R_std
  REAL(r_size),SAVE   :: PARA_DA_LPF_Neff_para,PARA_DA_LPF_Neff,PARA_DA_LPF_alpha
  REAL(r_size)        :: PARA_DA_infl,PARA_DA_local_scale,PARA_DA_hybrid_beta !<-- CAUTION (not SAVE) !!
  LOGICAL,SAVE        :: PARA_DA_FDEnVAR_HTL_FLAG
  LOGICAL,SAVE        :: PARA_DA_LPF_adaptive_stop
  REAL(r_size),SAVE   :: PARA_DA_LPF_ini_minres, PARA_DA_LPF_SWT_rate

  CHARACTER(256),SAVE :: PARA_DA_MODE
  LOGICAL,SAVE        :: PARA_DA_adap_infl,PARA_DA_localize
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:),SAVE :: PARA_DA_Q,PARA_DA_R

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_L96DA
  IMPLICIT NONE
  INTEGER :: i
  namelist / DIR_FILE_list         / OUT_DIR,TRUTH_DIR,OBS_DIR,INI_DIR
  namelist / L96DA_para_list_L96   / PARA_L96_J,PARA_L96_F,PARA_L96_P,PARA_L96_dt,     &
                                   & PARA_L96_calyear,PARA_L96_spinup,PARA_L96_oneday, &
                                   & PARA_L96_oneyear,PARA_L96_outper,PARA_L96_bias,   &
                                   & PARA_L96_NL_M,PARA_L96_TL_M
  namelist / L96DA_para_list_truth / PARA_truth_create
  namelist / L96DA_para_list_obs   / PARA_obs_num,PARA_obs_bias, &
                                   & PARA_obs_stable,PARA_obs_create,PARA_obs_nonlinear
  namelist / L96DA_para_list_DA    / PARA_DA_MODE,PARA_DA_calyear,PARA_DA_spinupday, &
                                   & PARA_DA_MEM,                               &
                                   & PARA_DA_VAR_B,PARA_DA_tmp_Q,PARA_DA_R_std, &
                                   & PARA_DA_infl,PARA_DA_adap_infl,            &
                                   & PARA_DA_localize,PARA_DA_local_scale,      &
                                   & PARA_DA_LPF_Neff_para, PARA_DA_LPF_alpha,  &
                                   & PARA_DA_LPF_ini_minres, PARA_DA_LPF_SWT_rate


  !-----------------------------------------------------------------------
  !--- READ NAMELIST ---
  OPEN(1,FILE='L96DA_para.cnf')
  READ(1,NML=DIR_FILE_list)
  READ(1,NML=L96DA_para_list_L96)
  READ(1,NML=L96DA_para_list_truth)
  READ(1,NML=L96DA_para_list_obs)
  READ(1,NML=L96DA_para_list_DA)
  CLOSE(1)
  !- L96 -
  PARA_L96_NWRT   = int(PARA_L96_outper/PARA_L96_dt)
  PARA_L96_NTIMES = int(dble(PARA_L96_oneyear*PARA_L96_calyear)*PARA_L96_oneday/PARA_L96_dt)
  PARA_L96_AD_M   = PARA_L96_TL_M
  !- OBS -
  PARA_obs_avail_hr = PARA_L96_outper
  !- DA -
  PARA_DA_TIMES = int(dble(PARA_L96_oneyear*PARA_DA_calyear)*PARA_L96_oneday/PARA_L96_outper) 
  PARA_DA_DAwindow_hr = dble(PARA_DA_LAG)*PARA_L96_outper
  PARA_DA_DAwindow_obstimes = int(PARA_DA_DAwindow_hr/PARA_OBS_avail_hr)    ! including t=0
  PARA_DA_DAwindow_times    = PARA_DA_DAwindow_obstimes*PARA_L96_NWRT       ! including t=0
  ALLOCATE(PARA_DA_Q(PARA_L96_J,PARA_L96_J)); PARA_DA_Q = PARA_DA_tmp_Q
  ALLOCATE(PARA_DA_R(PARA_obs_num,PARA_obs_num)); PARA_DA_R = 0.0d0

  IF (PARA_obs_nonlinear < 90) THEN
    DO i = 1,PARA_obs_num
      PARA_DA_R(i,i) = 1.0d0*PARA_DA_R_std**2
    ENDDO
  ELSEIF (PARA_obs_nonlinear==91) THEN
    DO i = 1,20
      PARA_DA_R(i,i) = 1.0d0
    ENDDO
  ELSEIF (PARA_obs_nonlinear==92) THEN
    DO i = 1,10
      PARA_DA_R(i,i) = 1.0d0
    ENDDO
    DO i = 11,20
      PARA_DA_R(i,i) = 1.0d0*0.10d0**2
    ENDDO
  ELSEIF (PARA_obs_nonlinear==93) THEN
    DO i = 1,10
      PARA_DA_R(i,i) = 1.0d0
    ENDDO
    DO i = 11,20
      PARA_DA_R(i,i) = 1.0d0*0.10d0**2
    ENDDO
  ENDIF

  PARA_DA_spinupstep = int(PARA_L96_oneday/PARA_L96_outper)*dble(PARA_DA_spinupday) ! spinup tekitou
  PARA_DA_LPF_Neff = PARA_DA_LPF_Neff_para*PARA_DA_MEM

  IF (PARA_DA_LPF_ini_minres > 1.0d0) THEN
    PARA_DA_LPF_adaptive_stop = .TRUE.
  ELSE
    PARA_DA_LPF_adaptive_stop = .FALSE.
  ENDIF

  RETURN
END SUBROUTINE set_common_L96DA

!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_para
  !- monitor -
  PRINT '(A)',' '
  PRINT '(A)','--------------------------------'
  PRINT '(A)','L96'
  PRINT '(A,I8)'  ,' PARA_L96_J                = ',PARA_L96_J
  PRINT '(A,F8.3)',' PARA_L96_F                = ',PARA_L96_F
  PRINT '(A,F8.3)',' PARA_L96_P                = ',PARA_L96_P
  PRINT '(A,F8.3)',' PARA_L96_dt               = ',PARA_L96_dt
  PRINT '(A,I8)'  ,' PARA_L96_calyear          = ',PARA_L96_calyear
  PRINT '(A,I8)'  ,' PARA_L96_spinup           = ',PARA_L96_spinup
  PRINT '(A,F8.3)',' PARA_L96_oneday           = ',PARA_L96_oneday
  PRINT '(A,I8)'  ,' PARA_L96_oneyear          = ',PARA_L96_oneyear
  PRINT '(A,F8.3)',' PARA_L96_outper           = ',PARA_L96_outper
  PRINT '(A,I8)'  ,' PARA_L96_NWRT             = ',PARA_L96_NWRT
  PRINT '(A,I8)'  ,' PARA_L96_NTIMES           = ',PARA_L96_NTIMES
  PRINT '(A,F8.3)',' PARA_L96_bias             = ',PARA_L96_bias
  PRINT '(A,A)'   ,' PARA_L96_NL_M             = ',trim(PARA_L96_NL_M)
  PRINT '(A,A)'   ,' PARA_L96_TL_M             = ',trim(PARA_L96_TL_M)
  PRINT '(A,A)'   ,' PARA_L96_AD_M             = ',trim(PARA_L96_AD_M)
  PRINT '(A)','--------------------------------'
  PRINT '(A)','TRUTH'
  PRINT *,'PARA_truth_create         = ',PARA_truth_create
  PRINT '(A)','--------------------------------'
  PRINT '(A)','OBS'
  PRINT '(A,I8)'  ,' PARA_obs_num              = ',PARA_obs_num
  PRINT '(A,F8.3)',' PARA_obs_avail_hr         = ',PARA_obs_avail_hr
  PRINT '(A,F8.3)',' PARA_obs_bias             = ',PARA_obs_bias
  PRINT *,'PARA_obs_stable           = ',PARA_obs_stable
  PRINT *,'PARA_obs_create           = ',PARA_obs_create
  PRINT '(A,I8)'  ,' PARA_obs_nonlinear        = ',PARA_obs_nonlinear
  PRINT '(A)','--------------------------------'
  PRINT '(A)','DA'
  PRINT '(A,I8)'  ,' PARA_DA_calyear           = ',PARA_DA_calyear
  PRINT '(A,I8)'  ,' PARA_DA_spinupday         = ',PARA_DA_spinupday
  PRINT '(A,I8)'  ,' PARA_DA_spinupstep        = ',PARA_DA_spinupstep
  PRINT '(A,I8)'  ,' PARA_DA_TIMES             = ',PARA_DA_TIMES
  PRINT '(A,I8)'  ,' PARA_DA_MEM               = ',PARA_DA_MEM
  PRINT '(A,F8.3)',' PARA_DA_VAR_B             = ',PARA_DA_VAR_B
  PRINT '(A,F8.3)',' PARA_DA_Q                 = ',PARA_DA_tmp_Q ! <-- CAUTION !!
  PRINT '(A,F8.3)',' PARA_DA_R_std             = ',PARA_DA_R_std
  PRINT '(A,F8.3)',' PARA_DA_R                 = ',PARA_DA_R_std ! <-- CAUTION !!
  PRINT '(A,F8.3)',' PARA_DA_infl              = ',PARA_DA_infl
  PRINT '(A,F8.3)',' PARA_DA_local_scale       = ',PARA_DA_local_scale
  PRINT '(A,F8.3)',' PARA_DA_LPF_Neff          = ',PARA_DA_LPF_Neff
  PRINT '(A,F8.3)',' PARA_DA_LPF_alpha         = ',PARA_DA_LPF_alpha
  PRINT *,'PARA_DA_adap_infl         = ',PARA_DA_adap_infl
  PRINT *,'PARA_DA_localize          = ',PARA_DA_localize
  
  RETURN
END SUBROUTINE monit_para

!=======================================================================
! Localization
!=======================================================================
!-----------------------------------------------------------------------
! Schur Product ! from letkf-master common_enkf.f90
!-----------------------------------------------------------------------
SUBROUTINE enkf_schur(scale,dist,factor)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: scale
  REAL(r_size),INTENT(IN) :: dist
  REAL(r_size),INTENT(OUT) :: factor
  REAL(r_size) :: a,b

  a = scale * SQRT(10.0d0/3.0d0)
  b = dist / a

  IF( dist <= a ) THEN
    factor = 1.0d0 -0.25d0*b**5 + 0.5d0*b**4 + 5.0d0/8.0d0*b**3 &
      & - 5.0d0/3.0d0*b**2
  ELSE IF( dist <= 2*a ) THEN
    factor = 1.0d0/12.0d0*b**5 - 0.5d0*b**4 + 5.0d0/8.0d0*b**3 &
      & + 5.0d0/3.0d0*b**2 - 5.0d0*b + 4.0d0 - 2.0d0/3.0d0/b
  ELSE
    factor = 0.0d0
  END IF

  RETURN
END SUBROUTINE enkf_schur

!=======================================================================
! Obsope
!=======================================================================
!-----------------------------------------------------------------------
! KK_FUNC_SET_OBS form MATLAB
!-----------------------------------------------------------------------
SUBROUTINE SET_OBS(obs_var,obs_loc,Yo,H,obsnum)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN)  :: obs_var(PARA_obs_num)
  INTEGER     ,INTENT(IN)  :: obs_loc(PARA_obs_num)
  REAL(r_size),INTENT(OUT) :: Yo(PARA_obs_num)
  INTEGER     ,INTENT(OUT) :: H(PARA_obs_num,PARA_L96_J)
  INTEGER     ,INTENT(OUT) :: obsnum
  INTEGER :: i
  INTEGER :: tmp_num
 
  !--- set Yo & H ---
  Yo     = obs_var
  H(:,:) = 0.0d0
  obsnum = 0
  DO i = 1,PARA_OBS_num
    tmp_num = obs_loc(i)
    IF (tmp_num>0) THEN
      H(i,tmp_num) = 1
      obsnum = obsnum + 1
    ENDIF
  ENDDO

  RETURN
END SUBROUTINE SET_OBS

SUBROUTINE infl_RTPP(xf,infl,xa)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN)    :: xf  (PARA_L96_J,PARA_DA_MEM)
  REAL(r_size),INTENT(IN)    :: infl
  REAL(r_size),INTENT(INOUT) :: xa  (PARA_L96_J,PARA_DA_MEM)
!  REAL(r_size),INTENT(OUT)   :: Pa  (PARA_L96_J,PARA_L96_J)
  REAL(r_size) :: xf_mean(PARA_L96_J)
  REAL(r_size) :: xa_mean(PARA_L96_J)
  REAL(r_size) :: Ef     (PARA_L96_J,PARA_DA_MEM)
  REAL(r_size) :: Ea     (PARA_L96_J,PARA_DA_MEM)
  REAL(r_size) :: Einfl  (PARA_L96_J,PARA_DA_MEM)
  REAL(r_size) :: dxa    (PARA_L96_J,PARA_DA_MEM)
  INTEGER :: i,j,k,l,m

  DO j = 1,PARA_L96_J
    CALL com_mean(PARA_DA_MEM,xf(j,:),xf_mean(j))
    CALL com_mean(PARA_DA_MEM,xa(j,:),xa_mean(j))
  ENDDO
  DO j = 1,PARA_DA_MEM
    Ef(:,j) = (xf(:,j)-xf_mean)/sqrt(dble(PARA_DA_MEM-1))
    Ea(:,j) = (xa(:,j)-xa_mean)/sqrt(dble(PARA_DA_MEM-1))
  ENDDO
  DO j = 1,PARA_L96_J
    DO m = 1,PARA_DA_MEM
      Einfl(j,m) = Ea(j,m)*(1.0d0-infl)+Ef(j,m)*infl
      xa   (j,m) = xa_mean(j) + Einfl(j,m)*sqrt(dble(PARA_DA_MEM-1))
    ENDDO
  ENDDO
  DO j = 1,PARA_DA_MEM
    dxa(:,j) = xa(:,j) - xa_mean
  ENDDO
!  Pa(:,:) = matmul(dxa,transpose(dxa))/dble(PARA_DA_MEM-1)
return
END SUBROUTINE infl_RTPP

END MODULE common_L96DA
