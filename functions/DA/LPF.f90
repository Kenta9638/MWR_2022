MODULE LPF
!$use omp_lib

  USE common
  USE common_mtx
  USE common_L96DA
  USE common_L96
 ! use f95_lapack
  USE LPF_tools

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: LPF_cycle

CONTAINS
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE LPF_cycle(DATA_truth,DATA_obs_var,DATA_obs_loc,DATA_initial, & ! IN
                     FD_FLAG,                                           & ! OUT
                     OUT_RMSE_AN_F,OUT_SPRD_AN_F,                       & ! OUT
                     OUT_min_res_mdl,OUT_min_res_obs,OUT_rank_hstg,     & ! OUT
                     OUT_min_res_mdl_1st,OUT_min_res_obs_1st,           & ! OUT
                     OUT_spread_mdl_1st,OUT_spread_obs_1st)               ! OUT
  IMPLICIT NONE
  REAL(r_size),INTENT(IN)  :: DATA_truth  (PARA_L96_J,PARA_DA_TIMES) 
  REAL(r_size),INTENT(IN)  :: DATA_obs_var(PARA_obs_num,PARA_DA_TIMES)
  INTEGER     ,INTENT(IN)  :: DATA_obs_loc(PARA_obs_num,PARA_DA_TIMES)
  REAL(r_size),INTENT(IN)  :: DATA_initial(PARA_L96_J,PARA_DA_MEM)
  LOGICAL     ,INTENT(OUT) :: FD_FLAG ! Filter Divergence FLAG
  REAL(r_size) :: OUT_RMSE_FG_F(PARA_DA_TIMES) ! FG, Filter
  REAL(r_size) :: OUT_SPRD_FG_F(PARA_DA_TIMES) ! FG, Filter
  REAL(r_size),INTENT(OUT),OPTIONAL :: OUT_RMSE_AN_F(PARA_DA_TIMES) ! AN, Filter
  REAL(r_size),INTENT(OUT),OPTIONAL :: OUT_SPRD_AN_F(PARA_DA_TIMES) ! AN, Filter
  REAL(r_size) :: OUT_RMSE_AN_S(PARA_DA_TIMES) ! AN, Smoother
  REAL(r_size) :: OUT_SPRD_AN_S(PARA_DA_TIMES) ! AN, Smoother
!  REAL(r_size) :: OUT_MEAN_AN_S(PARA_L96_J,PARA_DA_TIMES)

  !--- LPF ---
  REAL(r_size) :: tmp_DA_Yo     (PARA_obs_num)
  REAL(r_size) :: tmp_dxf1      (PARA_L96_J,PARA_DA_MEM)
  REAL(r_size) :: tmp_dxf2      (PARA_L96_J,PARA_DA_MEM)
  REAL(r_size) :: tmp_dxa       (PARA_L96_J,PARA_DA_MEM)
  INTEGER      :: tmp_DA_H      (PARA_obs_num,PARA_L96_J)
  REAL(r_size) :: tmp_DA_hx     (PARA_obs_num,PARA_DA_MEM)
  REAL(r_size) :: tmp_DA_HC     (PARA_obs_num,PARA_L96_J)
  REAL(r_size) :: tmp_DA_HCH    (PARA_obs_num,PARA_obs_num)
  INTEGER      :: tmp_DA_obsnum 
  INTEGER      :: tmp_DA_TIMES
  INTEGER      :: tmp_tmp_DA_counter
  INTEGER      :: tmp_tmp_DA_H(PARA_L96_J)
  REAL(r_size) :: tmp_tmp_DA_Yo,tmp_tmp_DA_R
  REAL(r_size) :: OUT_xa(PARA_L96_J,PARA_DA_MEM)
  INTEGER      :: LAG,Na

  REAL(r_size),DIMENSION(PARA_L96_J,PARA_DA_MEM,PARA_DA_LAG+1) :: x_pri,x_pos  
  REAL(r_size),DIMENSION(PARA_L96_J,PARA_DA_MEM) :: tmp_DA_xf,tmp_DA_xa
  REAL(r_size),DIMENSION(PARA_L96_J,PARA_DA_MEM) :: tmp_KS_xf,tmp_DA_xf_smoother
  REAL(r_size),DIMENSION(PARA_L96_J,PARA_L96_J)  :: tmp_KS_Pf,tmp_DA_Pf

  !- adaptive hybrid
  REAL(r_size) :: min_res_mdl(PARA_L96_J), min_res_obs(PARA_obs_num)
  REAL(r_size),INTENT(OUT),OPTIONAL :: OUT_min_res_mdl(PARA_DA_TIMES) 
  REAL(r_size),INTENT(OUT),OPTIONAL :: OUT_min_res_obs(PARA_DA_TIMES) 
  REAL(r_size),INTENT(OUT),OPTIONAL :: OUT_min_res_mdl_1st(PARA_DA_TIMES) 
  REAL(r_size),INTENT(OUT),OPTIONAL :: OUT_min_res_obs_1st(PARA_DA_TIMES) 
  REAL(r_size),INTENT(OUT),OPTIONAL :: OUT_spread_mdl_1st(PARA_DA_TIMES) 
  REAL(r_size),INTENT(OUT),OPTIONAL :: OUT_spread_obs_1st(PARA_DA_TIMES) 
  !- rank histogram
  REAL(r_size) :: work1(PARA_DA_MEM), work2
  INTEGER :: work3
  INTEGER,INTENT(OUT),OPTIONAL :: OUT_rank_hstg(PARA_DA_MEM+1) 

!!  !--- LPS ---
!!  LOGICAL      :: SAVE_smoother_obs_ind  (PARA_DA_LAG+1,PARA_obs_num)
!!  INTEGER      :: SAVE_smoother_niter    (PARA_DA_LAG+1)
!!  REAL(r_size) :: SAVE_smoother_omega_mdl(PARA_DA_LAG+1,PARA_L96_J  ,PARA_DA_MEM,99,PARA_obs_num)
!!  REAL(r_size) :: SAVE_smoother_omega_obs(PARA_DA_LAG+1,PARA_obs_num,PARA_DA_MEM,99,PARA_obs_num) 
!!  INTEGER      :: SAVE_smoother_LPF_ind  (PARA_DA_LAG+1,PARA_DA_MEM,99,PARA_obs_num)
!!  REAL(r_size) :: SAVE_smoother_hx       (PARA_DA_LAG+1,PARA_obs_num,PARA_DA_MEM)
!!  REAL(r_size) :: SAVE_smoother_xf       (PARA_DA_LAG+1,PARA_L96_J  ,PARA_DA_MEM)
!!  REAL(r_size) :: SAVE_smoother_xa       (PARA_DA_LAG+1,PARA_L96_J  ,PARA_DA_MEM)
!!  INTEGER      :: SAVE_smoother_H        (PARA_DA_LAG+1,PARA_obs_num,PARA_L96_J)
!!  REAL(r_size) :: SAVE_smoother_Yo       (PARA_DA_LAG+1,PARA_obs_num)
!!  INTEGER :: latest_smoother_num, DAwin_counter 
!!  LOGICAL :: LPS_flag

  REAL(r_size) :: tmp_DA_fact         (PARA_L96_J,PARA_L96_J)
  INTEGER :: dist

  REAL(r_size),ALLOCATABLE :: tmp_matrix1(:,:)
  REAL(r_size),ALLOCATABLE :: tmp_matrix2(:)
  REAL(r_size),ALLOCATABLE :: tmp_matrix3(:)
  REAL(r_size),ALLOCATABLE :: tmp_matrix4(:,:)
  REAL(r_size),ALLOCATABLE :: tmp_matrix5(:,:)
  REAL(r_size) :: tmp_mean(PARA_L96_J)
  REAL(r_size) :: tmp_var1,tmp_var2,tmp_var3
  REAL(r_size) :: tmp_xf (PARA_L96_J,PARA_DA_MEM)
  REAL(r_size) :: save_xf(PARA_L96_J,PARA_DA_MEM)

  INTEGER :: i,j,k,l,m,n
  INTEGER :: tmp_num1,tmp_num2,tmp_num3
  REAL(r_size) :: RMSE,threshold

  CHARACTER(256) :: OUT_FILE1,OUT_FILE2,OUT_FILE3,OUT_FILE4
  CHARACTER(256) :: OUT_FILE5,OUT_FILE6,OUT_FILE7,OUT_FILE8
  CHARACTER(256) :: OUT_FILE9,OUT_FILE10,OUT_FILE11,OUT_FILE12
  CHARACTER(256) :: OUT_FILE13,OUT_FILE14,OUT_FILE15,OUT_FILE16
  LOGICAL :: EXIT_FLAG

  WRITE(*,*) ' '
  WRITE(*,*) '----------------------------------------'
  WRITE(*,*) '   Hello from LPF !!'
  WRITE(*,*) '----------------------------------------'
  
  !=============================================
  ! 1. SETTING
  !=============================================
  !--- 1.1. Setting ---
  FD_FLAG      = .FALSE.
  EXIT_FLAG    = .FALSE.
  tmp_DA_TIMES = PARA_L96_NWRT ! <-- CAUTION !!
  tmp_DA_xf    = DATA_initial
  OUT_RMSE_AN_F = UNDEF
  OUT_SPRD_AN_F = UNDEF
  OUT_RMSE_AN_S = UNDEF
  OUT_SPRD_AN_S = UNDEF
!!  latest_smoother_num = 1
  threshold = 6.0d0
  OUT_rank_hstg = 0

!!  !--- 1.2. Smoother ---
!!  LAG          = PARA_DA_LAG
!!  latest_smoother_num = 1

  !--- 1.4. Localization Factor ---
  DO i = 1,PARA_L96_J
    DO j = 1,PARA_L96_J
      dist = min(abs(i-j),PARA_L96_J-abs(i-j))
      CALL enkf_schur(PARA_DA_local_scale,dble(dist),tmp_DA_fact(i,j))
    ENDDO
  ENDDO

!kld = UNDEF
!DO i = 1,PARA_L96_J
!  CALL com_kld(PARA_DA_MEM,tmp_DA_xf(i,:),kld(i))
!ENDDO
!write(*,*) kld
!print *
!read *
!DO j = 1,10
!  tmp_DA_xf(i,j) = dble(j)+dble(i)*0.1
!ENDDO
!ENDDO

  !=============================================
  ! 2. MAIN LOOP
  !=============================================
  DO i = 1,PARA_DA_TIMES

    DO j = 1,PARA_L96_J
    DO k = 1,PARA_DA_MEM
     IF ((isnan(tmp_DA_xf(j,k))) .OR. tmp_DA_xf(j,k)>30.0d0) THEN
       EXIT_FLAG = .TRUE.
       FD_FLAG   = .TRUE.
       exit
     ENDIF
    ENDDO
    ENDDO
       

    IF (EXIT_FLAG) EXIT  ! Filter Divergence 

    !--- 2.3. Obs Data ---
    CALL SET_OBS(DATA_obs_var(:,i),DATA_obs_loc(:,i), & ! IN
               & tmp_DA_Yo                          , & ! OUT
               & tmp_DA_H                           , & ! OUT
               & tmp_DA_obsnum)                         ! OUT

!    tmp_DA_H = 0
!    tmp_DA_H(1, 6) = 1
!    tmp_DA_H(2, 9) = 1
!    tmp_DA_H(3,10) = 1
!    tmp_DA_H(4,19) = 1
!    tmp_DA_H(5,21) = 1
    
    !-- hx
    IF     (PARA_obs_nonlinear==1) THEN
      tmp_DA_hx = matmul(tmp_DA_H,tmp_DA_xf**2)
    ELSEIF (PARA_obs_nonlinear==2) THEN
      tmp_DA_hx = matmul(tmp_DA_H,log(abs(tmp_DA_xf)))
    ELSEIF (PARA_obs_nonlinear==0) THEN
      tmp_DA_hx = matmul(tmp_DA_H,tmp_DA_xf)
    ELSEIF (PARA_obs_nonlinear==91) THEN
      tmp_DA_hx( 1:10,:) = matmul(tmp_DA_H( 1:10,:),tmp_DA_xf)
      tmp_DA_hx(11:20,:) = matmul(tmp_DA_H(11:20,:),tmp_DA_xf**2)
    ELSEIF (PARA_obs_nonlinear==92) THEN
      tmp_DA_hx( 1:10,:) = matmul(tmp_DA_H( 1:10,:),tmp_DA_xf)
      tmp_DA_hx(11:20,:) = matmul(tmp_DA_H(11:20,:),log(abs(tmp_DA_xf)))
    ELSEIF (PARA_obs_nonlinear==93) THEN
      tmp_DA_hx( 1:10,:) = matmul(tmp_DA_H( 1:10,:),tmp_DA_xf**2)
      tmp_DA_hx(11:20,:) = matmul(tmp_DA_H(11:20,:),log(abs(tmp_DA_xf)))
    ENDIF
    !-- HC & HCH
    tmp_DA_HC  = matmul(dble(tmp_DA_H),tmp_DA_fact)
    tmp_DA_HCH = matmul(tmp_DA_HC,transpose(dble(tmp_DA_H)))

    !-- CORE 
    CALL LPF_CORE(tmp_DA_xf,tmp_DA_fact,                          & ! IN    (standard LPF) 
                  tmp_DA_Yo,PARA_DA_R,                            & ! IN    (standard LPF)
                  tmp_DA_H,tmp_DA_hx,tmp_DA_HC,tmp_DA_HCH,        & ! IN    (standard LPF)
                  OUT_xa,                                         & ! OUT   (standard LPF)          
                  min_res_mdl,min_res_obs)                          ! OUT (for adaptive hybrid)

    !--- score (filter) ---
    DO j = 1,PARA_L96_J
      CALL com_mean(PARA_DA_MEM,OUT_xa(j,:),tmp_mean(j))
    ENDDO
    CALL com_rms(PARA_L96_J,tmp_mean-DATA_truth(:,i),OUT_RMSE_AN_F(i))
    CALL com_sprd(PARA_L96_J,PARA_DA_MEM,OUT_xa,OUT_SPRD_AN_F(i))
    !--- min_res ---
    CALL com_mean(PARA_L96_J,min_res_mdl,OUT_min_res_mdl(i))
    CALL com_mean(PARA_obs_num-5,min_res_obs,OUT_min_res_obs(i))  ! <-- !!
    !--- for KP2022 ---
    OUT_min_res_mdl_1st(i) = min_res_mdl(1)
    OUT_min_res_obs_1st(i) = min_res_obs(1)
    CALL com_sprd(1,PARA_DA_MEM,tmp_DA_xf(1,:),OUT_spread_mdl_1st(i))
    CALL com_sprd(1,PARA_DA_MEM,tmp_DA_hx(1,:),OUT_spread_obs_1st(i))
    !--- rank hist ---
    work1 = OUT_xa(1,:)
    CALL qsort(PARA_DA_MEM,work1)
    work2 = DATA_truth(1,i)
    work3 = 0
    DO j = 1,PARA_DA_MEM
      IF (work1(j)<work2) work3 = work3+1
    ENDDO
    work3 = work3 + 1
    OUT_rank_hstg(work3) = OUT_rank_hstg(work3) + 1

    IF (isnan(OUT_RMSE_AN_F(i))) THEN ! Filter Divergence
      OUT_RMSE_FG_F(i:PARA_DA_TIMES) = UNDEF
      OUT_RMSE_AN_F(i:PARA_DA_TIMES) = UNDEF
      OUT_RMSE_AN_S(i:PARA_DA_TIMES) = UNDEF
      OUT_SPRD_FG_F(i:PARA_DA_TIMES) = UNDEF
      OUT_SPRD_AN_F(i:PARA_DA_TIMES) = UNDEF
      OUT_SPRD_AN_S(i:PARA_DA_TIMES) = UNDEF
      EXIT_FLAG = .TRUE.
      FD_FLAG   = .TRUE.
    ENDIF

    !--- forecast ---
    ALLOCATE(tmp_matrix1(PARA_L96_J,tmp_DA_TIMES+1))
    DO k = 1,PARA_DA_MEM
      CALL RUN_L96_NL(OUT_xa(:,k),tmp_DA_TIMES,tmp_matrix1)
      tmp_DA_xf(:,k) = tmp_matrix1(:,tmp_DA_TIMES+1)
    ENDDO
    DEALLOCATE(tmp_matrix1)

    !--- 2.5. Monitor ---
    IF (mod(i,int(4.0d0/PARA_L96_outper)) == 1) THEN
      WRITE(*,*) ' '
      WRITE(*,*) ' ', i, ' / ', PARA_DA_TIMES
      WRITE(*,*) '     RMSE Filter   AN:      ', OUT_RMSE_AN_F(i)
      IF (i>=PARA_DA_spinupstep) THEN
        CALL com_mean   (i-PARA_DA_spinupstep+1,                  OUT_RMSE_AN_F  (PARA_DA_spinupstep:i),                  tmp_var1)
        CALL com_mean   (i-PARA_DA_spinupstep+1,                  OUT_min_res_mdl(PARA_DA_spinupstep:i),                  tmp_var2)
        CALL com_mean   (i-PARA_DA_spinupstep+1,                  OUT_min_res_obs(PARA_DA_spinupstep:i),                  tmp_var3)
        WRITE(*,*) 'MEAN RMSE Filter   AN:    ', tmp_var1
        WRITE(*,*) 'MEAN MINRES_MDL:    ', tmp_var2
        WRITE(*,*) 'MEAN MINRES_OBS:    ', tmp_var3
        RMSE = tmp_var1
        IF (RMSE>3.5d0 .OR. isnan(RMSE)) THEN ! Filter Divergence
       ! IF (OUT_RMSE_AN_F(i)>3.0d0 .OR. isnan(RMSE)) THEN ! Filter Divergence
      !  IF (isnan(RMSE)) THEN ! Filter Divergence
          OUT_RMSE_FG_F(i:PARA_DA_TIMES) = UNDEF
          OUT_RMSE_AN_F(i:PARA_DA_TIMES) = UNDEF
          OUT_RMSE_AN_S(i:PARA_DA_TIMES) = UNDEF
          OUT_SPRD_FG_F(i:PARA_DA_TIMES) = UNDEF
          OUT_SPRD_AN_F(i:PARA_DA_TIMES) = UNDEF
          OUT_SPRD_AN_S(i:PARA_DA_TIMES) = UNDEF
          EXIT_FLAG = .TRUE.
          FD_FLAG   = .TRUE.
        ENDIF
      ENDIF
      WRITE(*,*) '------------------' 
    ENDIF

    ! test
   ! OUT_RMSE_AN_F = 1.55d0
    OUT_min_res_mdl = 0.080d0
    EXIT_FLAG = .TRUE.

  ENDDO
!  write(*,*) OUT_min_res_mdl_1st(1:5)-OUT_min_res_obs_1st(1:5)
!  write(*,*) OUT_SPREAD_MDL_1st(1:5)-OUT_SPREAD_OBS_1st(1:5)

  WRITE(*,*) ''
  WRITE(*,*) '----------------------------------------'
  WRITE(*,*) '   Finish LPF !!'
  WRITE(*,*) 'MEAN RMSE: '  ,RMSE
  WRITE(*,*) 'MEAN MINRES_MDL: '  ,tmp_var2
  WRITE(*,*) 'MEAN MINRES_OBS: '  ,tmp_var3
  WRITE(*,*) 'EXIT_FLAG: ',EXIT_FLAG
  WRITE(*,*) 'FD_FLAG: '  ,FD_FLAG
  WRITE(*,*) '----------------------------------------'


  RETURN
  !=============================================
  ! 3. OUT
  !=============================================
  !--- Setting ---
  OUT_FILE1 =trim(OUT_DIR)//'/'//'OUT_RMSE.txt'
  OUT_FILE2 =trim(OUT_DIR)//'/'//'OUT_SPRD.txt'
  OUT_FILE3 =trim(OUT_DIR)//'/'//'RMSE_FG_F.dat'
  OUT_FILE4 =trim(OUT_DIR)//'/'//'RMSE_AN_F.dat'
  OUT_FILE5 =trim(OUT_DIR)//'/'//'RMSE_AN_S.dat'
  OUT_FILE6 =trim(OUT_DIR)//'/'//'SPRD_FG_F.dat'
  OUT_FILE7 =trim(OUT_DIR)//'/'//'SPRD_AN_F.dat'
  OUT_FILE8 =trim(OUT_DIR)//'/'//'SPRD_AN_S.dat'
  OUT_FILE9 =trim(OUT_DIR)//'/'//'OUT_min_res.txt'
  OUT_FILE10=trim(OUT_DIR)//'/'//'OUT_min_res_mdl.dat'
  OUT_FILE11=trim(OUT_DIR)//'/'//'OUT_min_res_obs.dat'
  OUT_FILE12=trim(OUT_DIR)//'/'//'OUT_rankhist.txt'
  OUT_FILE13=trim(OUT_DIR)//'/'//'OUT_MIN_RES_MDL_1st.dat'
  OUT_FILE14=trim(OUT_DIR)//'/'//'OUT_MIN_RES_OBS_1st.dat'
  OUT_FILE15=trim(OUT_DIR)//'/'//'OUT_SPREAD_MDL_1st.dat'
  OUT_FILE16=trim(OUT_DIR)//'/'//'OUT_SPREAD_OBS_1st.dat'
  write(*,*) OUT_FILE1
  !--- Text ---
  OPEN(71,FILE=OUT_FILE1, STATUS='replace')
  OPEN(72,FILE=OUT_FILE2, STATUS='replace')
  OPEN(73,FILE=OUT_FILE9, STATUS='replace')
  OPEN(74,FILE=OUT_FILE12,STATUS='replace')
  DO i = 1,PARA_DA_TIMES
    WRITE(71,'(f10.6,f10.6,f10.6)') OUT_RMSE_FG_F(i),OUT_RMSE_AN_F(i),OUT_RMSE_AN_S(i)
    WRITE(72,'(f10.6,f10.6,f10.6)') OUT_SPRD_FG_F(i),OUT_SPRD_AN_F(i),OUT_SPRD_AN_S(i)
    WRITE(73,'(f10.6,f10.6)') OUT_min_res_mdl(i),OUT_min_res_obs(i)
  ENDDO
  DO i = 1,PARA_DA_MEM+1
    WRITE(74,'(I5)') OUT_rank_hstg(i)
  ENDDO
  CLOSE(71)
  CLOSE(72)
  CLOSE(73)
  CLOSE(74)
  !--- Binary ---
  OPEN(71,FILE=OUT_FILE3, FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(72,FILE=OUT_FILE4, FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(73,FILE=OUT_FILE5, FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(74,FILE=OUT_FILE6, FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(75,FILE=OUT_FILE7, FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(76,FILE=OUT_FILE8, FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(77,FILE=OUT_FILE10,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(78,FILE=OUT_FILE11,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(79,FILE=OUT_FILE13,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(80,FILE=OUT_FILE14,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(81,FILE=OUT_FILE15,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(82,FILE=OUT_FILE16,FORM='unformatted',ACTION='write',STATUS='replace')
  WRITE(71) REAL(OUT_RMSE_FG_F)
  WRITE(72) REAL(OUT_RMSE_AN_F)
  WRITE(73) REAL(OUT_RMSE_AN_S)
  WRITE(74) REAL(OUT_SPRD_FG_F)
  WRITE(75) REAL(OUT_SPRD_AN_F)
  WRITE(76) REAL(OUT_SPRD_AN_S)
  WRITE(77) REAL(OUT_min_res_mdl)
  WRITE(78) REAL(OUT_min_res_obs)
  WRITE(79) REAL(OUT_min_res_mdl_1st)
  WRITE(80) REAL(OUT_min_res_obs_1st)
  WRITE(81) REAL(OUT_spread_mdl_1st)
  WRITE(82) REAL(OUT_spread_obs_1st)
  CLOSE(71)
  CLOSE(72)
  CLOSE(73)
  CLOSE(74)
  CLOSE(75)
  CLOSE(76)
  CLOSE(77)
  CLOSE(78)
  CLOSE(79)
  CLOSE(80)
  CLOSE(81)
  CLOSE(82)

  RETURN
END SUBROUTINE LPF_cycle


END MODULE LPF
