PROGRAM L96DA_main
!$ use omp_lib

!======================================================================================
  USE common
  USE common_L96
  USE common_L96DA
  USE LPF

  IMPLICIT NONE

  !--- Truth ---
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: Initial 
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: NatureRun_ALL,NatureRun
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: Truth
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: DATA_truth
  INTEGER :: NTIMES_1yr 

  !--- Obs ---
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: rand_num
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: OBS_org,OBS_var
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: OBS_org_mix1, OBS_org_mix2
  INTEGER     ,ALLOCATABLE,DIMENSION(:,:) :: OBS_loc
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: DATA_obs_org,DATA_obs_var
  INTEGER     ,ALLOCATABLE,DIMENSION(:,:) :: DATA_obs_loc

  !-- Initial ---
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: tmp_matrix1,tmp_matrix2
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_matrix3
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: DATA_initial

  !-- OUT ---
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: RMSE_FG_F,RMSE_AN_F,RMSE_AN_S
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: SPRD_FG_F,SPRD_AN_F,SPRD_AN_S
  CHARACTER(256) :: OUT_FILE1,OUT_FILE2,OUT_FILE3,OUT_FILE4
  CHARACTER(256) :: OUT_FILE5,OUT_FILE6,OUT_FILE7,OUT_FILE8,OUT_FILE9
  CHARACTER(256) :: TEXT1,TEXT2,TEXT3,TEXT4,TEXT5,TEXT6,TEXT7,TEXT8,TEXT9
  CHARACTER(256) :: TEXT99
  CHARACTER(256) :: TRUTH_FILE,OBS_FILE1,OBS_FILE2,OBS_FILE3,INI_FILE
  
  !--- Filter Divergence ---
  LOGICAL :: FD_FLAG  
  INTEGER :: FD_COUNTER 

  INTEGER,PARAMETER :: LOOP_NUM = 20
  INTEGER :: i,j,k,iii
  INTEGER :: irec1,irec2,irec3
  REAL :: time1

  !--- OUT ---
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: OUT_RMSE, OUT_SPRD
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: OUT_minres_mdl, OUT_minres_obs
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: OUT_minres_mdl_1st, OUT_minres_obs_1st
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: OUT_spread_mdl_1st, OUT_spread_obs_1st
  INTEGER     ,ALLOCATABLE,DIMENSION(:,:) :: OUT_rank_histg

!======================================================================================
  !-----------------------------------------------------------------------
  ! Read namelists 
  !-----------------------------------------------------------------------
  CALL set_common_L96DA
  CALL monit_para

  write(TEXT1,'(f4.2)') PARA_L96_outper
  write(TEXT2,'(i0)')   PARA_L96_calyear
  write(TEXT3,'(i0)')   PARA_OBS_num
  write(TEXT4,'(f4.2)') PARA_obs_avail_hr
  write(TEXT7,'(i0)')   PARA_DA_MEM
  IF (PARA_obs_nonlinear==0) THEN ! linear
    write(TEXT6,'(a)') '/linear/'
  ELSEIF (PARA_obs_nonlinear==1) THEN ! nonlinear, Y=H*x^2
    write(TEXT6,'(a)') '/nonlinear1/'
  ELSEIF (PARA_obs_nonlinear==2) THEN ! nonlinear, Y=H*log(abs(x))
    write(TEXT6,'(a)') '/nonlinear2/'
  ELSEIF (PARA_obs_nonlinear==91) THEN 
    write(TEXT6,'(a)') '/mix1/'
  ELSEIF (PARA_obs_nonlinear==92) THEN 
    write(TEXT6,'(a)') '/mix2/'
  ELSEIF (PARA_obs_nonlinear==93) THEN 
    write(TEXT6,'(a)') '/mix3/'
  ENDIF
  TEXT99 = trim(TEXT6)//'/obsnum'//trim(TEXT3)//'_obshr'//trim(TEXT4)
    
  TRUTH_FILE=trim(TRUTH_DIR)//'/L96outper'//trim(TEXT1)//'_L96calyr'//trim(TEXT2) //'_TRUTH.dat'
  OBS_FILE2 =trim(OBS_DIR)//trim(TEXT99)//'_OBSloc.dat'
  OBS_FILE3 =trim(OBS_DIR)//trim(TEXT99)//'_OBSvar.dat'
  INI_FILE  =trim(INI_DIR)//'/ini_'//trim(TEXT7)//'mem.dat'
  irec1 = 8*PARA_L96_J*(PARA_L96_NTIMES-NTIMES_1yr)
  irec2 = 8*PARA_obs_num*(PARA_L96_NTIMES-NTIMES_1yr)
  irec3 = 8*PARA_L96_J*PARA_DA_MEM
  write(*,*) OBS_FILE3
  write(*,*) OUT_DIR

  !-----------------------------------------------------------------------
  ! Truth 
  !-----------------------------------------------------------------------
  ALLOCATE(Truth(PARA_L96_J,PARA_L96_NTIMES-NTIMES_1yr)); Truth = UNDEF
  IF (PARA_truth_create) THEN
    !- Run L96 -
    ALLOCATE(Initial      (PARA_L96_J));                                   Initial       = 0.0d0
    ALLOCATE(NatureRun_ALL(PARA_L96_J,PARA_L96_NTIMES+1));                 NatureRun_ALL = UNDEF
    ALLOCATE(NatureRun    (PARA_L96_J,(PARA_L96_NTIMES/PARA_L96_NWRT)+1)); NatureRun     = UNDEF
    CALL RUN_L96_NL(Initial,PARA_L96_NTIMES,NatureRun_ALL)
    NatureRun(:,:) = NatureRun_ALL(:,1:PARA_L96_NTIMES+1:PARA_L96_NWRT)
    !- Cut Spinup -
    NTIMES_1yr = int(dble(PARA_L96_oneyear)*PARA_L96_oneday/PARA_L96_outper)
    Truth(:,:) = NatureRun(:,1+NTIMES_1yr+1:)
    !- Save -
    OPEN(71,FILE=TRUTH_FILE,FORM='unformatted',ACCESS='direct',RECL=irec1,ACTION='write',STATUS='replace')
    DO i = 1,PARA_L96_J
      WRITE(71,REC=i) Truth(i,:)
    ENDDO
    CLOSE(71)
    DEALLOCATE(Initial,NatureRun_ALL,NatureRun)
  ELSE
    i = PARA_L96_J*PARA_L96_NTIMES-NTIMES_1yr
    !- Load -
    OPEN(71,FILE=TRUTH_FILE,FORM='unformatted',ACCESS='direct',RECL=irec1,ACTION='read')
    DO i = 1,PARA_L96_J
      READ(71,REC=i) Truth(i,:)
    ENDDO
    CLOSE(71)
  ENDIF
  ALLOCATE(DATA_truth(PARA_L96_J,PARA_DA_TIMES)); DATA_truth = UNDEF
  DATA_truth(:,:) = Truth(:,1:PARA_DA_TIMES)
  !- Error Trap -
  IF ( minval(DATA_truth) < -999 ) THEN
    WRITE(*,*) "ERROR IN TRUTH !!"
    STOP
  ENDIF

  CALL CPU_TIME(time1)
  write(*,*) ' '
  write(*,*) ' FINISH TRUTH !! '
  write(*,*) time1, ' sec'
  write(*,*) '==============================================='

  !-----------------------------------------------------------------------
  ! Observation 
  !-----------------------------------------------------------------------
  ALLOCATE(OBS_org      (PARA_L96_J  ,PARA_L96_NTIMES-NTIMES_1yr)); OBS_org  = UNDEF
  ALLOCATE(OBS_org_mix1 (PARA_L96_J  ,PARA_L96_NTIMES-NTIMES_1yr)); OBS_org_mix1 = UNDEF
  ALLOCATE(OBS_org_mix2 (PARA_L96_J  ,PARA_L96_NTIMES-NTIMES_1yr)); OBS_org_mix2 = UNDEF
  ALLOCATE(OBS_loc (PARA_obs_num,PARA_L96_NTIMES-NTIMES_1yr)); OBS_loc  = -9999
  ALLOCATE(OBS_var (PARA_obs_num,PARA_L96_NTIMES-NTIMES_1yr)); OBS_var  = UNDEF
  IF (PARA_obs_create) THEN
    ALLOCATE(rand_num(PARA_L96_J  ,PARA_L96_NTIMES-NTIMES_1yr)); rand_num = UNDEF
    !- create noise -
    DO i = 1,PARA_L96_NTIMES-NTIMES_1yr
      CALL com_randn(PARA_L96_J,rand_num(:,i))
    ENDDO
    !- create observation by adding the noise -
    SELECT CASE (PARA_obs_nonlinear)
      CASE (0) ! linear
      rand_num = rand_num*PARA_DA_R_std
      OBS_org = Truth + rand_num + PARA_obs_bias
      CASE (1) ! nonlinear, Y=H*x^2
      rand_num = rand_num*PARA_DA_R_std
      OBS_org = Truth**2 + rand_num + PARA_obs_bias
      CASE (2) ! nonlinear, Y=H*log(abs(x))
      rand_num = rand_num*PARA_DA_R_std
      OBS_org = log(abs(Truth + rand_num)) + PARA_obs_bias
      CASE (91) ! mix1 (Y=H,Y=Hx^2)
      rand_num = rand_num*1.0d0
      OBS_org_mix1 = Truth + rand_num + PARA_obs_bias
      rand_num = rand_num*1.0d0
      OBS_org_mix2 = Truth**2 + rand_num + PARA_obs_bias
     ! OBS_org_mix2 = Truth + rand_num + PARA_obs_bias
      CASE (92) ! mix2 (Y=H,Y=H*log(abs(x))
      rand_num = rand_num*1.0d0
      OBS_org_mix1 = Truth + rand_num + PARA_obs_bias
      rand_num = rand_num*0.10d0
      OBS_org_mix2 = log(abs(Truth + rand_num)) + PARA_obs_bias
      CASE (93) ! mix3 (Y=Hx^2,Y=H*log(abs(x))
      rand_num = rand_num*1.0d0
      OBS_org_mix1 = Truth**2 + rand_num + PARA_obs_bias
      rand_num = rand_num*0.10d0
      OBS_org_mix2 = log(abs(Truth + rand_num)) + PARA_obs_bias
    ENDSELECT
    !- observation operator -
    DO i = 1,PARA_L96_NTIMES-NTIMES_1yr
      
      IF (PARA_obs_nonlinear > 90) THEN
        IF (PARA_obs_num/=20) THEN
          WRITE(*,*) 'ERROR in OBS!!'
          PRINT *
          READ *
        ENDIF
       ! OBS_loc(:,i) = [1,2,3,4,5,6,7,8,9,10,21,22,23,24,25,26,27,28,29,30]
        OBS_loc(:,i) = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]
        DO j = 1,10
          OBS_var(j,i) = OBS_org_mix1(OBS_loc(j,i),i)
        ENDDO
        DO j = 11,20
          OBS_var(j,i) = OBS_org_mix2(OBS_loc(j,i),i)
        ENDDO
      ELSE
        IF (PARA_obs_stable .AND. PARA_obs_num==20 )THEN
          OBS_loc(:,i) = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39]
       !   OBS_loc(:,i) = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
       ! OBS_loc(:,i) = [1,2,3,4,5,11,12,13,14,15,21,22,23,24,35,31,32,33,34,35]
       ! OBS_loc(:,i) = [1,2,3,4,5,6,7,8,9,10,21,22,23,34,35,36,37,38,39,40]
       ! OBS_loc(:,i) = [1,2,3,4,5,6,7,8,9,10,11,12,23,24,25,26,27,28,29,30]
        ELSEIF (PARA_obs_stable .AND. PARA_obs_num==15 )THEN
          OBS_loc(:,i) = [1,5,7,9,13,15,17,21,23,25,29,31,33,37,39]
        ELSEIF (PARA_obs_stable .AND. PARA_obs_num==10 )THEN
          OBS_loc(:,i) = [1,5,9,13,17,21,25,29,33,37]
        ELSEIF (PARA_obs_stable .AND. PARA_obs_num==5 )THEN
          OBS_loc(:,i) = [1,9,17,25,33]
        ELSE
          CALL com_randperm(PARA_L96_J,PARA_obs_num,OBS_loc(:,i))
        ENDIF
        OBS_var(:,i) = OBS_org(OBS_loc(:,i),i)
      ENDIF
    ENDDO
    !- Save -
    OPEN(72,FILE=OBS_FILE2,FORM='unformatted',ACCESS='direct',RECL=irec2,ACTION='write',STATUS='replace')
    OPEN(73,FILE=OBS_FILE3,FORM='unformatted',ACCESS='direct',RECL=irec2,ACTION='write',STATUS='replace')
    DO i = 1,PARA_obs_num
      WRITE(72,REC=i) OBS_loc(i,:)
      WRITE(73,REC=i) OBS_var(i,:)
    ENDDO
    CLOSE(72)
    CLOSE(73)
    DEALLOCATE(rand_num)
  ELSE
    !- Load -
    OPEN(72,FILE=OBS_FILE2,FORM='unformatted',ACCESS='direct',RECL=irec2,ACTION='read')
    OPEN(73,FILE=OBS_FILE3,FORM='unformatted',ACCESS='direct',RECL=irec2,ACTION='read')
    DO i = 1,PARA_obs_num
      READ(72,REC=i) OBS_loc(i,:)
      READ(73,REC=i) OBS_var(i,:)
    ENDDO
    CLOSE(72)
    CLOSE(73)
  ENDIF

  !- Finalize -
  ALLOCATE(DATA_obs_var(PARA_obs_num,PARA_DA_TIMES)); DATA_obs_var = UNDEF
  ALLOCATE(DATA_obs_loc(PARA_obs_num,PARA_DA_TIMES)); DATA_obs_loc = -9999
  DATA_obs_var = OBS_var(:,1:PARA_DA_TIMES)
  DATA_obs_loc = OBS_loc(:,1:PARA_DA_TIMES)
  !- Error Trap -
  IF (( minval(DATA_obs_var) < -999 )  .OR. &
    & ( minval(DATA_obs_loc) < -999 )) THEN
    WRITE(*,*) "ERROR IN OBSERVATION !!"
    write(*,*) minval(DATA_obs_var)
    write(*,*) minval(DATA_obs_loc)
    STOP
  ENDIF

  !- Deallocate -
  DEALLOCATE(Truth)
  DEALLOCATE(OBS_org,OBS_loc,OBS_var)
  DEALLOCATE(OBS_org_mix1,OBS_org_mix2)
  
  CALL CPU_TIME(time1)
  write(*,*) ' '
  write(*,*) ' FINISH OBSERVATION !! '
  write(*,*) time1, ' sec'
  write(*,*) '==============================================='

  !-----------------------------------------------------------------------
  ! Jan.31,2022
  !-----------------------------------------------------------------------
  ALLOCATE(OUT_RMSE(LOOP_NUM,PARA_DA_TIMES)); OUT_RMSE = UNDEF
  ALLOCATE(OUT_SPRD(LOOP_NUM,PARA_DA_TIMES)); OUT_SPRD = UNDEF
  ALLOCATE(OUT_minres_mdl(LOOP_NUM,PARA_DA_TIMES)); OUT_minres_mdl = UNDEF
  ALLOCATE(OUT_minres_obs(LOOP_NUM,PARA_DA_TIMES)); OUT_minres_obs = UNDEF
  ALLOCATE(OUT_minres_mdl_1st(LOOP_NUM,PARA_DA_TIMES)); OUT_minres_mdl_1st = UNDEF
  ALLOCATE(OUT_minres_obs_1st(LOOP_NUM,PARA_DA_TIMES)); OUT_minres_obs_1st = UNDEF
  ALLOCATE(OUT_spread_mdl_1st(LOOP_NUM,PARA_DA_TIMES)); OUT_spread_mdl_1st = UNDEF
  ALLOCATE(OUT_spread_obs_1st(LOOP_NUM,PARA_DA_TIMES)); OUT_spread_obs_1st = UNDEF
  ALLOCATE(OUT_rank_histg(LOOP_NUM,PARA_DA_TIMES)); OUT_rank_histg = -9999
  !$omp parallel do private(iii,FD_COUNTER,DATA_initial,tmp_matrix1,tmp_matrix2,tmp_matrix3,i,time1)
  DO iii = 1,LOOP_NUM 
      write(*,*) iii ! <- don't remove this ???

    !-----------------------------------------------------------------------
    ! RECREATE INITIAL ENSEMBLE WHEN FILTER DIVERGENCE OCCURED
    !-----------------------------------------------------------------------
    FD_COUNTER = 0
    DO WHILE (FD_COUNTER < 1)
      FD_COUNTER = FD_COUNTER + 1
  
      !-----------------------------------------------------------------------
      ! Initial (tekitou)
      !-----------------------------------------------------------------------
      ALLOCATE(DATA_initial(PARA_L96_J,PARA_DA_MEM)); DATA_initial = UNDEF
      ALLOCATE(tmp_matrix1 (PARA_L96_J))            ; tmp_matrix1  = UNDEF
      ALLOCATE(tmp_matrix2 (PARA_L96_J))            ; tmp_matrix2  = UNDEF
      ALLOCATE(tmp_matrix3 (PARA_L96_J,101))        ; tmp_matrix3  = UNDEF
      DO i = 1,PARA_DA_MEM
        CALL com_randn(PARA_L96_J,tmp_matrix1)
        tmp_matrix2 = DATA_truth(:,PARA_DA_TIMES) + tmp_matrix1
        CALL RUN_L96_NL(tmp_matrix2,100,tmp_matrix3)
        DATA_initial(:,i) = tmp_matrix3(:,101)
      ENDDO
      DEALLOCATE(tmp_matrix1,tmp_matrix2,tmp_matrix3)
      
      !- Error Trap -
      IF ( minval(DATA_initial) < -999 ) THEN
        WRITE(*,*) "ERROR IN INITIAL !!"
        STOP
      ENDIF
    
      CALL CPU_TIME(time1)
      write(*,*) ' '
      write(*,*) ' FINISH INITIAL !! '
      write(*,*) time1, ' sec'
      write(*,*) '==============================================='
  
    !-----------------------------------------------------------------------
    ! Data Assimilation
    !-----------------------------------------------------------------------
      SELECT CASE (PARA_DA_MODE)
        CASE ('LPF') 
        CALL LPF_cycle                 (DATA_truth,DATA_obs_var,DATA_obs_loc,DATA_initial, & ! IN
                                        FD_FLAG,                                           & ! OUT
                                        OUT_RMSE(iii,:),OUT_SPRD(iii,:),                   & ! OUT 
                                        OUT_minres_mdl(iii,:),OUT_minres_obs(iii,:),OUT_rank_histg(iii,:), & ! OUT
                                        OUT_minres_mdl_1st(iii,:),OUT_minres_obs_1st(iii,:), &   ! OUT
                                        OUT_spread_mdl_1st(iii,:),OUT_spread_obs_1st(iii,:))     ! OUT
      ENDSELECT
  
      DEALLOCATE(DATA_initial)
      IF (FD_FLAG) THEN
        CONTINUE
      ELSE
        EXIT
      ENDIF

    ENDDO ! while end
  ENDDO ! LOOP_NUM end
  !$omp end parallel do

  !- Binary out -
  OUT_FILE1 =trim(OUT_DIR)//'/'//'RMSE.dat'
  OUT_FILE2 =trim(OUT_DIR)//'/'//'SPRD.dat'
  OUT_FILE3 =trim(OUT_DIR)//'/'//'MINRES_MDL.dat'
  OUT_FILE4 =trim(OUT_DIR)//'/'//'MINRES_OBS.dat'
  OUT_FILE5 =trim(OUT_DIR)//'/'//'RANK_HISTG.dat'
  OUT_FILE6 =trim(OUT_DIR)//'/'//'MINRES_MDL_1st.dat'
  OUT_FILE7 =trim(OUT_DIR)//'/'//'MINRES_OBS_1st.dat'
  OUT_FILE8 =trim(OUT_DIR)//'/'//'SPREAD_MDL_1st.dat'
  OUT_FILE9 =trim(OUT_DIR)//'/'//'SPREAD_OBS_1st.dat'
  OPEN(71,FILE=OUT_FILE1,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(72,FILE=OUT_FILE2,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(73,FILE=OUT_FILE3,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(74,FILE=OUT_FILE4,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(75,FILE=OUT_FILE5,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(76,FILE=OUT_FILE6,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(77,FILE=OUT_FILE7,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(78,FILE=OUT_FILE8,FORM='unformatted',ACTION='write',STATUS='replace')
  OPEN(79,FILE=OUT_FILE9,FORM='unformatted',ACTION='write',STATUS='replace')
  WRITE(71) REAL(OUT_RMSE)
  WRITE(72) REAL(OUT_SPRD)
  WRITE(73) REAL(OUT_minres_mdl)
  WRITE(74) REAL(OUT_minres_obs)
  WRITE(75) INT(OUT_rank_histg)
  WRITE(76) REAL(OUT_minres_mdl_1st)
  WRITE(77) REAL(OUT_minres_obs_1st)
  WRITE(78) REAL(OUT_spread_mdl_1st)
  WRITE(79) REAL(OUT_spread_obs_1st)
  CLOSE(71)
  CLOSE(72)
  CLOSE(73)
  CLOSE(74)
  CLOSE(75)
  CLOSE(76)
  CLOSE(77)
  CLOSE(78)
  CLOSE(79)

  CALL CPU_TIME(time1)
  write(*,*) ' '
  write(*,*) ' FINISH DA !! '
  write(*,*) ' FD_COUNTER: ', FD_COUNTER
  write(*,*) time1, ' sec'
  write(*,*) '==============================================='
 
stop

  



  

END
