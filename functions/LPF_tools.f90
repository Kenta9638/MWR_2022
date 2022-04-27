!
! Jonathan Poterjoy, (June 2020)
!
MODULE LPF_tools
!$use omp_lib

  USE common
  USE common_mtx
  USE common_mpi
  USE common_rankine
  use f95_lapack
  USE Serial_EnSRF_tools
  USE SW_test
  USE kk_netcdf_common
  USE kk_netcdf_tools

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: LPF_CORE, CAL_MIN_RES_SW_TEST

  CONTAINS

!-----------------------------------------
  SUBROUTINE LPF_CORE(org_xf,loc_fact,org_Yo,org_R,                    & ! IN
                      org_H,org_hx,org_HC,org_HCH,                     & ! IN
                      LPF_adaptive_hybrid,test_min_res,minres_dim,     & ! IN
                      xa,                                              & ! OUT
                      minres_sample,minres_hx,minres_loc,mpi_loop_num, & ! IN  20210720
                      out_min_res_mdl,out_min_res_obs)                   ! OUT  
    IMPLICIT NONE
    ! in & out
    REAL(r_size),INTENT(IN)  :: org_xf   (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)  :: loc_fact (PARA_L96_J,PARA_L96_J)
    REAL(r_size),INTENT(IN)  :: org_Yo   (PARA_obs_num)
    REAL(r_size),INTENT(IN)  :: org_R    (PARA_obs_num,PARA_obs_num)
    INTEGER     ,INTENT(IN)  :: org_H    (PARA_obs_num,PARA_L96_J)
    REAL(r_size),INTENT(IN)  :: org_hx   (PARA_obs_num,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)  :: org_HC   (PARA_obs_num,PARA_L96_J)
    REAL(r_size),INTENT(IN)  :: org_HCH  (PARA_obs_num,PARA_obs_num)
    LOGICAL     ,INTENT(IN)  :: LPF_adaptive_hybrid
    REAL(r_size),INTENT(IN)  :: test_min_res
    INTEGER     ,INTENT(IN)  :: minres_dim                                 !    20210720
    REAL(r_size),INTENT(OUT) :: xa       (PARA_L96_J,PARA_DA_MEM)

    REAL(r_size),INTENT(IN) ,OPTIONAL  :: minres_sample (minres_dim  ,PARA_DA_MEM)   !    20210720
    REAL(r_size),INTENT(IN) ,OPTIONAL  :: minres_hx     (PARA_obs_num,PARA_DA_MEM)   !    20210720
    REAL(r_size),INTENT(IN) ,OPTIONAL  :: minres_loc    (minres_dim   ,minres_dim)   !    20210720
    INTEGER     ,INTENT(IN) ,OPTIONAL  :: mpi_loop_num    !    20211121
    REAL(r_size),INTENT(OUT),OPTIONAL  :: out_min_res_mdl(minres_dim), out_min_res_obs(PARA_obs_num)
    INTEGER     ,PARAMETER   :: max_iterations = 3 !99 ! number of allowed iterations
    REAL(r_size)  :: max_res 
    REAL(r_size)  :: beta_max
    INTEGER       :: niter
    LOGICAL       :: qcpass(PARA_obs_num)
    ! LPF 
    INTEGER :: OBS_num
    REAL(r_size) :: x(PARA_L96_J,PARA_DA_MEM), fg
    INTEGER     ,ALLOCATABLE,DIMENSION(:,:) :: H
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: hx, R
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: HC, HCH, wo
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: SAVE_x, SAVE_hx
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_omega_mdl,tmp_omega_obs
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: omega_mdl,omega_obs
    REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: Yo, pf_infl, var_infl
    REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: res_mdl, res_obs
    REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: beta_mdl, beta_obs
    REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: ws_mdl, ws_obs
    REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: xmpf, hxmpf
    REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: var_a_mdl, var_a_obs
    REAL(r_size),DIMENSION(PARA_DA_MEM) :: tmp_d, wop, wt, dum, tmp_w
    INTEGER :: ind(PARA_DA_MEM)
    LOGICAL :: EXIT_flag
    ! Adaptive Hybrid
    REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: min_res_mdl, min_res_obs 
    REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: tmp_min_res_mdl, tmp_min_res_obs 
    REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: tmp_dis, tmp_loc
    REAL(r_size) :: dist_floor  
    REAL(r_size) :: dist_cap    
    REAL(r_size) :: dist_length 
    LOGICAL :: Adaptive_Hybrid 
    ! others
    INTEGER :: i,j,k,l,m,n
    INTEGER :: counter1,counter2,counter3
    INTEGER                                     :: tmp_num1,tmp_num2,tmp_num3,tmp_num4
    REAL(r_size)                                :: tmp_var1,tmp_var2,tmp_var3,tmp_var4
    REAL(r_size),ALLOCATABLE,DIMENSION(:)       :: tmp_d1_mtrx1,tmp_d1_mtrx2,tmp_d1_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: tmp_d2_mtrx1,tmp_d2_mtrx2,tmp_d2_mtrx3,tmp_d2_mtrx4
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: tmp_d3_mtrx1,tmp_d3_mtrx2,tmp_d3_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:,:) :: tmp_d4_mtrx1,tmp_d4_mtrx2,tmp_d4_mtrx3

    CHARACTER(256) :: INP_FILE1
    REAL(r_sngl), DIMENSION(1)       :: dummy_dim1 = UNDEF
    REAL(r_sngl), DIMENSION(1,1)     :: dummy_dim2 = UNDEF
    REAL(r_sngl), DIMENSION(1,1,1)   :: dummy_dim3 = UNDEF
    REAL(r_sngl), DIMENSION(1,1,1,1) :: dummy_dim4 = UNDEF
    REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:) :: tmp_minres

    !----------------------------------------
    ! SETTING 
    !----------------------------------------
    !--- Parameters ---
    max_res    = 1.0d0      ! initial max residual
    beta_max   = abs(UNDEF) ! maximum allowed tempering coefficient
    niter      = 0          ! initial iteration count
    x          = org_xf
 
    !--- Obs-space priors ---
    ALLOCATE(tmp_d1_mtrx1(PARA_obs_num));             tmp_d1_mtrx1 = UNDEF  ! mean_hx
    ALLOCATE(tmp_d1_mtrx2(PARA_obs_num));             tmp_d1_mtrx2 = UNDEF  ! variance
    DO i = 1,PARA_obs_num
      CALL com_mean (PARA_DA_MEM,org_hx(i,:),tmp_d1_mtrx1(i))
      CALL com_stdev(PARA_DA_MEM,org_hx(i,:),tmp_var1)
      tmp_d1_mtrx2(i) = tmp_var1**2 ! variance
    ENDDO

    !------------------------------
    ! Filter
    !------------------------------
      !--- QC ---
      qcpass  = .TRUE.
      OBS_num = PARA_obs_num 
      DO i = 1,PARA_obs_num
        tmp_var1 = org_Yo(i) - tmp_d1_mtrx1(i) ! d
        IF (tmp_var1 > 4.0d0*sqrt(tmp_d1_mtrx2(i)+org_R(i,i))) THEN
          qcpass(i) = .FALSE.
          OBS_num = OBS_num - 1 
        ENDIF
      ENDDO
      ALLOCATE(hx (OBS_num,PARA_DA_MEM))
      ALLOCATE(Yo (OBS_num))
      ALLOCATE(H  (OBS_num,PARA_L96_J))
      ALLOCATE(R  (OBS_num,OBS_num));   R = 0.0d0
      ALLOCATE(HC (OBS_num,PARA_L96_J))
      ALLOCATE(HCH(OBS_num,OBS_num)); HCH = UNDEF
      counter1 = 0
      DO i = 1,PARA_obs_num
        IF (qcpass(i)) THEN
          counter1 = counter1 + 1
          hx(counter1,:)        = org_hx(i,:)
          Yo(counter1)          = org_Yo(i)
          H (counter1,:)        = org_H (i,:)
          HC(counter1,:)        = org_HC(i,:)
          R (counter1,counter1) = org_R (i,i)
        ENDIF
      ENDDO
      IF (counter1/=OBS_num) WRITE(*,*) 'ERROR 1 !!'
      !--- HC & HCH ---
    !  ALLOCATE(HC (OBS_num,PARA_L96_J))
    !  ALLOCATE(HCH(OBS_num,OBS_num))
    !  HC  = matmul(dble(H),loc_fact)
      HCH = matmul(HC,transpose(dble(H)))
      !--- Deallocate ---
      DEALLOCATE(tmp_d1_mtrx1,tmp_d1_mtrx2)
      !--- min_res ---
      ALLOCATE(min_res_mdl(PARA_L96_J)); min_res_mdl = test_min_res
      ALLOCATE(min_res_obs(OBS_num));    min_res_obs = test_min_res

      !-----------------------------------
      ! Adaptive Hybrid
      !-----------------------------------
      IF (LPF_adaptive_hybrid) THEN 
       
        out_min_res_mdl = UNDEF
        out_min_res_obs = UNDEF

        
        !-- temporary 2021.1.21
        IF (PARA_L96_J==minres_dim) THEN
          CALL CAL_MIN_RES_SW_TEST(PARA_L96_J,PARA_MEM,x, loc_fact, min_res_mdl) ! model space
          CALL CAL_MIN_RES_SW_TEST(OBS_num,   PARA_MEM,hx, HCH,     min_res_obs) ! obs space
        ELSEIF (PARA_L96_J/2==minres_dim) THEN
          ALLOCATE(tmp_min_res_mdl(minres_dim)); tmp_min_res_mdl = test_min_res
          ALLOCATE(tmp_min_res_obs(OBS_num));    tmp_min_res_obs = test_min_res

          !--- model space
          IF (FLAG_READ_MINRES) THEN
            INP_FILE1 = trim(MINRES_NCDIR)//'/'//'OUT_minres.nc'
            ALLOCATE(tmp_minres(91,91,100))
            CALL kk_netcdf_get_var_single(trim(INP_FILE1),'min_res_mdl',3,91,91,100,9999, &
                                          dummy_dim1,dummy_dim2,tmp_minres,dummy_dim4)
            counter1 = 0
            DO i = 1,91
              DO j = 1,91
                counter1 = counter1 + 1
                tmp_min_res_mdl(counter1) = dble(tmp_minres(j,i,mpi_loop_num)) ! <-- i,j !!
              ENDDO
            ENDDO
            IF (minres_dim /= counter1) THEN
              write(*,*) 'ERROR in FLAG_READ_MINRES !!'
              print *
              read *
            ENDIF
            DEALLOCATE(tmp_minres)
          ELSE
            CALL CAL_MIN_RES_SW_TEST(minres_dim,PARA_MEM,minres_sample, minres_loc, tmp_min_res_mdl) ! model space
          ENDIF ! FLAG_READ_MINRES

          min_res_mdl(1:minres_dim)  = tmp_min_res_mdl
          min_res_mdl(minres_dim+1:) = tmp_min_res_mdl
          !--- obs space
          CALL CAL_MIN_RES_SW_TEST(OBS_num,   PARA_MEM,minres_hx,     HCH,        tmp_min_res_obs) ! obs space
          min_res_obs = tmp_min_res_obs
        ELSE
          WRITE(*,*) 'ERROR (Adaptive Hybrid) !!'
          PRINT *
          READ *
        ENDIF

        !--- convert ---
        CALL CONV_MINRES(PARA_L96_J,min_res_mdl,min_res_mdl)
        CALL CONV_MINRES(OBS_num   ,min_res_obs,min_res_obs)

        out_min_res_mdl = min_res_mdl
        out_min_res_obs = min_res_obs
      ENDIF ! adaptive_hybrid end 

      !--- out ---
    !  CALL com_mean (PARA_L96_J,min_res_mdl,mean_min_res_mdl)
    !  CALL com_mean (OBS_num   ,min_res_obs,mean_min_res_obs)

      !--- Initialize variables for iterations ---
      ALLOCATE(res_mdl (PARA_L96_J)); res_mdl  = 1.0d0 ! initial obs-space   residuals
      ALLOCATE(res_obs (OBS_num));    res_obs  = 1.0d0 ! initial model-space residuals
      ALLOCATE(pf_infl (OBS_num));    pf_infl  = UNDEF ! obs inflation factors
      ALLOCATE(beta_mdl(PARA_L96_J)); beta_mdl = UNDEF
      ALLOCATE(beta_obs(OBS_num));    beta_obs = UNDEF


      !--- Include min_res in residual ---
      res_mdl = res_mdl - min_res_mdl!KLD2!min_res
      res_obs = res_obs - min_res_obs!matmul(dble(H),KLD2)!min_res
      max_res = max(maxval(res_mdl),maxval(res_obs))

      !----------------------------------------
      ! MAIN LOOP
      !----------------------------------------
      DO WHILE (max_res > 0) 
        niter = niter + 1

        !--- Allocate ---
        ALLOCATE(SAVE_x (PARA_L96_J,PARA_DA_MEM))
        ALLOCATE(SAVE_hx(OBS_num   ,PARA_DA_MEM))
        ALLOCATE(wo     (OBS_num   ,PARA_DA_MEM))

        !--- Store copy of original particles ---
        SAVE_x  = x
        SAVE_hx = hx

        !--- Obs space weights ---
        ALLOCATE(tmp_d1_mtrx1(PARA_DA_MEM)); tmp_d1_mtrx1 = UNDEF  ! d
        DO i = 1,OBS_num
          ! Squared innovations
          tmp_d1_mtrx1 = (Yo(i)-SAVE_hx(i,:))**2/(2.0d0*R(i,i)) ! d
          tmp_d1_mtrx1 = tmp_d1_mtrx1-minval(tmp_d1_mtrx1)

          ! Obse space weights need to be capped
          fg = max(1.0d0,maxval(tmp_d1_mtrx1)/200.0d0)
          CALL find_obs_infl(tmp_d1_mtrx1/fg,pf_infl(i))
          pf_infl(i) = pf_infl(i)*fg

          ! Weight calculations
          tmp_d1_mtrx1 = tmp_d1_mtrx1/pf_infl(i)
       !   tmp_d1_mtrx1 = tmp_d1_mtrx1-minval(tmp_d1_mtrx1)
          wo(i,:)      = exp(-tmp_d1_mtrx1)
          wo(i,:)      = wo(i,:)/sum(wo(i,:))
        ENDDO
        DEALLOCATE(tmp_d1_mtrx1)

        !--- Calculate regularization coefficients ---
        beta_mdl = UNDEF
        beta_obs = UNDEF
        CALL get_reg(PARA_L96_J,OBS_num,HC ,wo,res_mdl,beta_max,beta_mdl)
        CALL get_reg(OBS_num   ,OBS_num,HCH,wo,res_obs,beta_max,beta_obs)

        !--- Initialize weighting matrices ---
        ALLOCATE(tmp_omega_mdl(PARA_L96_J,PARA_DA_MEM)); tmp_omega_mdl = 0.0d0
        ALLOCATE(tmp_omega_obs(OBS_num   ,PARA_DA_MEM)); tmp_omega_obs = 0.0d0
        ALLOCATE(omega_mdl    (PARA_L96_J,PARA_DA_MEM)); omega_mdl     = 1.0d0
        ALLOCATE(omega_obs    (OBS_num   ,PARA_DA_MEM)); omega_obs     = 1.0d0

        !----------------------------------------
        ! Observation LOOP
        !----------------------------------------
        DO i = 1,OBS_num

          IF (1.0d0 > 0.98d0*dble(PARA_DA_MEM)*sum(wo(i,:)**2) ) CYCLE

          !--- Initialize ---
          tmp_d = 0.0d0
!          wop   = 1.0d0
          wt    = dble(PARA_DA_MEM)*wo(i,:)-1.0d0
          !--- Model-space localized weighting vectors ---
          DO j = 1,PARA_L96_J
            IF (beta_mdl(j) /= beta_max) THEN
              IF (HC(i,j) == 1) THEN
                dum = log(dble(PARA_DA_MEM)*wo(i,:))
              ELSE
                dum = wt*HC(i,j)
                DO k = 1,PARA_DA_MEM
                  IF (abs(dum(k))>0.30d0) dum(k) = log(dum(k)+1)
                ENDDO
              ENDIF
            ENDIF
            tmp_omega_mdl(j,:) = tmp_omega_mdl(j,:) - dum
            tmp_omega_mdl(j,:) = tmp_omega_mdl(j,:) - minval(tmp_omega_mdl(j,:))
          ENDDO
          !--- Obs-space localized weighting vectors ---
          DO j = 1,OBS_num
            IF (beta_obs(j) /= beta_max) THEN
              IF (HCH(i,j) == 1) THEN
                dum = log(dble(PARA_DA_MEM)*wo(i,:))
              ELSE
                dum = wt*HCH(i,j)
                DO k = 1,PARA_DA_MEM
                  IF (abs(dum(k))>0.30d0) dum(k) = log(dum(k)+1)
                ENDDO
              ENDIF
            ENDIF
            tmp_omega_obs(j,:) = tmp_omega_obs(j,:) - dum
            tmp_omega_obs(j,:) = tmp_omega_obs(j,:) - minval(tmp_omega_obs(j,:))
          ENDDO
   !!       !--- Innovations from recently updated particles in sequence ---
   !!       tmp_d = tmp_d + (Yo(i)-hx(i,:))**2/(2*R(i,i)*pf_infl(i)*beta_obs(i))
   !!       !--- Product of obs-space weights needed for ws ---
 ! !!        DO j = 1,PARA_DA_MEM
 ! !!          wop(j) = wop(j)*wo(i,j)
 ! !!        ENDDO
   !!       !--- Additional Inflation ---
   !!       ! Additional inflation factor for sampled particles to cope with
   !!       ! sampling errors that accumulate through each serial assimilation
   !!       ! step.
   !!  !     tmp_d = tmp_d/tmp_infl
   !!       !--- Sampling weights ---
   !!       tmp_w = exp(tmp_d*(-1.0d0))
   !!  !     tmp_w = exp(tmp_omega_obs(i,:)*(-1.0d0)/beta_obs(i))
   !!       IF (sum(tmp_w)==0.0d0) THEN
   !!         tmp_num1 = minval(minloc(tmp_d))
   !!         tmp_w(tmp_num1) = 1.0d0
   !!       ENDIF
   !!       tmp_w = tmp_w/sum(tmp_w)
          !--- Normalization for update equations ---
       !   ALLOCATE(ws_mdl(PARA_L96_J)) 
       !   ALLOCATE(ws_obs(OBS_num))
       !   ws_mdl = UNDEF
       !   ws_obs = UNDEF
        !  ws_mdl = 1.0d0
        !  ws_obs = 1.0d0
      !    wop    = wop/sum(wop)
      !    DO j = 1,PARA_L96_J
      !      ws_mdl(j) = dot_product(omega_mdl(j,:),wop)
      !    ENDDO
      !    DO j = 1,OBS_num
      !      ws_obs(j) = dot_product(omega_obs(j,:),wop)
      !    ENDDO
          !--- Normalize weights ---
          DO j = 1,PARA_DA_MEM
            DO k = 1,PARA_L96_J
              omega_mdl(k,j) = exp(tmp_omega_mdl(k,j)*(-1.0d0)/beta_mdl(k))
            ENDDO
            DO k = 1,OBS_num
              omega_obs(k,j) = exp(tmp_omega_obs(k,j)*(-1.0d0)/beta_obs(k))
            ENDDO
          ENDDO
          ALLOCATE(tmp_d1_mtrx1(PARA_L96_J))
          ALLOCATE(tmp_d1_mtrx2(OBS_num))
          tmp_d1_mtrx1 = sum(omega_mdl,2)
          tmp_d1_mtrx2 = sum(omega_obs,2)

          ! Skip update step if few particles removed
          tmp_w = omega_obs(i,:)/tmp_d1_mtrx2(i)
          IF (1.0d0 > 0.98d0*dble(PARA_DA_MEM)*sum(tmp_w**2) ) THEN
            DEALLOCATE(tmp_d1_mtrx1,tmp_d1_mtrx2)
            CYCLE
          ENDIF

          !--- Get posterior mean ---
          ALLOCATE( xmpf(PARA_L96_J));  xmpf = 0.0d0
          ALLOCATE(hxmpf(OBS_num));    hxmpf = 0.0d0
          DO j = 1,PARA_DA_MEM
            DO k = 1,PARA_L96_J
              omega_mdl(k,j) = omega_mdl(k,j)/tmp_d1_mtrx1(k)
              xmpf(k)        =  xmpf(k) + omega_mdl(k,j)*SAVE_x(k,j)
            ENDDO
            DO k = 1,OBS_num
              omega_obs(k,j) = omega_obs(k,j)/tmp_d1_mtrx2(k)
              hxmpf(k)       = hxmpf(k) + omega_obs(k,j)*SAVE_hx(k,j)
            ENDDO
          ENDDO
          DEALLOCATE(tmp_d1_mtrx1,tmp_d1_mtrx2)
          !--- Variance ---
          ! Calculate posterior variance using vector weights and original
          ! prior particles
          ALLOCATE(var_a_mdl(PARA_L96_J)); var_a_mdl = 0.0d0
          ALLOCATE(var_a_obs(OBS_num));    var_a_obs = 0.0d0
          DO j = 1,PARA_DA_MEM
            DO k = 1,PARA_L96_J
              var_a_mdl(k) = var_a_mdl(k) + omega_mdl(k,j)*(SAVE_x (k,j)-xmpf(k))**2
            ENDDO
            DO k = 1,OBS_num
              var_a_obs(k) = var_a_obs(k) + omega_obs(k,j)*(SAVE_hx(k,j)-hxmpf(k))**2
            ENDDO
          ENDDO
          DO k = 1,PARA_L96_J
            var_a_mdl(k) = var_a_mdl(k)/(1-sum(omega_mdl(k,:)**2)) 
          ENDDO
          DO k = 1,OBS_num
            var_a_obs(k) = var_a_obs(k)/(1-sum(omega_obs(k,:)**2))
          ENDDO
          !--- Apply deterministic resampling ---
          ind = 9999
          CALL sampling(SAVE_hx(i,:),omega_obs(i,:),ind,EXIT_flag)
          IF (EXIT_flag) THEN ! (KK 2020.10.19)
            xa = org_xf
            RETURN
          ENDIF

          !--- Merge steps ---
          ALLOCATE(tmp_d2_mtrx1(PARA_L96_J,PARA_DA_MEM)) 
          ALLOCATE(tmp_d2_mtrx2(OBS_num   ,PARA_DA_MEM))
          ALLOCATE(tmp_d2_mtrx3(PARA_L96_J,PARA_DA_MEM)) !  xa
          ALLOCATE(tmp_d2_mtrx4(OBS_num   ,PARA_DA_MEM)) ! hxa
          DO j = 1,PARA_DA_MEM
          !  tmp_d2_mtrx1(:,j) =  x(:,ind(j))
          !  tmp_d2_mtrx2(:,j) = hx(:,ind(j))
            tmp_d2_mtrx1(:,j) = SAVE_x (:,ind(j))
            tmp_d2_mtrx2(:,j) = SAVE_hx(:,ind(j))
          ENDDO
          CALL pf_merge( x,tmp_d2_mtrx1,HC (i,:), xmpf,var_a_mdl,PARA_L96_J,tmp_d2_mtrx3)
          CALL pf_merge(hx,tmp_d2_mtrx2,HCH(i,:),hxmpf,var_a_obs,OBS_num   ,tmp_d2_mtrx4)
          x  = tmp_d2_mtrx3
          hx = tmp_d2_mtrx4
          DEALLOCATE(tmp_d2_mtrx1,tmp_d2_mtrx2,tmp_d2_mtrx3,tmp_d2_mtrx4)
          !--- Deallocate
!          DEALLOCATE(ws_mdl,ws_obs)
          DEALLOCATE(xmpf,hxmpf)
          DEALLOCATE(var_a_mdl,var_a_obs)

        ENDDO ! obs loop end

        !--- Deallocate
        DEALLOCATE(SAVE_x,SAVE_hx)
        DEALLOCATE(wo)
        DEALLOCATE(tmp_omega_mdl,tmp_omega_obs)
        DEALLOCATE(omega_mdl,omega_obs)

        !--- Break --- 
        IF (niter == max_iterations) THEN
          EXIT
        ENDIF

        !--- Get max of residuals for all states ---
        max_res = max(maxval(res_mdl),maxval(res_obs))
      !  max_res = maxval(res_obs)

      ENDDO ! while end

    xa = x

    !-----------------------------------------------
    ! EnKF step for remaining obs error factor
    !-----------------------------------------------
    IF (maxval(min_res_obs)>0.0d0) THEN
      DO i = 1,OBS_num
        IF (min_res_obs(i) < 0.001d0) min_res_obs(i) = 0.001d0
      ENDDO
      ALLOCATE(var_infl(OBS_num)); var_infl = 1.0d0
      var_infl = var_infl/min_res_obs
      CALL Serial_EnSRF_CORE_tempered(x,loc_fact,OBS_num,H,hx,Yo,R,var_infl,xa)          
      CALL infl_RTPP(x,0.3d0,xa) 
    ENDIF

    RETURN
  END SUBROUTINE LPF_CORE

  SUBROUTINE find_beta(sum_exp,Neff,beta)
    IMPLICIT NONE
    ! in & out
    REAL(r_size),INTENT(IN)  :: sum_exp(PARA_DA_MEM)
    REAL(r_size),INTENT(IN)  :: Neff
    REAL(r_size),INTENT(OUT) :: beta
    REAL(r_size),PARAMETER   :: tol = 1e-5 ! Apply bisection method to find k
    REAL(r_size) :: beta_max
    REAL(r_size) :: w(PARA_DA_MEM)
    REAL(r_size) :: Neff_init, Nf_new
    REAL(r_size) :: var_k (3) ! start, middle, end
    REAL(r_size) :: var_fk(3) ! start, middle, end

    INTEGER :: i,j,k,l,m,n
    INTEGER :: counter1,counter2,counter3
    INTEGER                                     :: tmp_num1,tmp_num2,tmp_num3,tmp_num4
    REAL(r_size)                                :: tmp_var1,tmp_var2,tmp_var3,tmp_var4
    REAL(r_size),ALLOCATABLE,DIMENSION(:)       :: tmp_d1_mtrx1,tmp_d1_mtrx2,tmp_d1_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: tmp_d2_mtrx1,tmp_d2_mtrx2,tmp_d2_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: tmp_d3_mtrx1,tmp_d3_mtrx2,tmp_d3_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:,:) :: tmp_d4_mtrx1,tmp_d4_mtrx2,tmp_d4_mtrx3


    IF (Neff==1) THEN
      beta = 1.0d0
 !     write(*,*) 'find_beta QC 1 !!'
      RETURN
    ENDIF
   
  !  beta_max = maxval(sum_exp)**2
    beta_max = max(1e5,maxval(sum_exp)**4)
    w        = exp(-sum_exp)
    tmp_var1 = sum(w)
    IF (sum(w) > 0) THEN
      w = w/sum(w)
      Neff_init = 1/sum(w**2)
    ELSE
      Neff_init = 1.0d0
 !     write(*,*) 'find_beta QC 2 !!'
    ENDIF


    IF ((Neff_init < Neff) .OR. tmp_var1 == 0) THEN

      ! Initial start and end bounds
      var_k(1) = 1          ! start
      var_k(3) = beta_max   ! end

      DO i = 1,100
        var_k(2) = (var_k(1)+var_k(3))/2 ! middle

        ! Evaluate function at each point
        ! start
        w = exp(-sum_exp/var_k(1))
        w = w/sum(w)
        var_fk(1) = Neff - 1/sum(w**2)
        ! middle 
        w = exp(-sum_exp/var_k(2))
        w = w/sum(w)
        var_fk(2) = Neff - 1/sum(w**2)
        ! end
        w = exp(-sum_exp/var_k(3))
        w = w/sum(w)
        var_fk(3) = Neff - 1/sum(w**2)

        ! Exit criteria
        IF ((var_k(3)-var_k(1))/2 < tol) THEN
          EXIT
        ENDIF

        IF (var_fk(1)*var_fk(2) > 0) THEN
          var_k(1) = var_k(2)
        ELSE
          var_k(3) = var_k(2)
        ENDIF
      ENDDO

      ! Get beta from k
      beta = var_k(2)
      w = exp(-sum_exp/beta)
      w = w/sum(w)
      Nf_new = 1/sum(w**2)

      IF ((Nf_new <= Neff-1.0d0) .OR. (isnan(Nf_new))) THEN
        beta = beta_max
        write(*,*) 'find_beta QC 3 !!'
      ENDIF

    ELSE
      beta = 1.0d0
    ENDIF

    RETURN
  END SUBROUTINE find_beta 

  SUBROUTINE get_reg(in_num1,Ny,C,hw,res,beta_max,beta)
    IMPLICIT NONE
    ! in & out
    INTEGER     ,INTENT(IN)    :: in_num1, Ny
    REAL(r_size),INTENT(IN)    :: C   (Ny,in_num1)
    REAL(r_size),INTENT(IN)    :: hw  (Ny,PARA_DA_MEM)
    REAL(r_size),INTENT(INOUT) :: res (in_num1)
    REAL(r_size),INTENT(IN)    :: beta_max
    REAL(r_size),INTENT(OUT)   :: beta(in_num1)

    REAL(r_size)  :: wo (PARA_DA_MEM)
    REAL(r_size)  :: dum(PARA_DA_MEM)

    ! others
    INTEGER :: i,j,k,l,m,n
    INTEGER :: counter1,counter2,counter3
    INTEGER                                     :: tmp_num1,tmp_num2,tmp_num3,tmp_num4
    REAL(r_size)                                :: tmp_var1,tmp_var2,tmp_var3,tmp_var4
    REAL(r_size),ALLOCATABLE,DIMENSION(:)       :: tmp_d1_mtrx1,tmp_d1_mtrx2,tmp_d1_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: tmp_d2_mtrx1,tmp_d2_mtrx2,tmp_d2_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: tmp_d3_mtrx1,tmp_d3_mtrx2,tmp_d3_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:,:) :: tmp_d4_mtrx1,tmp_d4_mtrx2,tmp_d4_mtrx3

    DO j = 1,in_num1

      !--- QC ---
      IF (res(j) <= 0) THEN
        beta(j) = beta_max
      !  write(*,*) 'get_reg QC 1 !!'
        CYCLE
      ENDIF

      !--- Obs loop ---
      wo = 0.0d0
      DO i = 1,Ny
        IF (C(i,j)==1) THEN
          dum = log(PARA_DA_MEM*hw(i,:))
        ELSE
          dum = (PARA_DA_MEM*hw(i,:)-1)*C(i,j)
          DO k = 1,PARA_DA_MEM
            IF (abs(dum(k))>0.30d0) dum(k) = log(dum(k)+1)
          ENDDO
        ENDIF
        wo  = wo - dum
        wo  = wo - minval(wo)
      ENDDO

      !--- Calculate beta from log of weights ---
      CALL find_beta(wo,PARA_DA_LPF_Neff,beta(j))

      !--- Fix beta if value exceeds residual ---
      IF (res(j) < 1.0d0/beta(j)) THEN
        beta(j) = 1.0d0/res(j)
        res(j)  = 0.0d0
      ELSE
        res(j)  = res(j) - 1.0d0/beta(j)
      ENDIF

      !--- Store residuals ---
      beta(j) = min(beta(j),beta_max)

    ENDDO

    RETURN
  END SUBROUTINE get_reg

  SUBROUTINE sampling(x,w,ind,EXIT_flag)
    IMPLICIT NONE
    ! in & out
    REAL(r_size),INTENT(IN)    :: x  (PARA_DA_MEM)
    REAL(r_size),INTENT(IN)    :: w  (PARA_DA_MEM)
    INTEGER     ,INTENT(OUT)   :: ind(PARA_DA_MEM)
    LOGICAL     ,INTENT(OUT)   :: EXIT_flag
    INTEGER     ,DIMENSION(PARA_DA_MEM)   :: tmp_ind, tmp_tmp_ind, tmp_nums
    REAL(r_size),DIMENSION(PARA_DA_MEM)   :: tmp_x, tmp_w
    REAL(r_size),DIMENSION(PARA_DA_MEM+1) :: cum_weight
    REAL(r_size) :: base,frac
    LOGICAL :: flag
    ! others
    INTEGER :: i,j,k,l,m,n
    INTEGER :: counter1,counter2,counter3
    INTEGER                                     :: tmp_num1,tmp_num2,tmp_num3,tmp_num4
    REAL(r_size)                                :: tmp_var1,tmp_var2,tmp_var3,tmp_var4
    REAL(r_size),ALLOCATABLE,DIMENSION(:)       :: tmp_d1_mtrx1,tmp_d1_mtrx2,tmp_d1_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: tmp_d2_mtrx1,tmp_d2_mtrx2,tmp_d2_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: tmp_d3_mtrx1,tmp_d3_mtrx2,tmp_d3_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:,:) :: tmp_d4_mtrx1,tmp_d4_mtrx2,tmp_d4_mtrx3
   
    !--- Initialize ---
    ind         = 9999
    tmp_ind     = 9999
    tmp_tmp_ind = 9999
    tmp_x       = x
    cum_weight  = 0.0d0
    EXIT_flag = .FALSE.

    tmp_w = w
    tmp_w = tmp_w**0.5d0
    tmp_w = tmp_w/sum(tmp_w)

    !--- Sort ---
    ALLOCATE(tmp_d1_mtrx1(PARA_DA_MEM)); tmp_d1_mtrx1 = UNDEF
    tmp_nums = 9999
    DO i = 1,PARA_DA_MEM
      tmp_nums(i)         = minval(minloc(tmp_x))
      tmp_d1_mtrx1(i)     = minval(tmp_x)
      tmp_x(tmp_nums(i))  = 9999
    ENDDO
    DEALLOCATE(tmp_d1_mtrx1)

    !--- Apply deterministic sampling by taking value at every 1/Ne quantile ---
    ALLOCATE(tmp_d1_mtrx1(PARA_DA_MEM)); tmp_d1_mtrx1 = UNDEF
    DO i = 1,PARA_DA_MEM
      tmp_d1_mtrx1(i) = tmp_w(tmp_nums(i))
      cum_weight(i+1) = sum(tmp_d1_mtrx1(1:i))
    ENDDO
    DEALLOCATE(tmp_d1_mtrx1)

    base = 1/dble(PARA_DA_MEM)/2
    k    = 2
    DO i = 1,PARA_DA_MEM
      frac = base + dble(i-1)/dble(PARA_DA_MEM)
      flag = .FALSE.
      DO WHILE (flag .eqv. .FALSE.)
        IF ((cum_weight(k-1) < frac) .AND. (frac <= cum_weight(k))) THEN
          tmp_tmp_ind(i) = k-1
          flag = .TRUE.
        ELSE
          k = k+1
          IF (k>PARA_DA_MEM+1) THEN
            k = k-1
            tmp_tmp_ind(i) = k-1
            flag = .TRUE.
            EXIT_flag = .TRUE.
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !--- EXIT (KK 2020.10.19) ---
    IF (EXIT_flag) THEN
      ind = -9999
    ENDIF

    !--- Unsort indice ---
    DO i = 1,PARA_DA_MEM
      tmp_ind(i) = tmp_nums(tmp_tmp_ind(i))
    ENDDO

    !--- Replace removed particles with duplicated particles ---
    ! This part is neccesary for term3 in pf_merge.m !! KK
    DO i = 1,PARA_DA_MEM
      flag = .FALSE.
      DO j = 1,PARA_DA_MEM
        IF (tmp_ind(j)==i) THEN
          tmp_num1 = j
          flag = .TRUE.
          EXIT
        ENDIF
      ENDDO
      IF (flag) THEN
        ind(i) = i
        tmp_ind(tmp_num1) = 9999
      ENDIF
    ENDDO
    DO i = 1,PARA_DA_MEM
      IF (tmp_ind(i) /= 9999) THEN
        DO j = 1,PARA_DA_MEM
          IF (ind(j) == 9999) THEN
            ind(j) = tmp_ind(i)
            tmp_ind(i) = 9999
            EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE sampling

  SUBROUTINE pf_merge(xf,xs,loc,xmpf,var_a,in_num,xa)
    IMPLICIT NONE
    ! in & out
    REAL(r_size),INTENT(IN)    :: xf   (in_num,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)    :: xs   (in_num,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)    :: loc  (in_num)
    REAL(r_size),INTENT(IN)    :: xmpf (in_num)
    REAL(r_size),INTENT(IN)    :: var_a(in_num)
    INTEGER     ,INTENT(IN)    :: in_num
    REAL(r_size),INTENT(OUT)   :: xa   (in_num,PARA_DA_MEM)
    REAL(r_size) :: alpha
    REAL(r_size),DIMENSION(in_num) :: c, c2
    REAL(r_size),DIMENSION(in_num) :: term1, term2, term3
    REAL(r_size),DIMENSION(in_num) :: r1, r2
    REAL(r_size),DIMENSION(in_num) :: tmp_A, tmp_B, tmp_C
    REAL(r_size),DIMENSION(in_num) :: alpha2
    REAL(r_size),DIMENSION(in_num) :: pfm
    REAL(r_size),DIMENSION(in_num) :: m1,m2
    ! others
    INTEGER :: i,j,k,l,m,n
    INTEGER :: counter1,counter2,counter3
    INTEGER                                     :: tmp_num1,tmp_num2,tmp_num3,tmp_num4
    REAL(r_size)                                :: tmp_var1,tmp_var2,tmp_var3,tmp_var4
    REAL(r_size),ALLOCATABLE,DIMENSION(:)       :: tmp_d1_mtrx1,tmp_d1_mtrx2,tmp_d1_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: tmp_d2_mtrx1,tmp_d2_mtrx2,tmp_d2_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: tmp_d3_mtrx1,tmp_d3_mtrx2,tmp_d3_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:,:) :: tmp_d4_mtrx1,tmp_d4_mtrx2,tmp_d4_mtrx3
   
    !--- Setting ---
    alpha = PARA_DA_LPF_alpha
    !--- coefficient for merging step ---
  !  c = (1.0d0-loc)/(loc*dble(PARA_DA_MEM)*ws)
    c = (1.0d0-loc)/(loc)
    !--- Calculate r1 and r2 coefficients for weight update equation. ---
    term1 = 0.0d0; term2 = 0.0d0; term3 = 0.0d0
    DO i = 1,PARA_DA_MEM
      term1 = term1 + (xs(:,i)-xmpf)**2
      term2 = term2 + (xf(:,i)-xmpf)**2
      term3 = term3 + (xs(:,i)-xmpf)*(xf(:,i)-xmpf)
    ENDDO
    c2 = c*c
    r1 = term1 + c2*term2 + 2.0d0*c*term3
    r2 = c2/r1
    r1 = alpha*sqrt(dble(PARA_DA_MEM-1)*var_a/r1)
    r2 =       sqrt(dble(PARA_DA_MEM-1)*var_a*r2)
    !--- Calculate alpha2 which is added to r2 to fit posterior variance ---
    DO i = 1,in_num
      CALL com_mean (PARA_DA_MEM,xs(i,:)-xmpf(i),m1(i))
      CALL com_mean (PARA_DA_MEM,xf(i,:)-xmpf(i),m2(i))
    ENDDO

    term1 = term1 - PARA_DA_MEM*m1*m1!dot_product(m1,m1)
    term2 = term2 - PARA_DA_MEM*m2*m2!dot_product(m2,m2)
    term3 = term3 - PARA_DA_MEM*m1*m2!dot_product(m1,m2)
    tmp_A = term2
    tmp_B = 2.0d0*(r1*term3 + r2*term2)
    tmp_C = term1*r1**2 + term2*r2**2 + 2.0d0*term3*r1*r2 - dble(PARA_DA_MEM-1)*var_a
    alpha2 = (-tmp_B + sqrt(tmp_B**2-4.0d0*tmp_A*tmp_C))/(2.0d0*tmp_A)
    r2 = r2+alpha2
    !--- Generate localized posterior particles and calculate mean ---
    pfm = 0.0d0
    DO i = 1,PARA_DA_MEM
      xa(:,i) = xmpf + r1*(xs(:,i)-xmpf) + r2*(xf(:,i)-xmpf)
      pfm = pfm + xa(:,i)
    ENDDO
    pfm = pfm/dble(PARA_DA_MEM)
    !--- Recenter on weight-estimated mean, since sampling error in merging ---
    !--- step can cause deviations from this estimate. ---
    DO i = 1,PARA_DA_MEM
      xa(:,i) = xmpf + (xa(:,i) - pfm)
    ENDDO
    !--- Finalize ---
    DO i = 1,PARA_DA_MEM
      DO j = 1,in_num
        IF (isnan(xa(j,i))) xa(j,i) = xf(j,i)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE pf_merge

  SUBROUTINE find_obs_infl(lw,alpha)
    IMPLICIT NONE
    ! in & out
    REAL(r_size),INTENT(IN)  :: lw(PARA_DA_MEM)
    REAL(r_size),INTENT(OUT) :: alpha
    REAL(r_size),PARAMETER   :: tol = 1e-3 ! Apply bisection method to find k
    REAL(r_size) :: minwt
    REAL(r_size) :: minw
    REAL(r_size) :: beta_max
    REAL(r_size) :: w(PARA_DA_MEM)
    REAL(r_size) :: Neff_init, Nf_new
    REAL(r_size) :: var_k (3) ! start, middle, end
    REAL(r_size) :: var_fk(3) ! start, middle, end

    INTEGER :: i,j,k,l,m,n
    INTEGER :: counter1,counter2,counter3
    INTEGER                                     :: tmp_num1,tmp_num2,tmp_num3,tmp_num4
    REAL(r_size)                                :: tmp_var1,tmp_var2,tmp_var3,tmp_var4
    REAL(r_size),ALLOCATABLE,DIMENSION(:)       :: tmp_d1_mtrx1,tmp_d1_mtrx2,tmp_d1_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: tmp_d2_mtrx1,tmp_d2_mtrx2,tmp_d2_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: tmp_d3_mtrx1,tmp_d3_mtrx2,tmp_d3_mtrx3
    REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:,:) :: tmp_d4_mtrx1,tmp_d4_mtrx2,tmp_d4_mtrx3

    minwt = (1e-10)/PARA_DA_MEM
    w  = exp(-lw)
    IF (sum(w) > 0) THEN
      w    = w/sum(w)
      minw = minval(w)
    ELSE
      minw = 0.0d0
    ENDIF

    IF (minw > minwt) THEN
      alpha = 1.0d0
      RETURN
    ENDIF

    ! Initial start and end bounds
    var_k(1) = 1            ! start
    var_k(3) = maxval(lw)   ! end

    DO i = 1,1000
      var_k(2) = (var_k(1)+var_k(3))/2 ! middle

      ! Evaluate function at each point
      ! start
      w = exp(-lw/var_k(1))
      w = w/sum(w)
      var_fk(1) = minval(w) - minwt
      ! middle 
      w = exp(-lw/var_k(2))
      w = w/sum(w)
      var_fk(2) = minval(w) - minwt
      ! end
      w = exp(-lw/var_k(3))
      w = w/sum(w)
      var_fk(3) = minval(w) - minwt

      ! Exit criteria
      IF ((var_k(3)-var_k(1))/2 < tol) THEN
        EXIT
      ENDIF

      IF (var_fk(1)*var_fk(2) > 0) THEN
        var_k(1) = var_k(2)
      ELSE
        var_k(3) = var_k(2)
      ENDIF
    ENDDO

    alpha = var_k(2)

    RETURN
  END SUBROUTINE find_obs_infl

!-----------------------------------------------------------------------
! The Shapiro-Wilk test (Shapiro and Wilk, 1965)
! Kenta Kurosawa, June 21, 2021
!-----------------------------------------------------------------------
SUBROUTINE CAL_MIN_RES_SW_TEST(dim1,SAMPLE_SIZE,xf,loc,min_res)
  IMPLICIT NONE
  !-- in & out
  INTEGER     ,INTENT(IN)  :: dim1
  INTEGER     ,INTENT(IN)  :: SAMPLE_SIZE !  ens size
  REAL(r_size),INTENT(IN)  :: xf     (dim1,SAMPLE_SIZE)
  REAL(r_size),INTENT(IN)  :: loc    (dim1,dim1)
  REAL(r_size),INTENT(OUT) :: min_res(dim1)
  !-- joint pdf
  INTEGER     ,ALLOCATABLE,DIMENSION(:)   :: JPDF_cnt
  INTEGER     ,ALLOCATABLE,DIMENSION(:,:) :: JPDF_ind
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: JPDF_loc, JPDF_dat
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: tmp_dat, tmp_loc
  !-- 
  INTEGER :: i,j,k,l,m,n
  INTEGER :: counter1,counter2,counter3
  INTEGER                                     :: tmp_num1,tmp_num2,tmp_num3,tmp_num4
  REAL(r_size)                                :: tmp_var1,tmp_var2,tmp_var3,tmp_var4,tmp_var5
  REAL(r_size),ALLOCATABLE,DIMENSION(:)       :: work13,work14,work15,work16
  REAL(r_size),ALLOCATABLE,DIMENSION(:)       :: tmp_d1_mtrx1,tmp_d1_mtrx2,tmp_d1_mtrx3
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: tmp_d2_mtrx1,tmp_d2_mtrx2,tmp_d2_mtrx3
  INTEGER :: my_num


  !--- ALLOCATE
  !- 1. The Shapiro-Wilk test (Shapiro and Wilk, 1965) (bivariational case) -
  !$omp parallel do private(i,tmp_d1_mtrx1,work13,work14,tmp_var5,tmp_d2_mtrx1,j)
  DO i = 1,dim1
    !--- loc
    ALLOCATE(tmp_d1_mtrx1(dim1))          ! tmp_loc_fact
    ALLOCATE(work13(dim1)); work13 = 0.0d0 ! loc_fact
    ALLOCATE(work14(dim1)); work14 = UNDEF ! 0 or 1
    tmp_d1_mtrx1 = loc(:,i)
!    counter1 = 0
    DO j = 1,dim1
      IF (tmp_d1_mtrx1(j) > 0.001d0) THEN
!        counter1 = counter1 + 1
        ALLOCATE(tmp_d2_mtrx1(2,SAMPLE_SIZE)) ! tmp_xf
        work13(j) = tmp_d1_mtrx1(j)
        tmp_d2_mtrx1(1,:) = xf(i,:)
        tmp_d2_mtrx1(2,:) = xf(j,:)
        CALL SW_ROYSTON_TEST(tmp_d2_mtrx1,SAMPLE_SIZE,2,work14(j)) ! Royston's Multivariate Normality Test
        DEALLOCATE(tmp_d2_mtrx1)
      ENDIF
    ENDDO
!    write(*,*) counter1
!    print *
!    read *
    tmp_var5 = dot_product(work13,work14)/sum(work13)
!    tmp_var5 = max(tmp_var5,0.05d0)
!    tmp_var5 = min(tmp_var5,0.95d0)
    min_res(i) = tmp_var5
!!!    my_num = omp_get_thread_num()
    IF (myrank==0) THEN
      WRITE(*,*) i,'/',dim1, min_res(i)
    ENDIF

    DEALLOCATE(tmp_d1_mtrx1)
    DEALLOCATE(work13,work14)
  ENDDO
  !$omp end parallel do

END SUBROUTINE CAL_MIN_RES_SW_TEST


SUBROUTINE CONV_MINRES(dim1,in_minres,out_minres)
  IMPLICIT NONE
  !-- in & out
  INTEGER,INTENT(IN)  :: dim1
  REAL(r_size),INTENT(IN)  :: in_minres(dim1)
  REAL(r_size),INTENT(OUT) :: out_minres(dim1)
  REAL(r_size) :: conv_para1,conv_para2,conv_para3,conv_para4
  INTEGER  :: i
  conv_para1 = 0.3d0 ! x
  conv_para2 = 0.7d0 ! x
  conv_para3 = 0.3d0 ! y
  conv_para4 = 0.9d0 ! y
  DO i = 1,dim1
    IF (in_minres(i) <= conv_para1) THEN
      out_minres(i) = conv_para3
    ELSEIF (in_minres(i) >= conv_para2) THEN
      out_minres(i) = conv_para4
    ELSE
      out_minres(i) = (conv_para4-conv_para3)/(conv_para2-conv_para1)*(in_minres(i)-conv_para1)+conv_para3
    ENDIF
  ENDDO
END SUBROUTINE CONV_MINRES
!-----------------------------------------------------------------------
! Royston's Multivariate Normality Test.
! Kenta Kurosawa, June 21, 2021
! ORIGINAL: Roystest.m
!           https://jp.mathworks.com/matlabcentral/fileexchange/17811-roystest
!-----------------------------------------------------------------------
! ... The Shapiro-Wilk test (Shapiro and Wilk, 1965), is generally considered
! to be an excellent test of univariate normality. It is only natural to
! extend it to the multivariate case, as done by Royston (1982). ... 
!-----------------------------------------------------------------------
!!!!!!SUBROUTINE SW_ROYSTON_TEST(IN_X,OUT_RESULT)
!!!!!!  IMPLICIT NONE
!!!!!!  !-- in & out
!!!!!!  REAL(r_size),INTENT(IN)  :: IN_X(2,PARA_DA_MEM)
!!!!!!  REAL(r_size),INTENT(OUT)  :: OUT_RESULT      ! 1 or 0,  1:not rejected (gaussian), 0:rejected (non-gaussian)
!!!!!!  !-- Shapiro-Wilk test
!!!!!!  REAL(r_size),PARAMETER :: alpha = 0.05d0     ! significance level (default = 0.05)
!!!!!!  REAL(r_size) :: tmp_W(2), tmp_Z(2), tmp_R(2) ! bivriational
!!!!!!  REAL(r_size) :: x, g, m, s
!!!!!!  REAL(r_size) :: u, v, l 
!!!!!!  REAL(r_size) :: T, mC, e, H, P
!!!!!!  REAL(r_size),DIMENSION(2,2) :: C  ! correlation matrix
!!!!!!  REAL(r_size),DIMENSION(2,2) :: NC ! transformed correlation matrix
!!!!!!  !-- cdf 
!!!!!!  REAL(r_size) :: r8_normal_01_cdf_inverse, alnorm
!!!!!!  LOGICAL      :: upper
!!!!!!  REAL(r_size) :: tmp_cdf,bound
!!!!!!  INTEGER(kind=4) :: status
!!!!!!  !-- 
!!!!!!  INTEGER :: i,j,k
!!!!!!  REAL(r_size)                        :: work01, work02
!!!!!!  REAL(r_size),DIMENSION(PARA_DA_MEM) :: work11
!!!!!!  
!!!!!!  IF (PARA_DA_MEM<3) THEN
!!!!!!    WRITE(*,*) 'sample size is too small in SW_ROYSTON_TEST.'
!!!!!!    PRINT *
!!!!!!    READ *
!!!!!!  ELSEIF (PARA_DA_MEM >= 4 .AND. PARA_DA_MEM <= 11) THEN
!!!!!!    x = dble(PARA_DA_MEM)
!!!!!!    g = -2.273d0 + 0.459d0*x
!!!!!!    m = 0.5440d0 - 0.39978d0*x + 0.025054d0*x**2 - 0.0006714d0*x**3
!!!!!!    s = exp(1.3822d0 - 0.77857d0*x + 0.062767d0*x**2 - 0.0020322d0*x**3)
!!!!!!    DO j = 1,2 ! bivariational
!!!!!!      CALL ShaWilstat(IN_X(j,:),tmp_W(j))
!!!!!!      tmp_Z(j) = (-log(g-(log(1-tmp_W(j))))-m)/s
!!!!!!    ENDDO
!!!!!!  ELSEIF (PARA_DA_MEM >= 12 .AND. PARA_DA_MEM <= 2000) THEN
!!!!!!    x = log(dble(PARA_DA_MEM))
!!!!!!    g = 0.0d0
!!!!!!    m = -1.5861d0 - 0.31082d0*x - 0.083751d0*x**2 + 0.0038915d0*x**3
!!!!!!    s = exp(-0.4803d0 -0.082676d0*x + 0.0030302d0*x**2)
!!!!!!    DO j = 1,2 ! bivariational
!!!!!!      CALL ShaWilstat(IN_X(j,:),tmp_W(j))
!!!!!!      tmp_Z(j) = ((log(1-tmp_W(j)))+g-m)/s
!!!!!!      IF (tmp_W(j)>0.95d0 .AND. tmp_Z(j)>2.0d0) tmp_Z(j) = 0.0d0 ! KKUROSAWA
!!!!!!    ENDDO
!!!!!!  ELSE 
!!!!!!    WRITE(*,*) 'sample size is too large in SW_ROYSTON_TEST.'
!!!!!!    PRINT *
!!!!!!    READ *
!!!!!!  ENDIF
!!!!!!
!!!!!!  DO j = 1,2 ! bivariational
!!!!!!  !  work01 = alnorm(-tmp_Z(j),upper)
!!!!!!    CALL cal_alnorm(-tmp_Z(j),.TRUE.,work01)
!!!!!!    work01 = work01/2.0d0
!!!!!!    work02 = r8_normal_01_cdf_inverse(work01)
!!!!!!    work02 = work02**2
!!!!!!    tmp_R(j) = work02
!!!!!!  ENDDO
!!!!!!
!!!!!!  u = 0.715d0
!!!!!!  v = 0.21364d0 + 0.015124d0*(log(dble(PARA_DA_MEM)))**2 - 0.0018034d0*(log(dble(PARA_DA_MEM)))**3
!!!!!!  l = 5.0d0
!!!!!!  IF ((abs(minval(IN_X(1,:))-minval(IN_X(2,:)))<0.0001d0)  .AND. &
!!!!!!      (abs(maxval(IN_X(1,:))-maxval(IN_X(2,:)))<0.0001d0)) THEN
!!!!!!    C(1,2) = 1.00d0
!!!!!!  ELSE
!!!!!!    CALL com_correl(PARA_DA_MEM,IN_X(1,:),IN_X(2,:),C(1,2))
!!!!!!    C(1,2) = min(C(1,2),1.00d0)
!!!!!!  ENDIF
!!!!!!  C(1,1) = 1.0d0
!!!!!!  C(2,1) = C(1,2)
!!!!!!  C(2,2) = 1.0d0
!!!!!!
!!!!!!  CALL com_mean(2,tmp_W,work02)
!!!!!!  IF (C(1,2)>0.90d0 .AND. work02>0.90) THEN ! KKUROSAWA
!!!!!!    P = alpha
!!!!!!  ELSE
!!!!!!    work01 = 0.0d0
!!!!!!    DO i = 1,2
!!!!!!      DO j = 1,2
!!!!!!        NC(i,j) = (C(i,j)**l)*(1.0d0-(u*(1.0d0-C(i,j))**u)/v)
!!!!!!        work01  = work01 + NC(i,j)
!!!!!!      ENDDO
!!!!!!    ENDDO
!!!!!!    T  = work01-dble(2)         ! total
!!!!!!    mC = T/(dble(2)**2-dble(2)) ! average correlation
!!!!!!    e  = dble(2)/(1.0d0+(dble(2)-1.0d0)*mC) ! equivalent degrees of freedom
!!!!!!    H  = e*sum(tmp_R)/dble(2)
!!!!!!    IF (isnan(H) .OR. isnan(e)) THEN
!!!!!!      WRITE(*,*) 'ERROR in SW_ROYSTON_TEST!!'
!!!!!!      write(*,*) 'tmp_W', tmp_W
!!!!!!      write(*,*) 'tmp_Z', tmp_Z
!!!!!!      write(*,*) 'tmp_R', tmp_R
!!!!!!      write(*,*) 'IN_X(1,:)', IN_X(1,:)
!!!!!!      write(*,*) '---------------------------'
!!!!!!      write(*,*) 'IN_X(2,:)', IN_X(2,:)
!!!!!!      write(*,*) 'work01', work01
!!!!!!      write(*,*) 'NC', NC
!!!!!!      write(*,*) 'C', C
!!!!!!      write(*,*) 'T', T
!!!!!!      write(*,*) 'mC', mC
!!!!!!      PRINT *
!!!!!!      READ *
!!!!!!    ENDIF
!!!!!!    CALL cdfchi(1,tmp_cdf,work01,H,e,status,bound)
!!!!!!    P = work01
!!!!!!  ENDIF
!!!!!!  !--- Finalize ---
!!!!!!  IF (P>=alpha) THEN ! Data analyzed have a normal distribution.
!!!!!!    OUT_RESULT = 1.0d0
!!!!!!  ELSE ! Data analyzed do not have a normal distribution.
!!!!!!    OUT_RESULT = 0.0d0
!!!!!!  ENDIF
!!!!!!  RETURN 
!!!!!!END SUBROUTINE SW_ROYSTON_TEST
!!!!!!
!!!!!!!-----------------------------------------------------------------------
!!!!!!! Shapiro-Wilk' W statistic for assessing a sample normality
!!!!!!! Kenta Kurosawa, June 21, 2021
!!!!!!! ORIGINAL: ShaWilstat in Roystest.m
!!!!!!!           https://jp.mathworks.com/matlabcentral/fileexchange/17811-roystest
!!!!!!!-----------------------------------------------------------------------
!!!!!!SUBROUTINE ShaWilstat(IN_X,OUT_W)
!!!!!!  IMPLICIT NONE
!!!!!!  !-- in & out
!!!!!!  REAL(r_size),INTENT(IN)  :: IN_X(PARA_DA_MEM)
!!!!!!  REAL(r_size),INTENT(OUT) :: OUT_W   
!!!!!!  !-- Shapiro-Wilk test
!!!!!!  REAL(r_size) :: tmp_x(PARA_DA_MEM)
!!!!!!  REAL(r_size) :: tmp_m(PARA_DA_MEM)
!!!!!!  REAL(r_size) :: tmp_w(PARA_DA_MEM)
!!!!!!  REAL(r_size) :: tmp_c(PARA_DA_MEM)
!!!!!!  REAL(r_size) :: tmp_u, phi
!!!!!!  INTEGER      :: ct
!!!!!!  REAL(r_size) :: kurt
!!!!!!  REAL(r_size) :: W
!!!!!!  REAL(r_size),DIMENSION(6) :: p1, p2
!!!!!!  !--
!!!!!!  INTEGER :: i,j,k
!!!!!!  REAL(r_size)                        :: work01, work02
!!!!!!  REAL(r_size),DIMENSION(PARA_DA_MEM) :: work11
!!!!!!  !--
!!!!!!  real ( kind = 8 ) r8_normal_01_cdf_inverse
!!!!!!
!!!!!!  tmp_x = IN_X
!!!!!!  CALL qsort(PARA_DA_MEM,tmp_x)
!!!!!!  DO i = 1,PARA_DA_MEM
!!!!!!    work01 = (dble(i)-3.0d0/8.0d0)/(dble(PARA_DA_MEM)+0.25d0)
!!!!!!    tmp_m(i) = r8_normal_01_cdf_inverse(work01)
!!!!!!  ENDDO
!!!!!!  tmp_w = 0.0d0
!!!!!!  CALL CAL_KURTOSIS(PARA_DA_MEM,tmp_x,kurt)
!!!!!!  IF (kurt > 3.0d0) THEN ! %Shapiro-Francia test is better for leptokurtic samples
!!!!!!    work01   = dot_product(tmp_m,tmp_m)
!!!!!!    tmp_w    = 1/sqrt(work01)*tmp_m
!!!!!!    CALL com_mean(PARA_DA_MEM,tmp_x,work01)
!!!!!!    work11 = tmp_x - work01
!!!!!!    work01 = dot_product(work11,work11) 
!!!!!!    work02 = dot_product(tmp_w ,tmp_x) 
!!!!!!    work02 = work02**2
!!!!!!    W = work02/work01
!!!!!!  ELSE ! Shapiro-Wilk test is better for platykurtic samples
!!!!!!    work01 = dot_product(tmp_m,tmp_m)
!!!!!!    tmp_c  = 1/sqrt(work01)*tmp_m
!!!!!!    tmp_u  = 1/sqrt(dble(PARA_DA_MEM))
!!!!!!    p1 = (/-2.706056d0,4.434685d0,-2.071190d0,-0.147981d0,0.221157d0,tmp_c(PARA_DA_MEM)/)
!!!!!!    p2 = (/-3.582633d0,5.682633d0,-1.752461d0,-0.293762d0,0.042981d0,tmp_c(PARA_DA_MEM-1)/)
!!!!!!    tmp_w(PARA_DA_MEM) = p1(1)*tmp_u**5+p1(2)*tmp_u**4+p1(3)*tmp_u**3+p1(4)*tmp_u**2+p1(5)*tmp_u**1+p1(6)
!!!!!!    tmp_w(1) = -tmp_w(PARA_DA_MEM)
!!!!!!    IF (PARA_DA_MEM==3) THEN
!!!!!!      tmp_w(1) = 0.707106781d0
!!!!!!      tmp_w(PARA_DA_MEM) = -tmp_w(1)
!!!!!!    ENDIF
!!!!!!    IF (PARA_DA_MEM>=6) THEN
!!!!!!      ct = 3
!!!!!!      tmp_w(PARA_DA_MEM-1) = p2(1)*tmp_u**5+p2(2)*tmp_u**4+p2(3)*tmp_u**3+p2(4)*tmp_u**2+p2(5)*tmp_u**1+p2(6)
!!!!!!      tmp_w(2) = -tmp_w(PARA_DA_MEM-1)
!!!!!!      work01 = dot_product(tmp_m,tmp_m) - 2*tmp_m(PARA_DA_MEM)**2 - 2*tmp_m(PARA_DA_MEM-1)**2
!!!!!!      work02 = 1.0d0 - 2*tmp_w(PARA_DA_MEM)**2 - 2*tmp_w(PARA_DA_MEM-1)**2
!!!!!!      phi = work01/work02
!!!!!!    ELSE
!!!!!!      ct = 2
!!!!!!      work01 = dot_product(tmp_m,tmp_m) - 2*tmp_m(PARA_DA_MEM)**2
!!!!!!      work02 = 1.0d0 - 2*tmp_w(PARA_DA_MEM)**2
!!!!!!      phi = work01/work02
!!!!!!    ENDIF
!!!!!!    tmp_w(ct:PARA_DA_MEM-ct+1) = tmp_m(ct:PARA_DA_MEM-ct+1)/sqrt(phi)
!!!!!!    CALL com_mean(PARA_DA_MEM,tmp_x,work01)
!!!!!!    work11 = tmp_x - work01
!!!!!!    work01 = dot_product(work11,work11) 
!!!!!!    work02 = dot_product(tmp_w ,tmp_x) 
!!!!!!    work02 = work02**2
!!!!!!    W = work02/work01
!!!!!!  !    write(*,*) 'phi',IN_x
!!!!!!  ENDIF
!!!!!!
!!!!!!  OUT_W = W
!!!!!!  RETURN
!!!!!!
!!!!!!END SUBROUTINE ShaWilstat
!!!!!!
!!!!!!!-----------------------------------------------------------------------
!!!!!!! KURTOSIS
!!!!!!! Kenta Kurosawa, June 21, 2021
!!!!!!! ORIGINAL: kurtosis.m
!!!!!!!-----------------------------------------------------------------------
!!!!!!SUBROUTINE CAL_KURTOSIS(dim1,IN_X,OUT_KURT)
!!!!!!  IMPLICIT NONE
!!!!!!  !-- in & out
!!!!!!  INTEGER,INTENT(IN)       :: dim1
!!!!!!  REAL(r_size),INTENT(IN)  :: IN_X(dim1)
!!!!!!  REAL(r_size),INTENT(OUT) :: OUT_KURT
!!!!!!  !--    
!!!!!!  INTEGER :: j, n
!!!!!!  REAL(r_size) :: tmp_var1, tmp_var2(dim1)
!!!!!!  REAL(r_size) :: x0(dim1), s2, m4, k
!!!!!!
!!!!!!  CALL com_mean(dim1,IN_X,tmp_var1)
!!!!!!  x0 = IN_X - tmp_var1 
!!!!!!  CALL com_mean(dim1,x0**2,s2)
!!!!!!  CALL com_mean(dim1,x0**4,m4)
!!!!!!  k = m4/s2**2  
!!!!!!
!!!!!!  OUT_KURT = k
!!!!!!  RETURN
!!!!!!END SUBROUTINE CAL_KURTOSIS

END MODULE LPF_tools

