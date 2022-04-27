MODULE Serial_EnSRF_tools
  USE common
  USE common_mtx
  USE common_L96DA
  USE common_L96
  use f95_lapack

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Serial_EnSRF_CORE_tempered, Serial_EnSRF_CORE_tempered_v2, Serial_EnSRF_CORE_tempered_v3

  CONTAINS

!!!!!-----------------------------------------
!!!!  SUBROUTINE Serial_EnSRF_CORE(xf,loc_fact,obsnum,H,Yo,R,xa,xf_KS,MDA_alpha)
!!!!    IMPLICIT NONE
!!!!    REAL(r_size),INTENT(IN)  :: xf       (PARA_L96_J,PARA_DA_MEM)
!!!!    REAL(r_size),INTENT(IN)  :: loc_fact (PARA_L96_J,PARA_L96_J)
!!!!    INTEGER     ,INTENT(IN)  :: obsnum
!!!!    INTEGER     ,INTENT(IN)  :: H        (PARA_obs_num,PARA_L96_J)
!!!!    REAL(r_size),INTENT(IN)  :: Yo       (PARA_obs_num)
!!!!    REAL(r_size),INTENT(IN)  :: R        (PARA_obs_num,PARA_obs_num)
!!!!    REAL(r_size),INTENT(OUT) :: xa       (PARA_L96_J,PARA_DA_MEM)
!!!!    REAL(r_size),INTENT(IN),OPTIONAL  :: xf_KS  (PARA_L96_J,PARA_DA_MEM)
!!!!    REAL(r_size),INTENT(IN),OPTIONAL  :: MDA_alpha
!!!!    REAL(r_size) :: xm    (PARA_L96_J)
!!!!    REAL(r_size) :: xp    (PARA_L96_J,PARA_DA_MEM)
!!!!    REAL(r_size) :: hx    (PARA_obs_num,PARA_DA_MEM)
!!!!    REAL(r_size) :: hxm   (PARA_obs_num)
!!!!    REAL(r_size) :: hxp   (PARA_obs_num,PARA_DA_MEM)
!!!!    REAL(r_size) :: tmp_R (PARA_obs_num)
!!!!    REAL(r_size) :: inov
!!!!    REAL(r_size) :: hxo   (PARA_DA_MEM)
!!!!    REAL(r_size) :: tmp_matrix1
!!!!    REAL(r_size),ALLOCATABLE :: P (:)
!!!!    REAL(r_size),ALLOCATABLE :: KalmanGain    (:)
!!!!    REAL(r_size),ALLOCATABLE :: KalmanGain_new(:)
!!!!    REAL(r_size) :: tmp_fact
!!!!    REAL(r_size) :: alpha
!!!!    INTEGER :: k_posi,k_posi2
!!!!    INTEGER :: i,j,k
!!!!    INTEGER :: QC(PARA_obs_num) ! <-- !!
!!!!  
!!!!    IF     (PARA_obs_nonlinear==1) THEN
!!!!      hx = matmul(H,xf**2)
!!!!    ELSEIF (PARA_obs_nonlinear==2) THEN
!!!!      hx = matmul(H,log(abs(xf)))
!!!!    ELSE
!!!!      hx = matmul(H,xf)
!!!!    ENDIF
!!!!  
!!!!    !--- model ---
!!!!    DO i = 1,PARA_L96_J  
!!!!      CALL com_mean(PARA_DA_MEM,xf(i,:),xm(i))
!!!!    ENDDO
!!!!    DO i = 1,PARA_DA_MEM
!!!!      xp(:,i) = xf(:,i) - xm 
!!!!    ENDDO
!!!!  
!!!!    !--- obs ---
!!!!    DO i = 1,PARA_obs_num 
!!!!      CALL com_mean(PARA_DA_MEM,hx(i,:),hxm(i))
!!!!      tmp_R(i) = R(i,i)
!!!!    ENDDO
!!!!    DO i = 1,PARA_DA_MEM
!!!!      hxp(:,i) = hx(:,i) - hxm
!!!!    ENDDO
!!!!  
!!!!    !--- ESMDA ---
!!!!    IF (present(MDA_alpha)) THEN
!!!!    IF (PARA_DA_MODE == 'ESMDA' .OR. PARA_DA_DAwindow_obstimes>0 .OR. PARA_DA_LAG>0) THEN
!!!!      DO i = 1,PARA_L96_J  
!!!!        CALL com_mean(PARA_DA_MEM,xf_KS(i,:),xm(i))
!!!!      ENDDO
!!!!      DO i = 1,PARA_DA_MEM
!!!!        xp(:,i) = xf_KS(:,i) - xm 
!!!!      ENDDO
!!!!      tmp_R = tmp_R*MDA_alpha
!!!!    ENDIF
!!!!    ENDIF
!!!!
!!!!    !--- QC ---
!!!!    QC = 0
!!!!    DO i = 1,PARA_obs_num
!!!!      inov = Yo(i) - hxm(i)
!!!!      CALL com_stdev(PARA_DA_MEM,hx(i,:),tmp_matrix1)
!!!!      tmp_matrix1 = tmp_matrix1**2 ! variance
!!!!      IF (inov > 4.0d0*sqrt(tmp_matrix1+tmp_R(i))) THEN
!!!!        QC(i) = 1
!!!!      ENDIF
!!!!    ENDDO
!!!!      
!!!!    !--- MAIN ---
!!!!    DO i = 1,PARA_obs_num
!!!!
!!!!      IF (QC(i)==0) THEN
!!!!     
!!!!        ! Position 
!!!!        DO j = 1,PARA_L96_J
!!!!          IF (H(i,j)==1) THEN
!!!!            k_posi = j
!!!!            EXIT
!!!!          ENDIF
!!!!        ENDDO
!!!!    
!!!!        ! Innovations
!!!!        inov = Yo(i) - hxm(i)
!!!!        hxo  = hxp(i,:)
!!!!        tmp_matrix1 = 1/(dot_product(hxo,hxo)/dble(PARA_DA_MEM-1)+tmp_R(i)) ! var_den
!!!!    
!!!!        !--- 1. Model Space Update ---
!!!!        ! Kalman Gain
!!!!        ALLOCATE(P             (PARA_L96_J))
!!!!        ALLOCATE(KalmanGain    (PARA_L96_J))
!!!!        ALLOCATE(KalmanGain_new(PARA_L96_J))
!!!!        P = matmul(xp,hxo)/dble(PARA_DA_MEM-1)
!!!!        KalmanGain = P*tmp_matrix1
!!!!        ! Localization 
!!!!        DO j = 1,PARA_L96_J
!!!!          tmp_fact = loc_fact(k_posi,j)
!!!!          KalmanGain(j) = KalmanGain(j)*tmp_fact
!!!!        ENDDO
!!!!        ! xa mean
!!!!        xm = xm + KalmanGain*inov
!!!!        ! Perturbation Update
!!!!        alpha = 1/(1+sqrt(tmp_R(i)*tmp_matrix1))
!!!!        KalmanGain_new = KalmanGain*alpha
!!!!        DO j = 1,PARA_L96_J
!!!!          DO k = 1,PARA_DA_MEM
!!!!            xp(j,k) = xp(j,k) - KalmanGain_new(j)*hxo(k)
!!!!          ENDDO
!!!!        ENDDO
!!!!        DEALLOCATE(P,KalmanGain,KalmanGain_new)
!!!!  
!!!!        !--- 2. Obs Space Update ---
!!!!        ! Kalman Gain
!!!!        ALLOCATE(P             (PARA_obs_num))
!!!!        ALLOCATE(KalmanGain    (PARA_obs_num))
!!!!        ALLOCATE(KalmanGain_new(PARA_obs_num))
!!!!        P = matmul(hxp,hxo)/dble(PARA_DA_MEM-1)
!!!!        KalmanGain = P*tmp_matrix1
!!!!        ! Localization 
!!!!        DO j = 1,PARA_obs_num
!!!!          DO k = 1,PARA_L96_J
!!!!            IF (H(j,k)==1) THEN
!!!!              k_posi2 = k
!!!!              tmp_fact = loc_fact(k_posi,k_posi2)
!!!!              KalmanGain(j) = KalmanGain(j)*tmp_fact
!!!!            ENDIF
!!!!          ENDDO
!!!!        ENDDO
!!!!        ! xa mean
!!!!        hxm = hxm + KalmanGain*inov
!!!!        ! Perturbation Update
!!!!        KalmanGain_new = KalmanGain*alpha
!!!!        DO j = 1,PARA_obs_num
!!!!          DO k = 1,PARA_DA_MEM
!!!!            hxp(j,k) = hxp(j,k) - KalmanGain_new(j)*hxo(k)
!!!!          ENDDO
!!!!        ENDDO
!!!!        DEALLOCATE(P,KalmanGain,KalmanGain_new)
!!!!      
!!!!      ENDIF
!!!!    ENDDO
!!!!
!!!!    !--- Finalize ---
!!!!    DO i = 1,PARA_DA_MEM
!!!!      xa(:,i) = xm + xp(:,i)
!!!!    ENDDO
!!!!  
!!!!  return
!!!!  END SUBROUTINE Serial_EnSRF_CORE
!!!!
!-----------------------------------------
! only for LPF by Jon Poterjoy !!!
!-----------------------------------------
! Function for performing enkf update using Whitaker and Hamill filter with                             
! Anderson adaptive state-space inflation
!
!  INPUT 
!           x: prior ensemble (Nx x Ne)
!          hx: obs-space prior ensemble (Ny x Ne)
!           y: observation vector  (Ny x 1)
!       var_y: obs error variance
!          HC: matrix determining localization between obs- and model-space  (Ny x Nx)
!         HCH: matrix determining localization between obs- and obs-space  (Ny x Ny)
!    inf_flag: flag for inflation option
!         inf: prior mean of inf coefficients (Nx x 1)
!     var_inf: prior variance of inf coefficients
!       gamma: relaxation parameter for RTPP 
!
!  OUTPUT 
!          xm: posterior mean (Nx x 1)
!           x: posterior ensemble (Nx x Ne)
!         inf: posterior mean of inf coefficients (Nx x 1)
!      e_flag: error flag

!!!!  SUBROUTINE Serial_EnSRF_CORE_tempered(xf,loc_fact,obsnum,H,Yo,R,obs_infl,xa)
!!!!    IMPLICIT NONE
!!!!    REAL(r_size),INTENT(IN)  :: xf       (PARA_L96_J,PARA_DA_MEM)
!!!!    REAL(r_size),INTENT(IN)  :: loc_fact (PARA_L96_J,PARA_L96_J)
!!!!    INTEGER     ,INTENT(IN)  :: obsnum
!!!!    INTEGER     ,INTENT(IN)  :: H        (obsnum,PARA_L96_J)
!!!!    REAL(r_size),INTENT(IN)  :: Yo       (obsnum)
!!!!    REAL(r_size),INTENT(IN)  :: R        (obsnum,obsnum)
!!!!    REAL(r_size),INTENT(IN)  :: obs_infl (obsnum)
!!!!    REAL(r_size),INTENT(OUT) :: xa       (PARA_L96_J,PARA_DA_MEM)
!!!!    REAL(r_size) :: xm    (PARA_L96_J)
!!!!    REAL(r_size) :: xp    (PARA_L96_J,PARA_DA_MEM)
!!!!    REAL(r_size) :: hx    (obsnum,PARA_DA_MEM)
!!!!    REAL(r_size) :: hxm   (obsnum)
!!!!    REAL(r_size) :: hxp   (obsnum,PARA_DA_MEM)
!!!!    REAL(r_size) :: tmp_R (obsnum)
!!!!    REAL(r_size) :: inov
!!!!    REAL(r_size) :: hxo   (PARA_DA_MEM)
!!!!    REAL(r_size) :: tmp_matrix1
!!!!    REAL(r_size),ALLOCATABLE :: P (:)
!!!!    REAL(r_size),ALLOCATABLE :: KalmanGain    (:)
!!!!    REAL(r_size),ALLOCATABLE :: KalmanGain_new(:)
!!!!    REAL(r_size) :: tmp_fact
!!!!    REAL(r_size) :: alpha
!!!!    INTEGER :: k_posi,k_posi2
!!!!    INTEGER :: i,j,k
!!!!  
!!!!    IF     (PARA_obs_nonlinear==1) THEN
!!!!      hx = matmul(H,xf**2)
!!!!    ELSEIF (PARA_obs_nonlinear==2) THEN
!!!!      hx = matmul(H,log(abs(xf)))
!!!!    ELSE
!!!!      hx = matmul(H,xf)
!!!!    ENDIF
!!!!  
!!!!    !--- model ---
!!!!    DO i = 1,PARA_L96_J  
!!!!      CALL com_mean(PARA_DA_MEM,xf(i,:),xm(i))
!!!!    ENDDO
!!!!    DO i = 1,PARA_DA_MEM
!!!!      xp(:,i) = xf(:,i) - xm 
!!!!    ENDDO
!!!!  
!!!!    !--- obs ---
!!!!    DO i = 1,obsnum 
!!!!      CALL com_mean(PARA_DA_MEM,hx(i,:),hxm(i))
!!!!      tmp_R(i) = R(i,i)*obs_infl(i) ! <-- LPF !!!
!!!!    ENDDO
!!!!    DO i = 1,PARA_DA_MEM
!!!!      hxp(:,i) = hx(:,i) - hxm
!!!!    ENDDO
!!!!  
!!!!    !--- MAIN ---
!!!!    DO i = 1,obsnum
!!!!     
!!!!      IF (obs_infl(i) > 99999) CYCLE ! <-- LPF !!!
!!!!
!!!!      ! Position 
!!!!      DO j = 1,PARA_L96_J
!!!!        IF (H(i,j)==1) THEN
!!!!          k_posi = j
!!!!          EXIT
!!!!        ENDIF
!!!!      ENDDO
!!!!  
!!!!      ! Innovations
!!!!      inov = Yo(i) - hxm(i)
!!!!      hxo  = hxp(i,:)
!!!!      tmp_matrix1 = 1/(dot_product(hxo,hxo)/dble(PARA_DA_MEM-1)+tmp_R(i)) ! var_den
!!!!  
!!!!      !--- 1. Model Space Update ---
!!!!      ! Kalman Gain
!!!!      ALLOCATE(P             (PARA_L96_J))
!!!!      ALLOCATE(KalmanGain    (PARA_L96_J))
!!!!      ALLOCATE(KalmanGain_new(PARA_L96_J))
!!!!      P = matmul(xp,hxo)/dble(PARA_DA_MEM-1)
!!!!      KalmanGain = P*tmp_matrix1
!!!!      ! Localization 
!!!!      DO j = 1,PARA_L96_J
!!!!        tmp_fact = loc_fact(k_posi,j)
!!!!        KalmanGain(j) = KalmanGain(j)*tmp_fact
!!!!      ENDDO
!!!!      ! xa mean
!!!!      xm = xm + KalmanGain*inov
!!!!      ! Perturbation Update
!!!!      alpha = 1/(1+sqrt(tmp_R(i)*tmp_matrix1))
!!!!      KalmanGain_new = KalmanGain*alpha
!!!!      DO j = 1,PARA_L96_J
!!!!        DO k = 1,PARA_DA_MEM
!!!!          xp(j,k) = xp(j,k) - KalmanGain_new(j)*hxo(k)
!!!!        ENDDO
!!!!      ENDDO
!!!!      DEALLOCATE(P,KalmanGain,KalmanGain_new)
!!!!
!!!!      !--- 2. Obs Space Update ---
!!!!      ! Kalman Gain
!!!!      ALLOCATE(P             (obsnum))
!!!!      ALLOCATE(KalmanGain    (obsnum))
!!!!      ALLOCATE(KalmanGain_new(obsnum))
!!!!      P = matmul(hxp,hxo)/dble(PARA_DA_MEM-1)
!!!!      KalmanGain = P*tmp_matrix1
!!!!      ! Localization 
!!!!      DO j = 1,obsnum
!!!!        DO k = 1,PARA_L96_J
!!!!          IF (H(j,k)==1) THEN
!!!!            k_posi2 = k
!!!!            tmp_fact = loc_fact(k_posi,k_posi2)
!!!!            KalmanGain(j) = KalmanGain(j)*tmp_fact
!!!!          ENDIF
!!!!        ENDDO
!!!!      ENDDO
!!!!      ! xa mean
!!!!      hxm = hxm + KalmanGain*inov
!!!!      ! Perturbation Update
!!!!      KalmanGain_new = KalmanGain*alpha
!!!!      DO j = 1,obsnum
!!!!        DO k = 1,PARA_DA_MEM
!!!!          hxp(j,k) = hxp(j,k) - KalmanGain_new(j)*hxo(k)
!!!!        ENDDO
!!!!      ENDDO
!!!!      DEALLOCATE(P,KalmanGain,KalmanGain_new)
!!!!      
!!!!    ENDDO
!!!!  
!!!!
!!!!    !--- Finalize ---
!!!!    DO i = 1,PARA_DA_MEM
!!!!      xa(:,i) = xm + xp(:,i)
!!!!    ENDDO
!!!!  
!!!!  return
!!!!  END SUBROUTINE Serial_EnSRF_CORE_tempered

  SUBROUTINE Serial_EnSRF_CORE_tempered(xf,loc_fact,obsnum,H,hx,Yo,R,obs_infl,xa)
    IMPLICIT NONE
    REAL(r_size),INTENT(IN)  :: xf       (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)  :: loc_fact (PARA_L96_J,PARA_L96_J)
    INTEGER     ,INTENT(IN)  :: obsnum
    INTEGER     ,INTENT(IN)  :: H        (obsnum,PARA_L96_J)
    REAL(r_size),INTENT(IN)  :: hx(obsnum,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)  :: Yo       (obsnum)
    REAL(r_size),INTENT(IN)  :: R        (obsnum,obsnum)
    REAL(r_size),INTENT(IN)  :: obs_infl (obsnum)
    REAL(r_size),INTENT(OUT) :: xa       (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size) :: xm    (PARA_L96_J)
    REAL(r_size) :: xp    (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size) :: hxm   (obsnum)
    REAL(r_size) :: hxp   (obsnum,PARA_DA_MEM)
    REAL(r_size) :: tmp_R (obsnum)
    REAL(r_size) :: inov
    REAL(r_size) :: hxo   (PARA_DA_MEM)
    REAL(r_size) :: tmp_matrix1
    REAL(r_size),ALLOCATABLE :: P (:)
    REAL(r_size),ALLOCATABLE :: KalmanGain    (:)
    REAL(r_size),ALLOCATABLE :: KalmanGain_new(:)
    REAL(r_size) :: tmp_fact
    REAL(r_size) :: alpha
    INTEGER :: k_posi,k_posi2
    INTEGER :: i,j,k


    !--- model ---
    DO i = 1,PARA_L96_J
      CALL com_mean(PARA_DA_MEM,xf(i,:),xm(i))
    ENDDO
    DO i = 1,PARA_DA_MEM
      xp(:,i) = xf(:,i) - xm
    ENDDO

    !--- obs ---
    DO i = 1,obsnum
      CALL com_mean(PARA_DA_MEM,hx(i,:),hxm(i))
      tmp_R(i) = R(i,i)*obs_infl(i) ! <-- LPF !!!
    ENDDO
    DO i = 1,PARA_DA_MEM
      hxp(:,i) = hx(:,i) - hxm
    ENDDO

    !--- MAIN ---
    DO i = 1,obsnum

      IF (obs_infl(i) > 99999) CYCLE ! <-- LPF !!!

      ! Position
      DO j = 1,PARA_L96_J
        IF (H(i,j)==1) THEN
          k_posi = j
          EXIT
        ENDIF
      ENDDO

      ! Innovations
      inov = Yo(i) - hxm(i)
      hxo  = hxp(i,:)
      tmp_matrix1 = 1/(dot_product(hxo,hxo)/dble(PARA_DA_MEM-1)+tmp_R(i)) ! var_den

      !--- 1. Model Space Update ---
      ! Kalman Gain
      ALLOCATE(P             (PARA_L96_J))
      ALLOCATE(KalmanGain    (PARA_L96_J))
      ALLOCATE(KalmanGain_new(PARA_L96_J))
      P = matmul(xp,hxo)/dble(PARA_DA_MEM-1)
      KalmanGain = P*tmp_matrix1
      ! Localization
      DO j = 1,PARA_L96_J
        tmp_fact = loc_fact(k_posi,j)
        KalmanGain(j) = KalmanGain(j)*tmp_fact
      ENDDO
      ! xa mean
      xm = xm + KalmanGain*inov
      ! Perturbation Update
      alpha = 1/(1+sqrt(tmp_R(i)*tmp_matrix1))
      KalmanGain_new = KalmanGain*alpha
      DO j = 1,PARA_L96_J
        DO k = 1,PARA_DA_MEM
          xp(j,k) = xp(j,k) - KalmanGain_new(j)*hxo(k)
        ENDDO
      ENDDO
      DEALLOCATE(P,KalmanGain,KalmanGain_new)

      !--- 2. Obs Space Update ---
      ! Kalman Gain
      ALLOCATE(P             (obsnum))
      ALLOCATE(KalmanGain    (obsnum))
      ALLOCATE(KalmanGain_new(obsnum))
      P = matmul(hxp,hxo)/dble(PARA_DA_MEM-1)
      KalmanGain = P*tmp_matrix1
      ! Localization
      DO j = 1,obsnum
        DO k = 1,PARA_L96_J
          IF (H(j,k)==1) THEN
            k_posi2 = k
            tmp_fact = loc_fact(k_posi,k_posi2)
            KalmanGain(j) = KalmanGain(j)*tmp_fact
          ENDIF
        ENDDO
      ENDDO
      ! xa mean
      hxm = hxm + KalmanGain*inov
      ! Perturbation Update
      KalmanGain_new = KalmanGain*alpha
      DO j = 1,obsnum
        DO k = 1,PARA_DA_MEM
          hxp(j,k) = hxp(j,k) - KalmanGain_new(j)*hxo(k)
        ENDDO
      ENDDO
      DEALLOCATE(P,KalmanGain,KalmanGain_new)

    ENDDO


    !--- Finalize ---
    DO i = 1,PARA_DA_MEM
      xa(:,i) = xm + xp(:,i)
    ENDDO

  return
  END SUBROUTINE Serial_EnSRF_CORE_tempered

  ! only for adaptive stop
  SUBROUTINE Serial_EnSRF_CORE_tempered_v2(xf,loc_fact,obsnum,H,hx,Yo,R,mdl_infl,obs_infl,xa,out_hxa)
    IMPLICIT NONE
    REAL(r_size),INTENT(IN)  :: xf       (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)  :: loc_fact (PARA_L96_J,PARA_L96_J)
    INTEGER     ,INTENT(IN)  :: obsnum
    INTEGER     ,INTENT(IN)  :: H        (obsnum,PARA_L96_J)
    REAL(r_size),INTENT(IN)  :: hx(obsnum,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)  :: Yo       (obsnum)
    REAL(r_size),INTENT(IN)  :: R        (obsnum,obsnum)
    REAL(r_size),INTENT(IN)  :: mdl_infl (PARA_L96_J)
    REAL(r_size),INTENT(IN)  :: obs_infl (obsnum)
    REAL(r_size),INTENT(OUT) :: xa       (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size),INTENT(OUT),OPTIONAL :: out_hxa (obsnum,PARA_DA_MEM)
    REAL(r_size) :: hxa (obsnum,PARA_DA_MEM)
    REAL(r_size) :: xm    (PARA_L96_J)
    REAL(r_size) :: xp    (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size) :: hxm   (obsnum)
    REAL(r_size) :: hxp   (obsnum,PARA_DA_MEM)
    REAL(r_size) :: tmp_R (obsnum)
    REAL(r_size) :: inov
    REAL(r_size) :: hxo   (PARA_DA_MEM)
    REAL(r_size) :: tmp_matrix1
    REAL(r_size),ALLOCATABLE :: P (:)
    REAL(r_size),ALLOCATABLE :: KalmanGain    (:)
    REAL(r_size),ALLOCATABLE :: KalmanGain_new(:)
    REAL(r_size) :: tmp_fact
    REAL(r_size) :: alpha
    INTEGER :: k_posi,k_posi2
    INTEGER :: i,j,k


    !--- model ---
    DO i = 1,PARA_L96_J
      CALL com_mean(PARA_DA_MEM,xf(i,:),xm(i))
    ENDDO
    DO i = 1,PARA_DA_MEM
      xp(:,i) = xf(:,i) - xm
    ENDDO
    !--- obs ---
    DO i = 1,obsnum
      CALL com_mean(PARA_DA_MEM,hx(i,:),hxm(i))
    !  tmp_R(i) = R(i,i)*obs_infl(i) ! <-- LPF !!!
      tmp_R(i) = R(i,i)
    ENDDO
    DO i = 1,PARA_DA_MEM
      hxp(:,i) = hx(:,i) - hxm
    ENDDO

    !--- MAIN ---
    DO i = 1,obsnum

      IF (obs_infl(i) > 99999) CYCLE ! <-- LPF !!!

      ! Position
      DO j = 1,PARA_L96_J
        IF (H(i,j)==1) THEN
          k_posi = j
          EXIT
        ENDIF
      ENDDO

      ! Innovations
      inov = Yo(i) - hxm(i)
      hxo  = hxp(i,:)
      tmp_matrix1 = 1/(dot_product(hxo,hxo)/dble(PARA_DA_MEM-1)+tmp_R(i)) ! var_den


      !--- 1. Model Space Update ---
      ! Kalman Gain
      ALLOCATE(P             (PARA_L96_J))
      ALLOCATE(KalmanGain    (PARA_L96_J))
      ALLOCATE(KalmanGain_new(PARA_L96_J))
      P = matmul(xp,hxo)/dble(PARA_DA_MEM-1)
      KalmanGain = P*tmp_matrix1
      ! Localization
      DO j = 1,PARA_L96_J
        tmp_fact = loc_fact(k_posi,j)
      !  KalmanGain(j) = KalmanGain(j)*tmp_fact
        KalmanGain(j) = KalmanGain(j)*tmp_fact*mdl_infl(j) ! <- adaptive stop
      ENDDO
      ! xa mean
      xm = xm + KalmanGain*inov
      ! Perturbation Update
      alpha = 1/(1+sqrt(tmp_R(i)*tmp_matrix1))
      KalmanGain_new = KalmanGain*alpha
      DO j = 1,PARA_L96_J
        DO k = 1,PARA_DA_MEM
          xp(j,k) = xp(j,k) - KalmanGain_new(j)*hxo(k)
        ENDDO
      ENDDO
      DEALLOCATE(P,KalmanGain,KalmanGain_new)

      !--- 2. Obs Space Update ---
      ! Kalman Gain
      ALLOCATE(P             (obsnum))
      ALLOCATE(KalmanGain    (obsnum))
      ALLOCATE(KalmanGain_new(obsnum))
      P = matmul(hxp,hxo)/dble(PARA_DA_MEM-1)
      KalmanGain = P*tmp_matrix1
      ! Localization
      DO j = 1,obsnum
        DO k = 1,PARA_L96_J
          IF (H(j,k)==1) THEN
            k_posi2 = k
            tmp_fact = loc_fact(k_posi,k_posi2)
           ! KalmanGain(j) = KalmanGain(j)*tmp_fact
            KalmanGain(j) = KalmanGain(j)*tmp_fact*obs_infl(j) !<- adaptive stop
          ENDIF
        ENDDO
      ENDDO
      ! xa mean
      hxm = hxm + KalmanGain*inov
      ! Perturbation Update
      KalmanGain_new = KalmanGain*alpha
      DO j = 1,obsnum
        DO k = 1,PARA_DA_MEM
          hxp(j,k) = hxp(j,k) - KalmanGain_new(j)*hxo(k)
        ENDDO
      ENDDO
      DEALLOCATE(P,KalmanGain,KalmanGain_new)

    ENDDO


    !--- Finalize ---
    DO i = 1,PARA_DA_MEM
      xa (:,i) = xm  + xp (:,i)
      hxa(:,i) = hxm + hxp(:,i)
    ENDDO
    IF (PRESENT(out_hxa)) THEN
      out_hxa = hxa
    ENDIF


  return
  END SUBROUTINE Serial_EnSRF_CORE_tempered_v2

  ! only for adaptive stop
  SUBROUTINE Serial_EnSRF_CORE_tempered_v3(xf,loc_fact,obsnum,H,hx,Yo,R,mdl_infl,obs_infl,xa,out_hxa)
    IMPLICIT NONE
    REAL(r_size),INTENT(IN)  :: xf       (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)  :: loc_fact (PARA_L96_J,PARA_L96_J)
    INTEGER     ,INTENT(IN)  :: obsnum
    INTEGER     ,INTENT(IN)  :: H        (obsnum,PARA_L96_J)
    REAL(r_size),INTENT(IN)  :: hx(obsnum,PARA_DA_MEM)
    REAL(r_size),INTENT(IN)  :: Yo       (obsnum)
    REAL(r_size),INTENT(IN)  :: R        (obsnum,obsnum)
    REAL(r_size),INTENT(IN)  :: mdl_infl (PARA_L96_J)
    REAL(r_size),INTENT(IN)  :: obs_infl (obsnum)
    REAL(r_size),INTENT(OUT) :: xa       (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size),INTENT(OUT),OPTIONAL :: out_hxa (obsnum,PARA_DA_MEM)
    REAL(r_size) :: hxa (obsnum,PARA_DA_MEM)
    REAL(r_size) :: xm    (PARA_L96_J)
    REAL(r_size) :: xp    (PARA_L96_J,PARA_DA_MEM)
    REAL(r_size) :: hxm   (obsnum)
    REAL(r_size) :: hxp   (obsnum,PARA_DA_MEM)
    REAL(r_size) :: tmp_R (obsnum)
    REAL(r_size) :: inov
    REAL(r_size) :: hxo   (PARA_DA_MEM)
    REAL(r_size) :: tmp_matrix1(PARA_L96_J), tmp_matrix2(obsnum)
    REAL(r_size),ALLOCATABLE :: P (:)
    REAL(r_size),ALLOCATABLE :: KalmanGain    (:)
    REAL(r_size),ALLOCATABLE :: KalmanGain_new(:)
    REAL(r_size) :: tmp_fact
    REAL(r_size) :: alpha_mdl(PARA_L96_J), alpha_obs(obsnum)
    INTEGER :: k_posi,k_posi2
    INTEGER :: i,j,k


    !--- model ---
    DO i = 1,PARA_L96_J
      CALL com_mean(PARA_DA_MEM,xf(i,:),xm(i))
    ENDDO
    DO i = 1,PARA_DA_MEM
      xp(:,i) = xf(:,i) - xm
    ENDDO
    !--- obs ---
    DO i = 1,obsnum
      CALL com_mean(PARA_DA_MEM,hx(i,:),hxm(i))
    !  tmp_R(i) = R(i,i)*obs_infl(i) ! <-- LPF !!!
      tmp_R(i) = R(i,i)
    ENDDO
    DO i = 1,PARA_DA_MEM
      hxp(:,i) = hx(:,i) - hxm
    ENDDO

    !--- MAIN ---
    DO i = 1,obsnum

      IF (obs_infl(i) > 99999) CYCLE ! <-- LPF !!!

      ! Position
      DO j = 1,PARA_L96_J
        IF (H(i,j)==1) THEN
          k_posi = j
          EXIT
        ENDIF
      ENDDO

      ! Innovations
      inov = Yo(i) - hxm(i)
      hxo  = hxp(i,:)
!      tmp_matrix1 = 1/(dot_product(hxo,hxo)/dble(PARA_DA_MEM-1)+tmp_R(i)) ! var_den


      !--- 1. Model Space Update ---
      ! Kalman Gain
      ALLOCATE(P             (PARA_L96_J))
      ALLOCATE(KalmanGain    (PARA_L96_J))
      ALLOCATE(KalmanGain_new(PARA_L96_J))
    !!  P = matmul(xp,hxo)/dble(PARA_DA_MEM-1)
    !!  KalmanGain = P*tmp_matrix1

      DO j = 1,PARA_L96_J
        tmp_matrix1(j) = 1/(dot_product(hxo,hxo)/dble(PARA_DA_MEM-1)+tmp_R(i)*(1.0d0/mdl_infl(j))) ! var_den
        P(j) = dot_product(xp(j,:),hxo)/dble(PARA_DA_MEM-1)
        KalmanGain(j) = P(j)*tmp_matrix1(j)
      ENDDO


      ! Localization
      DO j = 1,PARA_L96_J
        tmp_fact = loc_fact(k_posi,j)
        KalmanGain(j) = KalmanGain(j)*tmp_fact
     !   KalmanGain(j) = KalmanGain(j)*tmp_fact*mdl_infl(j) ! <- adaptive stop
      ENDDO
      ! xa mean
      xm = xm + KalmanGain*inov
      ! Perturbation Update
      DO j = 1,PARA_L96_J
        alpha_mdl(j) = 1/(1+sqrt(tmp_R(i)*tmp_matrix1(j)))
        KalmanGain_new(j) = KalmanGain(j)*alpha_mdl(j)
      ENDDO
      DO j = 1,PARA_L96_J
        DO k = 1,PARA_DA_MEM
          xp(j,k) = xp(j,k) - KalmanGain_new(j)*hxo(k)
        ENDDO
      ENDDO
      DEALLOCATE(P,KalmanGain,KalmanGain_new)

      !--- 2. Obs Space Update ---
      ! Kalman Gain
      ALLOCATE(P             (obsnum))
      ALLOCATE(KalmanGain    (obsnum))
      ALLOCATE(KalmanGain_new(obsnum))
!      P = matmul(hxp,hxo)/dble(PARA_DA_MEM-1)
!      KalmanGain = P*tmp_matrix1

      DO j = 1,obsnum
        tmp_matrix2(j) = 1/(dot_product(hxo,hxo)/dble(PARA_DA_MEM-1)+tmp_R(i)*(1.0d0/mdl_infl(j))) ! var_den
        P(j) = dot_product(xp(j,:),hxo)/dble(PARA_DA_MEM-1)
        KalmanGain(j) = P(j)*tmp_matrix2(j)
      ENDDO

      ! Localization
      DO j = 1,obsnum
        DO k = 1,PARA_L96_J
          IF (H(j,k)==1) THEN
            k_posi2 = k
            tmp_fact = loc_fact(k_posi,k_posi2)
            KalmanGain(j) = KalmanGain(j)*tmp_fact
           ! KalmanGain(j) = KalmanGain(j)*tmp_fact*obs_infl(j) !<- adaptive stop
          ENDIF
        ENDDO
      ENDDO
      ! xa mean
      hxm = hxm + KalmanGain*inov
      ! Perturbation Update
      DO j = 1,obsnum
        alpha_obs(j) = 1/(1+sqrt(tmp_R(i)*tmp_matrix2(j)))
        KalmanGain_new(j) = KalmanGain(j)*alpha_mdl(j)
      ENDDO
      DO j = 1,obsnum
        DO k = 1,PARA_DA_MEM
          hxp(j,k) = hxp(j,k) - KalmanGain_new(j)*hxo(k)
        ENDDO
      ENDDO
      DEALLOCATE(P,KalmanGain,KalmanGain_new)

    ENDDO


    !--- Finalize ---
    DO i = 1,PARA_DA_MEM
      xa (:,i) = xm  + xp (:,i)
      hxa(:,i) = hxm + hxp(:,i)
    ENDDO
    IF (PRESENT(out_hxa)) THEN
      out_hxa = hxa
    ENDIF


  return
  END SUBROUTINE Serial_EnSRF_CORE_tempered_v3

END MODULE Serial_EnSRF_tools
