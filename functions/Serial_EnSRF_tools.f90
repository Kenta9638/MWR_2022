MODULE Serial_EnSRF_tools
  USE common
  USE common_mtx
  USE common_rankine
  use f95_lapack

  IMPLICIT NONE
  PRIVATE
!  PUBLIC :: Serial_EnSRF_CORE, Serial_EnSRF_CORE_tempered
  PUBLIC :: Serial_EnSRF_CORE_tempered

  CONTAINS

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

END MODULE Serial_EnSRF_tools
