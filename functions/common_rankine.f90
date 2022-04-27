MODULE common_rankine

!======================================================================================
  USE common
  IMPLICIT NONE
  PUBLIC

!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  !--- Parameter (DIRECTORY & FILE) ---
  CHARACTER(256),SAVE :: OUT_DIR
  CHARACTER(256),SAVE :: OBS_DIR
  CHARACTER(256),SAVE :: INI_DIR

  !--- Parameter (rankine vortex) ---
  INTEGER     ,SAVE   :: PARA_MEM, PARA_OBn
  REAL(r_size),SAVE   :: PARA_INF, PARA_LOC
  REAL(r_size),SAVE   :: PARA_Nef
  REAL(r_size),SAVE   :: PARA_RMm, PARA_RMs
  REAL(r_size),SAVE   :: PARA_VMm, PARA_VMs
  REAL(r_size),SAVE   :: PARA_POb, PARA_POs
  INTEGER     ,SAVE   :: PARA_Mra
  REAL(r_size),SAVE   :: PARA_OBs
  INTEGER     ,SAVE   :: PARA_RAf
  INTEGER     ,SAVE   :: PARA_IOt, PARA_JOt
  
  !--- Parameter (DA) ---
  INTEGER,SAVE        :: PARA_L96_J
  INTEGER,SAVE        :: PARA_DA_MEM
  REAL(r_size),SAVE   :: PARA_DA_LPF_Neff,PARA_DA_LPF_alpha
  INTEGER,SAVE        :: PARA_obs_num
  CHARACTER(256),SAVE :: MINRES_NCDIR

  !--- FLAG ---
  INTEGER,SAVE        :: PARA_LOOP_NUM
  LOGICAL,SAVE        :: FLAG_READ_MINRES

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_rankine
  IMPLICIT NONE
  INTEGER :: i
  namelist / DIR_FILE_list         / OUT_DIR,OBS_DIR,INI_DIR
  namelist / rankine_para_list     / PARA_MEM,PARA_INF,PARA_LOC,PARA_OBn, &
                                   & PARA_Nef,PARA_RMm,PARA_RMs, &
                                   & PARA_VMm,PARA_Vms,          &
                                   & PARA_POb,PARA_POs,          &
                                   & PARA_Mra,PARA_OBs,PARA_RAf, &
                                   & PARA_IOt,PARA_JOt,          &
                                   & PARA_LOOP_NUM,              &
                                   & FLAG_READ_MINRES, MINRES_NCDIR


  !-----------------------------------------------------------------------
  !--- READ NAMELIST ---
  OPEN(1,FILE='rankine_para.cnf')
  READ(1,NML=DIR_FILE_list)
  READ(1,NML=rankine_para_list)
  CLOSE(1)
!  PARA_POs = max(PARA_POs,0.001d0) ! <-- KK

  !--- DA PARA ---
  PARA_DA_MEM       = PARA_MEM
  PARA_DA_LPF_Neff  = PARA_Nef*dble(PARA_MEM)
  PARA_DA_LPF_alpha = PARA_INF
  PARA_obs_num      = PARA_OBn

  RETURN
END SUBROUTINE set_common_rankine

!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_para
  !- monitor -
  PRINT '(A)',' '
  PRINT '(A)','--------------------------------'
  PRINT '(A)','Rankine Vortex'
  PRINT '(A,I8)'  ,' PARA_MEM = ',PARA_MEM
  PRINT '(A,F8.3)',' PARA_INF = ',PARA_INF
  PRINT '(A,F8.3)',' PARA_LOC = ',PARA_LOC
  PRINT '(A,I8)'  ,' PARA_OBn = ',PARA_OBn
  PRINT '(A,F8.3)',' PARA_Nef = ',PARA_Nef
  PRINT '(A,F8.3)',' PARA_RMm = ',PARA_RMm
  PRINT '(A,F8.3)',' PARA_RMs = ',PARA_RMs
  PRINT '(A,F8.3)',' PARA_VMm = ',PARA_VMm
  PRINT '(A,F8.3)',' PARA_VMs = ',PARA_VMs
  PRINT '(A,F8.3)',' PARA_POb = ',PARA_POb
  PRINT '(A,F8.3)',' PARA_POs = ',PARA_POs
  PRINT '(A,I8)'  ,' PARA_Mra = ',PARA_Mra
  PRINT '(A,F8.3)',' PARA_OBs = ',PARA_OBs
  PRINT '(A,I8)'  ,' PARA_RAf = ',PARA_RAf
  PRINT '(A,I8)'  ,' PARA_IOt = ',PARA_IOt
  PRINT '(A,I8)'  ,' PARA_JOt = ',PARA_JOt
  
  RETURN
END SUBROUTINE monit_para

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

SUBROUTINE calc_vorticity(u,v,dim1,dim2,dx,vort)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN)    :: u   (dim1,dim2)
  REAL(r_size),INTENT(IN)    :: v   (dim1,dim2)
  REAL(r_size),INTENT(OUT)   :: vort(dim1,dim2)
  INTEGER     ,INTENT(IN)    :: dim1,dim2
  REAL(r_size),INTENT(IN)    :: dx
  REAL(r_size) :: dvx, duy
  INTEGER :: i,j
  DO j = 1,dim2
    DO i = 1,dim1
      IF     (i==1)      THEN; dvx = ( v(i+1,j) - v(i,j)   )/dx
      ELSEIF (i==2)      THEN; dvx = ( v(i+1,j) - v(i-1,j) )/2.0d0/dx
      ELSEIF (i==dim1-1) THEN; dvx = ( v(i+1,j) - v(i-1,j) )/2.0d0/dx
      ELSEIF (i==dim1)   THEN; dvx = ( v(i,j)   - v(i-1,j) )/dx
      ELSE;                    dvx = ( -v(i+2,j) + 8.0d0*v(i+1,j) - 8.0d0*v(i-1,j) + v(i-2,j) )/12.0d0/dx
      ENDIF

      IF     (j==1)      THEN; duy = ( u(i,j+1) - u(i,j)   )/dx
      ELSEIF (j==2)      THEN; duy = ( u(i,j+1) - u(i,j-1) )/2.0d0/dx
      ELSEIF (j==dim2-1) THEN; duy = ( u(i,j+1) - u(i,j-1) )/2.0d0/dx
      ELSEIF (j==dim2)   THEN; duy = ( u(i,j)   - u(i,j-1) )/dx
      ELSE;                    duy = ( -u(i,j+2) + 8.0d0*u(i,j+1) - 8.0d0*u(i,j-1) + u(i,j-2) )/12.0d0/dx
      ENDIF
  
      vort(i,j) = dvx - duy
    ENDDO
  ENDDO
return
END SUBROUTINE calc_vorticity
END MODULE common_rankine
