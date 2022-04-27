PROGRAM rankine_main
!$ use omp_lib

!======================================================================================
  USE common
  USE common_mpi
  USE common_rankine
  USE common_mtx
  use f95_lapack
  USE LPF_tools
  USE SW_test
  USE kk_netcdf_common
  USE kk_netcdf_tools

  IMPLICIT NONE

  !--- rankine vortex ---
  INTEGER                                 :: IM, Nx
  INTEGER     ,ALLOCATABLE,DIMENSION(:,:) :: xxx, yyy 
  INTEGER     ,ALLOCATABLE,DIMENSION(:)   :: xx, yy 
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: dx, dy, dd
  INTEGER     ,ALLOCATABLE,DIMENSION(:)   :: dif
  INTEGER     ,ALLOCATABLE,DIMENSION(:)   :: i_work11, i_work12
  INTEGER     ,ALLOCATABLE,DIMENSION(:,:) :: i_work21
  REAL(r_size)                            :: r_work01, r_work02, r_work03
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: r_work11
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: r_work21, r_work22
  INTEGER     ,ALLOCATABLE,DIMENSION(:)   :: ind
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: vt, phi, u, v, vr
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: u_true, v_true, velc_true, y
  REAL(r_size)                            :: drx, dry, dV, dR
  REAL(r_size)                            :: RM, VM
  REAL(r_size)                            :: IO, JO
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: hx, x, R_mtx
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: xa1, xa2, xa3, xa4
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: min_res_mdl, min_res_obs
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: mpi_gather_min_res_mdl, mpi_gather_min_res_obs
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:) :: tmp_OUT_xf0 
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:) :: tmp_OUT_xa1, tmp_OUT_xa2, tmp_OUT_xa3, tmp_OUT_xa4
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:) :: OUT_xf0 
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:) :: OUT_xa1, OUT_xa2, OUT_xa3, OUT_xa4
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa1_u_rmse,tmp_xa1_v_rmse,tmp_xa1_velc_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa2_u_rmse,tmp_xa2_v_rmse,tmp_xa2_velc_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa3_u_rmse,tmp_xa3_v_rmse,tmp_xa3_velc_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa4_u_rmse,tmp_xa4_v_rmse,tmp_xa4_velc_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa1_u_sprd,tmp_xa1_v_sprd,tmp_xa1_velc_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa2_u_sprd,tmp_xa2_v_sprd,tmp_xa2_velc_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa3_u_sprd,tmp_xa3_v_sprd,tmp_xa3_velc_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa4_u_sprd,tmp_xa4_v_sprd,tmp_xa4_velc_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa1_u,tmp_xa2_u,tmp_xa3_u,tmp_xa4_u
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa1_v,tmp_xa2_v,tmp_xa3_v,tmp_xa4_v
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: tmp_xa1_velc,tmp_xa2_velc,tmp_xa3_velc,tmp_xa4_velc
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: OUT_xa1_u_rmse,OUT_xa2_u_rmse,OUT_xa3_u_rmse,OUT_xa4_u_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: OUT_xa1_u_sprd,OUT_xa2_u_sprd,OUT_xa3_u_sprd,OUT_xa4_u_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: OUT_xa1_v_rmse,OUT_xa2_v_rmse,OUT_xa3_v_rmse,OUT_xa4_v_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: OUT_xa1_v_sprd,OUT_xa2_v_sprd,OUT_xa3_v_sprd,OUT_xa4_v_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: OUT_xa1_velc_rmse,OUT_xa2_velc_rmse,OUT_xa3_velc_rmse,OUT_xa4_velc_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: OUT_xa1_velc_sprd,OUT_xa2_velc_sprd,OUT_xa3_velc_sprd,OUT_xa4_velc_sprd

  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: mpi_xa1_u_rmse,mpi_xa1_v_rmse,mpi_xa1_velc_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: mpi_xa2_u_rmse,mpi_xa2_v_rmse,mpi_xa2_velc_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: mpi_xa3_u_rmse,mpi_xa3_v_rmse,mpi_xa3_velc_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: mpi_xa4_u_rmse,mpi_xa4_v_rmse,mpi_xa4_velc_rmse
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: mpi_xa1_u_sprd,mpi_xa1_v_sprd,mpi_xa1_velc_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: mpi_xa2_u_sprd,mpi_xa2_v_sprd,mpi_xa2_velc_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: mpi_xa3_u_sprd,mpi_xa3_v_sprd,mpi_xa3_velc_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: mpi_xa4_u_sprd,mpi_xa4_v_sprd,mpi_xa4_velc_sprd

  INTEGER     ,ALLOCATABLE,DIMENSION(:,:) :: H, H_extd
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: HC, HCH, loc_fact, loc_fact_uv

  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:) :: reshape_u, reshape_v, reshape_velc, reshape_vort
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: velc, hx_velc, loc_velc
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: vort, hx_vort, loc_vort

  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: SW_u_sum, SW_v_sum, SW_velc_sum, SW_vort_sum
!  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: SW_u_freq, SW_v_freq, SW_velc_freq, SW_vort_freq
  !, SW_v_mean, SW_velc_mean, SW_vort_mean

  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: SW_u_sum_W, SW_v_sum_W, SW_velc_sum_W, SW_vort_sum_W
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: SW_u_sum_K, SW_v_sum_K, SW_velc_sum_K, SW_vort_sum_K
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:) :: SW_u_sum_Z, SW_v_sum_Z, SW_velc_sum_Z, SW_vort_sum_Z
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: SW_u_mean_W, SW_v_mean_W, SW_velc_mean_W, SW_vort_mean_W
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: SW_u_mean_K, SW_v_mean_K, SW_velc_mean_K, SW_vort_mean_K
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: SW_u_mean_Z, SW_v_mean_Z, SW_velc_mean_Z, SW_vort_mean_Z

  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: xa_u,    xa_u_mean,    xa_u_rmse,    xa_u_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: xa_v,    xa_v_mean,    xa_v_rmse,    xa_v_sprd
  REAL(r_size),ALLOCATABLE,DIMENSION(:)   :: xa_velc, xa_velc_mean, xa_velc_rmse, xa_velc_sprd

  REAL(r_size),ALLOCATABLE,DIMENSION(:)       :: tmp_mpi_save_randn_all ! <-- !!
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: mpi_save_randn_all     ! <-- !!
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: mpi_save_randn         ! <-- !!
  !--- OUT 
  ! u
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: OUT_u
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: OUT_u_mean_SWres, OUT_u_mean_W, OUT_u_mean_K, OUT_u_mean_Z
  ! v
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: OUT_v
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: OUT_v_mean_SWres, OUT_v_mean_W, OUT_v_mean_K, OUT_v_mean_Z
  ! velc
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: OUT_velc
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: OUT_velc_mean_SWres, OUT_velc_mean_W, OUT_velc_mean_K, OUT_velc_mean_Z
  ! vort
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: OUT_vort
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:)     :: OUT_vort_mean_SWres, OUT_vort_mean_W, OUT_vort_mean_K, OUT_vort_mean_Z
  INTEGER :: i, j, k, n
  INTEGER :: iii
  INTEGER :: counter1, counter2
  REAL(r_size) :: tmp_var
  
  REAL(r_sngl), DIMENSION(1)       :: dummy_dim1 = UNDEF
  REAL(r_sngl), DIMENSION(1,1)     :: dummy_dim2 = UNDEF
  REAL(r_sngl), DIMENSION(1,1,1)   :: dummy_dim3 = UNDEF
  REAL(r_sngl), DIMENSION(1,1,1,1) :: dummy_dim4 = UNDEF


  !-- OUT ---
  CHARACTER(256) :: OUT_FILE1,OUT_FILE2,OUT_FILE3,OUT_FILE4
  CHARACTER(256) :: OUT_FILE5,OUT_FILE6,OUT_FILE7,OUT_FILE8
  CHARACTER(256) :: OUT_FILE9,OUT_FILE10,OUT_FILE11,OUT_FILE12,OUT_FILE13,OUT_FILE14
  CHARACTER(256) :: OUT_FILE15
  CHARACTER(256) :: TEXT1,TEXT2,TEXT3,TEXT4,TEXT5,TEXT6,TEXT7,TEXT8
  CHARACTER(256) :: TEXT99
  CHARACTER(256) :: TRUTH_FILE,OBS_FILE1,OBS_FILE2,OBS_FILE3,INI_FILE
  
  !--- Filter Divergence ---
  LOGICAL :: FD_FLAG  
  INTEGER :: FD_COUNTER 

  INTEGER :: irec1,irec2,irec3
  REAL :: time1

  LOGICAL :: upper
  REAL(r_size) :: r8_normal_01_cdf_inverse, alnorm

  !--- read minres
  CHARACTER(256) :: INP_FILE1
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:,:) :: mpi_xf_u_all
  REAL(r_sngl),ALLOCATABLE,DIMENSION(:,:,:,:) :: mpi_xf_v_all
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: mpi_xf_u_all_reshape
  REAL(r_size),ALLOCATABLE,DIMENSION(:,:,:)   :: mpi_xf_v_all_reshape

  !--- mpi ---
  INTEGER :: ierr!, nij1max, mpi_str, mpi_end
!  INTEGER :: MPI_PARA_LOOP_NUM

!======================================================================================
  !-----------------------------------------------------------------------
  ! mpi (1)
  !-----------------------------------------------------------------------
  CALL initialize_mpi

  !-----------------------------------------------------------------------
  ! Read namelists 
  !-----------------------------------------------------------------------
  CALL set_common_rankine
  CALL monit_para

  !-----------------------------------------------------------------------
  ! mpi (2)
  !-----------------------------------------------------------------------
  CALL get_str_end_num_mpi(PARA_LOOP_NUM)
!!  i = MOD(PARA_LOOP_NUM,nprocs)
!!  nij1max = (PARA_LOOP_NUM - i)/nprocs + 1
!!  IF(myrank < i) THEN
!!    MPI_PARA_LOOP_NUM = nij1max
!!    mpi_str = 1+nij1max*myrank
!!  ELSE
!!    MPI_PARA_LOOP_NUM = nij1max - 1
!!    mpi_str = nij1max*i+1+(myrank-i)*(nij1max-1)
!!  ENDIF
!!  mpi_end = mpi_str + MPI_PARA_LOOP_NUM-1
!!  WRITE(6,'(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of loops: MPI_PARA_LOOP_NUM= ',MPI_PARA_LOOP_NUM
!!  WRITE(*,*) 'MYRANK ',myrank,' mpi_str= ',mpi_str,' mpi_end= ',mpi_end
!!  WRITE(*,*) '-------------------------------------------------------------------------'

  !--- randn (mpi_scatter easy) ---
  ALLOCATE(tmp_mpi_save_randn_all(PARA_LOOP_NUM*PARA_MEM*4))
  ALLOCATE(    mpi_save_randn_all(PARA_LOOP_NUM,PARA_MEM,4))
  CALL com_randn(PARA_LOOP_NUM*PARA_MEM*4,tmp_mpi_save_randn_all)
  COUNTER1 = 0
  DO i = 1,PARA_LOOP_NUM
  DO j = 1,PARA_MEM
  DO k = 1,4
    COUNTER1 = COUNTER1 + 1
    mpi_save_randn_all(i,j,k) = tmp_mpi_save_randn_all(COUNTER1)
  ENDDO
  ENDDO
  ENDDO
  ALLOCATE(mpi_save_randn(MPI_PARA_LOOP_NUM,PARA_MEM,4))
  mpi_save_randn(:,:,:) = mpi_save_randn_all(mpi_str:mpi_end,:,:)
  DEALLOCATE(tmp_mpi_save_randn_all,mpi_save_randn_all)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !-----------------------------------------------------------------------
  ! Setting up experiment
  !-----------------------------------------------------------------------
  !--- define grid index matrices ---
  IM = PARA_Mra*2+1
  Nx = IM*IM
  ALLOCATE(xxx(IM,IM),yyy(IM,IM),xx(Nx),yy(Nx))
  counter1 = 0
  DO i = 1,IM
    DO j = 1,IM
      counter1 = counter1 + 1
      xx(counter1) = i
      yy(counter1) = j
    ENDDO
  ENDDO

  !--- define H ---
  ALLOCATE(dif(Nx)); dif = -9999.999d0
  counter1 = 0
  counter2 = 0
  DO i = 1,IM
    DO j = 1,IM
      counter1 = counter1 + 1
      r_work01 = sqrt(dble((i-PARA_IOt)**2+(j-PARA_JOt)**2))
      IF (r_work01<=30.0d0) THEN
        counter2 = counter2 + 1
        dif(counter2) = counter1
      ENDIF
    ENDDO
  ENDDO
  ALLOCATE(i_work11(counter2))
  ALLOCATE(i_work12(PARA_OBn))
  ALLOCATE(ind     (PARA_OBn))
  i_work11 = dif(1:counter2)
  CALL com_randperm(counter2,PARA_OBn,i_work12)
  DO i = 1,PARA_OBn
    ind(i) = i_work11(i_work12(i))
  ENDDO
  ALLOCATE(H(PARA_OBn,Nx)); H = 0
  DO i = 1,PARA_OBn
    H(i,ind(i)) = 1
  ENDDO
  DEALLOCATE(dif,i_work11,i_work12,ind)

  !--- DA ---
  PARA_L96_J = Nx*2
  ALLOCATE(R_mtx(PARA_OBn,PARA_OBn)); R_mtx = 0.0d0
  DO i = 1,PARA_OBn
  DO j = 1,PARA_OBn
    R_mtx(i,j) = PARA_OBs**2
  ENDDO
  ENDDO

  ALLOCATE(u_true(Nx),v_true(Nx),velc_true(Nx))
  ALLOCATE(hx(PARA_OBn,PARA_MEM))
  ALLOCATE(x (Nx*2,PARA_MEM)) ! u and v
  ALLOCATE(velc        (Nx,PARA_MEM))
  ALLOCATE(vort        (Nx,PARA_MEM))

  ALLOCATE(tmp_OUT_xf0(Nx*2,PARA_MEM,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_OUT_xa1(Nx*2,PARA_MEM,MPI_PARA_LOOP_NUM)) ! EnKF
  ALLOCATE(tmp_OUT_xa2(Nx*2,PARA_MEM,MPI_PARA_LOOP_NUM)) ! LoPF
  ALLOCATE(tmp_OUT_xa3(Nx*2,PARA_MEM,MPI_PARA_LOOP_NUM)) ! fx05
  ALLOCATE(tmp_OUT_xa4(Nx*2,PARA_MEM,MPI_PARA_LOOP_NUM)) ! adap

  ALLOCATE(min_res_mdl(Nx      ,MPI_PARA_LOOP_NUM)) ! vel
  ALLOCATE(min_res_obs(PARA_OBn,MPI_PARA_LOOP_NUM)) ! vel
  ALLOCATE(tmp_xa1_u_rmse   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa2_u_rmse   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa3_u_rmse   (Nx,MPI_PARA_LOOP_NUM))  
  ALLOCATE(tmp_xa4_u_rmse   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa1_u_sprd   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa2_u_sprd   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa3_u_sprd   (Nx,MPI_PARA_LOOP_NUM))  
  ALLOCATE(tmp_xa4_u_sprd   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa1_v_rmse   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa2_v_rmse   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa3_v_rmse   (Nx,MPI_PARA_LOOP_NUM))  
  ALLOCATE(tmp_xa4_v_rmse   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa1_v_sprd   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa2_v_sprd   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa3_v_sprd   (Nx,MPI_PARA_LOOP_NUM))  
  ALLOCATE(tmp_xa4_v_sprd   (Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa1_velc_rmse(Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa2_velc_rmse(Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa3_velc_rmse(Nx,MPI_PARA_LOOP_NUM))  
  ALLOCATE(tmp_xa4_velc_rmse(Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa1_velc_sprd(Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa2_velc_sprd(Nx,MPI_PARA_LOOP_NUM)) 
  ALLOCATE(tmp_xa3_velc_sprd(Nx,MPI_PARA_LOOP_NUM))  
  ALLOCATE(tmp_xa4_velc_sprd(Nx,MPI_PARA_LOOP_NUM)) 
  
  IF (FLAG_READ_MINRES) THEN
    INP_FILE1 = trim(MINRES_NCDIR)//'/'//'OUT_xf.nc'
    write(*,*) INP_FILE1
    ALLOCATE(mpi_xf_u_all(IM,IM,PARA_MEM,PARA_LOOP_NUM))
    ALLOCATE(mpi_xf_v_all(IM,IM,PARA_MEM,PARA_LOOP_NUM))
    ALLOCATE(mpi_xf_u_all_reshape(MPI_PARA_LOOP_NUM,PARA_MEM,Nx))
    ALLOCATE(mpi_xf_v_all_reshape(MPI_PARA_LOOP_NUM,PARA_MEM,Nx))
    CALL kk_netcdf_get_var_single(trim(INP_FILE1),'xf_u',4,IM,IM,PARA_MEM,PARA_LOOP_NUM, &
                                  dummy_dim1,dummy_dim2,dummy_dim3,mpi_xf_u_all)
    CALL kk_netcdf_get_var_single(trim(INP_FILE1),'xf_v',4,IM,IM,PARA_MEM,PARA_LOOP_NUM, &
                                  dummy_dim1,dummy_dim2,dummy_dim3,mpi_xf_v_all)
    counter1 = 0
    DO i = 1,IM
      DO j = 1,IM
        counter1 = counter1 + 1
        DO k = 1,MPI_PARA_LOOP_NUM
        DO n = 1,PARA_MEM
          mpi_xf_u_all_reshape(k,n,counter1) = dble(mpi_xf_u_all(j,i,n,mpi_str+k-1)) ! <-- i,j !!
          mpi_xf_v_all_reshape(k,n,counter1) = dble(mpi_xf_v_all(j,i,n,mpi_str+k-1)) ! <-- i,j !!
        ENDDO
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(mpi_xf_u_all,mpi_xf_v_all)
!    write(*,*) 'minval', minval(mpi_xf_u_all_reshape)
!    write(*,*) 'minval', minval(mpi_xf_v_all_reshape)
!    write(*,*) 'maxval', maxval(mpi_xf_u_all_reshape)
!    write(*,*) 'maxval', maxval(mpi_xf_v_all_reshape)
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !---------------------------------------------------------------------------
  ! Main Loop
  !---------------------------------------------------------------------------
  DO iii = 1,MPI_PARA_LOOP_NUM
  
    !--- Generate truth and prior samples ---
    ALLOCATE(dx(Nx),dy(Nx),dd(Nx))
    ALLOCATE(vt(Nx))
    ALLOCATE(phi(Nx))
    ALLOCATE(u(Nx),v(Nx),vr(Nx))

    DO n = 0,PARA_MEM

      !--- FLAG_READ_MINRES ---
      IF ((FLAG_READ_MINRES) .AND. (n/=0)) THEN
        u(:) = mpi_xf_u_all_reshape(iii,n,:)
        v(:) = mpi_xf_v_all_reshape(iii,n,:)
        !- Calculate obs-space priors
        dx = dble(xx)-dble(PARA_IOt)
        dy = dble(yy)-dble(PARA_JOt)
        DO i = 1,Nx
          phi(i) = atan(dy(i)/dx(i))
          IF (dx(i)<0) THEN
            phi(i) = phi(i) + pi
          ENDIF
          vr(i) = u(i)*cos(phi(i))+v(i)*sin(phi(i)) 
          IF (isnan(vr(i))) vr(i) = 0.0d0
        ENDDO
        hx(:,n)       = matmul(dble(H),vr)
        x (   1:Nx,n) = u 
        x (Nx+1:  ,n) = v 

      !--- orginal ---  
      ELSE
        !- Random perturbations for parameters
        IF (n == 0) THEN
          drx = 0.0d0
          dry = 0.0d0
          dV  = 0.0d0
          dR  = 0.0d0
        ELSE
          ALLOCATE(r_work11(4)) 
       !   CALL com_randn(4,r_work11)
          r_work11 = mpi_save_randn(iii,n,:)  ! <-- MPI !!
          drx = PARA_POs*r_work11(1)+PARA_POb
          dry = PARA_POs*r_work11(2)+PARA_POb
          dV  = PARA_VMs*r_work11(3)
          dR  = PARA_RMs*r_work11(4)
          DEALLOCATE(r_work11)
        ENDIF
        RM = max(PARA_RMm+dR,1.0d0)
        VM = max(PARA_VMm+dV,3.0d0)
        IO = dble(PARA_Mra)+drx
        JO = dble(PARA_Mra)+dry
         
        dx = dble(xx)-IO-1.0d0
        dy = dble(yy)-JO-1.0d0
        dd = sqrt(dx*dx+dy*dy)
        DO i = 1,Nx
          IF (dd(i)<=RM) THEN
            vt(i) = VM*dd(i)/RM
          ENDIF
          IF (dd(i)>RM) THEN
            vt(i) = VM*RM/dd(i)
          ENDIF
        ENDDO
    
        !- Convert to Cartesian coordinates
        DO i = 1,Nx
          phi(i) = atan(dy(i)/dx(i))
          IF (dx(i)<0) THEN
            phi(i) = phi(i) + pi
          ENDIF
          u(i) = -vt(i)*sin(phi(i)) 
          v(i) =  vt(i)*cos(phi(i)) 
          IF (isnan(u(i))) u(i) = 0.0d0
          IF (isnan(v(i))) v(i) = 0.0d0
        ENDDO
        
        !- Calculate obs-space priors
        dx = dble(xx)-dble(PARA_IOt)
        dy = dble(yy)-dble(PARA_JOt)
        DO i = 1,Nx
          phi(i) = atan(dy(i)/dx(i))
          IF (dx(i)<0) THEN
            phi(i) = phi(i) + pi
          ENDIF
          vr(i) = u(i)*cos(phi(i))+v(i)*sin(phi(i)) 
          IF (isnan(vr(i))) vr(i) = 0.0d0
        ENDDO
    
        IF (n == 0) THEN ! Save truth and create obs when n == 0
          u_true = u
          v_true = v
          velc_true = sqrt(u_true**2+v_true**2)
          ALLOCATE(r_work11(PARA_OBn)) 
          CALL com_randn(PARA_OBn,r_work11)
          y = matmul(dble(H),vr)+r_work11*PARA_OBs
          DEALLOCATE(r_work11)
        ELSE
          hx(:,n)       = matmul(dble(H),vr)
          x (   1:Nx,n) = u 
          x (Nx+1:  ,n) = v 
        ENDIF

      ENDIF ! FLAG_READ_MINRES
    ENDDO ! ensemble loop
    tmp_OUT_xf0(:,:,iii) = x

    DEALLOCATE(dx,dy,dd)
    DEALLOCATE(vt)
    DEALLOCATE(phi)
    DEALLOCATE(u,v,vr)

    !---------------------
    ! Velocity & Vorticity
    !---------------------
    ALLOCATE(reshape_u   (IM,IM,PARA_MEM))
    ALLOCATE(reshape_v   (IM,IM,PARA_MEM))
    ALLOCATE(reshape_velc(IM,IM,PARA_MEM))
    ALLOCATE(reshape_vort(IM,IM,PARA_MEM))
    DO k = 1,PARA_MEM
      counter1 = 0
      DO i = 1,IM
        DO j = 1,IM
          counter1 = counter1 + 1
          reshape_u   (xx(counter1),yy(counter1),k) = x(   counter1,k)
          reshape_v   (xx(counter1),yy(counter1),k) = x(Nx+counter1,k)
          reshape_velc(xx(counter1),yy(counter1),k) = sqrt(x(counter1,k)**2+x(Nx+counter1,k)**2)
        ENDDO
      ENDDO
      CALL calc_vorticity(reshape_u(:,:,k),reshape_v(:,:,k),IM,IM,1.0d0,reshape_vort(:,:,k))
    ENDDO
  
    DO k = 1,PARA_MEM
      counter1 = 0
      DO i = 1,IM
        DO j = 1,IM
          counter1 = counter1 + 1
          velc(counter1,k) = reshape_velc(xx(counter1),yy(counter1),k)
          vort(counter1,k) = reshape_vort(xx(counter1),yy(counter1),k)
        ENDDO
      ENDDO
    ENDDO
  
    DEALLOCATE(reshape_u,reshape_v,reshape_velc,reshape_vort)
  
    !---------------------
    ! Generate localization matrix 
    !---------------------
    ALLOCATE(dx(Nx),dy(Nx),dd(Nx))
    ALLOCATE(r_work21(Nx,Nx)); r_work21 = -9999.999d0
    DO i = 1,Nx
      dx = dble(xx-xx(i))
      dy = dble(yy-yy(i))
      dd = sqrt(dx*dx+dy*dy)
      r_work21(:,i) = exp(-dd**2/(2*PARA_LOC**2)) ! loc_fact
    ENDDO
    
    ALLOCATE(loc_fact(Nx,Nx)); loc_fact = -9999.999d0
    loc_fact = r_work21
    DEALLOCATE(dx,dy,dd)
    DEALLOCATE(r_work21)
  
    !---------------------
    ! MAIN
    !---------------------
    !--- LOC MATRIX ---
    ALLOCATE(r_work22(PARA_OBn,Nx))     ; r_work22    = -9999.999d0 
    ALLOCATE(HC         (PARA_OBn,Nx*2)); HC          = -9999.999d0
    ALLOCATE(loc_fact_uv(Nx*2    ,Nx*2)); loc_fact_uv = -9999.999d0
    r_work22 = matmul(dble(H),loc_fact) ! HC
    HC(:,   1:Nx) = r_work22
    HC(:,Nx+1:  ) = r_work22
    loc_fact_uv(   1:Nx,   1:Nx) = loc_fact
    loc_fact_uv(   1:Nx,Nx+1:  ) = loc_fact
    loc_fact_uv(Nx+1:  ,   1:Nx) = loc_fact
    loc_fact_uv(Nx+1:  ,Nx+1:  ) = loc_fact
    DEALLOCATE(r_work22)
    !--- SET ---
    ALLOCATE(HCH(PARA_OBn,PARA_OBn)); HCH = -9999.999d0
    ALLOCATE(H_extd(PARA_OBn,Nx*2))
    H_extd(:,   1:Nx) = H
    H_extd(:,Nx+1:  ) = 0 ! H
    HCH = matmul(HC,transpose(dble(H_extd)))
  
    ALLOCATE(xa1(Nx*2,PARA_MEM)) ! EnKF
    ALLOCATE(xa2(Nx*2,PARA_MEM)) ! LoPF
    ALLOCATE(xa3(Nx*2,PARA_MEM)) ! fx05
    ALLOCATE(xa4(Nx*2,PARA_MEM)) ! adap
    write(*,*) '--------------'
    write(*,*) 'EnKF', iii, myrank
    CALL LPF_CORE(x,loc_fact_uv,y,R_mtx,H_extd,hx,HC,HCH,.false.,1.0d0,-9999, &  ! IN
                  xa1)                                                           ! OUT
    write(*,*) '--------------'
    write(*,*) 'LoPF', iii, myrank
    CALL LPF_CORE(x,loc_fact_uv,y,R_mtx,H_extd,hx,HC,HCH,.false.,0.0d0,-9999, &  ! IN
                  xa2)                                                           ! OUT
    write(*,*) '--------------'
    write(*,*) 'fx05', iii, myrank
    CALL LPF_CORE(x,loc_fact_uv,y,R_mtx,H_extd,hx,HC,HCH,.false.,0.5d0,-9999, &  ! IN
                  xa3)                                                           ! OUT
    write(*,*) '--------------'
    write(*,*) 'adap', iii, myrank
    CALL LPF_CORE(x,loc_fact_uv,y,R_mtx,H_extd,hx,HC,HCH,.true.,-9999.9d0,Nx,  &  ! IN
                  xa4,                                                         &  ! OUT
                  velc,hx,loc_fact,mpi_str+iii-1,                              &  ! IN
                  min_res_mdl(:,iii),min_res_obs(:,iii))                          ! OUT
  
    ALLOCATE(tmp_xa1_u   (Nx,PARA_MEM)); tmp_xa1_u = xa1(1:Nx,:)
    ALLOCATE(tmp_xa2_u   (Nx,PARA_MEM)); tmp_xa2_u = xa2(1:Nx,:)
    ALLOCATE(tmp_xa3_u   (Nx,PARA_MEM)); tmp_xa3_u = xa3(1:Nx,:)
    ALLOCATE(tmp_xa4_u   (Nx,PARA_MEM)); tmp_xa4_u = xa4(1:Nx,:)
    ALLOCATE(tmp_xa1_v   (Nx,PARA_MEM)); tmp_xa1_v = xa1(Nx+1:,:)
    ALLOCATE(tmp_xa2_v   (Nx,PARA_MEM)); tmp_xa2_v = xa2(Nx+1:,:)
    ALLOCATE(tmp_xa3_v   (Nx,PARA_MEM)); tmp_xa3_v = xa3(Nx+1:,:)
    ALLOCATE(tmp_xa4_v   (Nx,PARA_MEM)); tmp_xa4_v = xa4(Nx+1:,:)
    ALLOCATE(tmp_xa1_velc(Nx,PARA_MEM)); tmp_xa1_velc = sqrt(tmp_xa1_u**2+tmp_xa1_v**2)
    ALLOCATE(tmp_xa2_velc(Nx,PARA_MEM)); tmp_xa2_velc = sqrt(tmp_xa2_u**2+tmp_xa2_v**2)
    ALLOCATE(tmp_xa3_velc(Nx,PARA_MEM)); tmp_xa3_velc = sqrt(tmp_xa3_u**2+tmp_xa3_v**2)
    ALLOCATE(tmp_xa4_velc(Nx,PARA_MEM)); tmp_xa4_velc = sqrt(tmp_xa4_u**2+tmp_xa4_v**2)
  
    DO j = 1,Nx
      CALL com_rms (PARA_MEM,tmp_xa1_u(j,:)-u_true(j),tmp_xa1_u_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa2_u(j,:)-u_true(j),tmp_xa2_u_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa3_u(j,:)-u_true(j),tmp_xa3_u_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa4_u(j,:)-u_true(j),tmp_xa4_u_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa1_v(j,:)-v_true(j),tmp_xa1_v_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa2_v(j,:)-v_true(j),tmp_xa2_v_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa3_v(j,:)-v_true(j),tmp_xa3_v_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa4_v(j,:)-v_true(j),tmp_xa4_v_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa1_velc(j,:)-velc_true(j),tmp_xa1_velc_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa2_velc(j,:)-velc_true(j),tmp_xa2_velc_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa3_velc(j,:)-velc_true(j),tmp_xa3_velc_rmse(j,iii))
      CALL com_rms (PARA_MEM,tmp_xa4_velc(j,:)-velc_true(j),tmp_xa4_velc_rmse(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa1_u(j,:),tmp_xa1_u_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa2_u(j,:),tmp_xa2_u_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa3_u(j,:),tmp_xa3_u_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa4_u(j,:),tmp_xa4_u_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa1_v(j,:),tmp_xa1_v_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa2_v(j,:),tmp_xa2_v_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa3_v(j,:),tmp_xa3_v_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa4_v(j,:),tmp_xa4_v_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa1_velc(j,:),tmp_xa1_velc_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa2_velc(j,:),tmp_xa2_velc_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa3_velc(j,:),tmp_xa3_velc_sprd(j,iii))
      CALL com_sprd(1,PARA_MEM,tmp_xa4_velc(j,:),tmp_xa4_velc_sprd(j,iii))
    ENDDO
      
    tmp_OUT_xa1(:,:,iii) = xa1
    tmp_OUT_xa2(:,:,iii) = xa2
    tmp_OUT_xa3(:,:,iii) = xa3
    tmp_OUT_xa4(:,:,iii) = xa4
  
    DEALLOCATE(xa1,xa2,xa3,xa4)
    DEALLOCATE(tmp_xa1_u,tmp_xa2_u,tmp_xa3_u,tmp_xa4_u)
    DEALLOCATE(tmp_xa1_v,tmp_xa2_v,tmp_xa3_v,tmp_xa4_v)
    DEALLOCATE(tmp_xa1_velc,tmp_xa2_velc,tmp_xa3_velc,tmp_xa4_velc)
  
    DEALLOCATE(loc_fact)
    DEALLOCATE(HC,loc_fact_uv,HCH,H_extd)
  
  ENDDO ! iii
  IF (FLAG_READ_MINRES) DEALLOCATE(mpi_xf_u_all_reshape,mpi_xf_v_all_reshape)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !--------------------------------------------------------
  ! FINALIZE
  !--------------------------------------------------------
  IF (myrank == 0) THEN
    ALLOCATE(OUT_xf0(Nx*2,PARA_MEM,PARA_LOOP_NUM)) 
    ALLOCATE(OUT_xa1(Nx*2,PARA_MEM,PARA_LOOP_NUM)) ! EnKF
    ALLOCATE(OUT_xa2(Nx*2,PARA_MEM,PARA_LOOP_NUM)) ! LoPF
    ALLOCATE(OUT_xa3(Nx*2,PARA_MEM,PARA_LOOP_NUM)) ! fx05
    ALLOCATE(OUT_xa4(Nx*2,PARA_MEM,PARA_LOOP_NUM)) ! adap
    ALLOCATE(mpi_gather_min_res_mdl(Nx      ,PARA_LOOP_NUM))
    ALLOCATE(mpi_gather_min_res_obs(PARA_OBn,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa1_u_rmse   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa2_u_rmse   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa3_u_rmse   (Nx,PARA_LOOP_NUM))  
    ALLOCATE(mpi_xa4_u_rmse   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa1_u_sprd   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa2_u_sprd   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa3_u_sprd   (Nx,PARA_LOOP_NUM))  
    ALLOCATE(mpi_xa4_u_sprd   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa1_v_rmse   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa2_v_rmse   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa3_v_rmse   (Nx,PARA_LOOP_NUM))  
    ALLOCATE(mpi_xa4_v_rmse   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa1_v_sprd   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa2_v_sprd   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa3_v_sprd   (Nx,PARA_LOOP_NUM))  
    ALLOCATE(mpi_xa4_v_sprd   (Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa1_velc_rmse(Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa2_velc_rmse(Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa3_velc_rmse(Nx,PARA_LOOP_NUM))  
    ALLOCATE(mpi_xa4_velc_rmse(Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa1_velc_sprd(Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa2_velc_sprd(Nx,PARA_LOOP_NUM)) 
    ALLOCATE(mpi_xa3_velc_sprd(Nx,PARA_LOOP_NUM))  
    ALLOCATE(mpi_xa4_velc_sprd(Nx,PARA_LOOP_NUM)) 
  ENDIF
  CALL MPI_GATHER(tmp_OUT_xf0,Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  OUT_xf0(:,:,mpi_str:mpi_end),Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_OUT_xa1,Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  OUT_xa1(:,:,mpi_str:mpi_end),Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_OUT_xa2,Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  OUT_xa2(:,:,mpi_str:mpi_end),Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_OUT_xa3,Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  OUT_xa3(:,:,mpi_str:mpi_end),Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_OUT_xa4,Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  OUT_xa4(:,:,mpi_str:mpi_end),Nx*2*PARA_MEM*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(min_res_mdl,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_gather_min_res_mdl(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(min_res_obs,PARA_OBn*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_gather_min_res_obs(:,mpi_str:mpi_end),PARA_OBn*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa1_u_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa1_u_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa2_u_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa2_u_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa3_u_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa3_u_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa4_u_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa4_u_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa1_u_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa1_u_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa2_u_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa2_u_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa3_u_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa3_u_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa4_u_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa4_u_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa1_v_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa1_v_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa2_v_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa2_v_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa3_v_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa3_v_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa4_v_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa4_v_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa1_v_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa1_v_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa2_v_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa2_v_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa3_v_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa3_v_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa4_v_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa4_v_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa1_velc_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa1_velc_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa2_velc_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa2_velc_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa3_velc_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa3_velc_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa4_velc_rmse,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa4_velc_rmse(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa1_velc_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa1_velc_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa2_velc_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa2_velc_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa3_velc_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa3_velc_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(tmp_xa4_velc_sprd,Nx*MPI_PARA_LOOP_NUM,MPI_REAL8, &
                  mpi_xa4_velc_sprd(:,mpi_str:mpi_end),Nx*MPI_PARA_LOOP_NUM,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  IF (myrank == 0) THEN

    !--- true ---
    OUT_FILE1 = trim(OUT_DIR)//'/'//'OUT_true.nc'
    CALL kk_netcdf_define_dims   (trim(OUT_FILE1),4,'xdim',IM,'ydim',IM,'Mem',PARA_MEM,'repetition',PARA_LOOP_NUM)
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'u_true','u_true','xdim','ydim','UNDEF','UNDEF', &
                                  IM,IM,9999,9999,dummy_dim1,sngl(u_true),dummy_dim3,dummy_dim4)
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'v_true','v_true','xdim','ydim','UNDEF','UNDEF', &
                                  IM,IM,9999,9999,dummy_dim1,sngl(v_true),dummy_dim3,dummy_dim4)
    !--- xf ---
    OUT_FILE1 = trim(OUT_DIR)//'/'//'OUT_xf.nc'
    CALL kk_netcdf_define_dims   (trim(OUT_FILE1),4,'xdim',IM,'ydim',IM,'Mem',PARA_MEM,'repetition',PARA_LOOP_NUM)
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xf_u','xf_u','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xf0(1:Nx,:,:)))
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xf_v','xf_v','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xf0(Nx+1:,:,:)))
    !--- xa ---
    OUT_FILE1 = trim(OUT_DIR)//'/'//'OUT_xa1.nc'
    CALL kk_netcdf_define_dims   (trim(OUT_FILE1),4,'xdim',IM,'ydim',IM,'Mem',PARA_MEM,'repetition',PARA_LOOP_NUM)
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xa1_u','xa_u(EnKF)','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xa1(1:Nx,:,:)))
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xa1_v','xa_v(EnKF)','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xa1(Nx+1:,:,:)))
    OUT_FILE1 = trim(OUT_DIR)//'/'//'OUT_xa2.nc'
    CALL kk_netcdf_define_dims   (trim(OUT_FILE1),4,'xdim',IM,'ydim',IM,'Mem',PARA_MEM,'repetition',PARA_LOOP_NUM)
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xa2_u','xa_u(LPF)','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xa2(1:Nx,:,:)))
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xa2_v','xa_v(LPF)','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xa2(Nx+1:,:,:)))
    OUT_FILE1 = trim(OUT_DIR)//'/'//'OUT_xa3.nc'
    CALL kk_netcdf_define_dims   (trim(OUT_FILE1),4,'xdim',IM,'ydim',IM,'Mem',PARA_MEM,'repetition',PARA_LOOP_NUM)
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xa3_u','xa_u(Hyb:0.5)','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xa3(1:Nx,:,:)))
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xa3_v','xa_v(Hyb:0.5)','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xa3(Nx+1:,:,:)))
    OUT_FILE1 = trim(OUT_DIR)//'/'//'OUT_xa4.nc'
    CALL kk_netcdf_define_dims   (trim(OUT_FILE1),4,'xdim',IM,'ydim',IM,'Mem',PARA_MEM,'repetition',PARA_LOOP_NUM)
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xa4_u','xa_u(Hyb:adapt)','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xa4(1:Nx,:,:)))
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),4,'xa4_v','xa_v(Hyb:adapt)','xdim','ydim','Mem','repetition', &
                                  IM,IM,PARA_MEM,PARA_LOOP_NUM,dummy_dim1,dummy_dim2,dummy_dim3,sngl(OUT_xa4(Nx+1:,:,:)))
    !--- minres ---
    OUT_FILE1 = trim(OUT_DIR)//'/'//'OUT_minres.nc'
    CALL kk_netcdf_define_dims(trim(OUT_FILE1),5,'xdim',IM,'ydim',IM,'Mem',PARA_MEM,'repetition',PARA_LOOP_NUM,'obs_num',PARA_OBn)
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),3,'min_res_mdl','min_res_mdl','xdim','ydim','repetition','UNDEF', &
                                  IM,IM,PARA_LOOP_NUM,9999,dummy_dim1,dummy_dim2,sngl(mpi_gather_min_res_mdl),dummy_dim4)
    CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'min_res_obs','min_res_obs','obs_num','repetition','UNDEF','UNDEF', &
                                  PARA_OBn,PARA_LOOP_NUM,9999,9999,dummy_dim1,sngl(mpi_gather_min_res_obs),dummy_dim3,dummy_dim4)
    !--- score ---
    ALLOCATE(OUT_xa1_u_rmse(Nx),OUT_xa1_v_rmse(Nx),OUT_xa1_velc_rmse(Nx))
    ALLOCATE(OUT_xa2_u_rmse(Nx),OUT_xa2_v_rmse(Nx),OUT_xa2_velc_rmse(Nx))
    ALLOCATE(OUT_xa3_u_rmse(Nx),OUT_xa3_v_rmse(Nx),OUT_xa3_velc_rmse(Nx))
    ALLOCATE(OUT_xa4_u_rmse(Nx),OUT_xa4_v_rmse(Nx),OUT_xa4_velc_rmse(Nx))
    ALLOCATE(OUT_xa1_u_sprd(Nx),OUT_xa1_v_sprd(Nx),OUT_xa1_velc_sprd(Nx))
    ALLOCATE(OUT_xa2_u_sprd(Nx),OUT_xa2_v_sprd(Nx),OUT_xa2_velc_sprd(Nx))
    ALLOCATE(OUT_xa3_u_sprd(Nx),OUT_xa3_v_sprd(Nx),OUT_xa3_velc_sprd(Nx))
    ALLOCATE(OUT_xa4_u_sprd(Nx),OUT_xa4_v_sprd(Nx),OUT_xa4_velc_sprd(Nx))
    DO j = 1,Nx
      CALL com_mean(PARA_LOOP_NUM,mpi_xa1_u_rmse(j,:),OUT_xa1_u_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa2_u_rmse(j,:),OUT_xa2_u_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa3_u_rmse(j,:),OUT_xa3_u_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa4_u_rmse(j,:),OUT_xa4_u_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa1_v_rmse(j,:),OUT_xa1_v_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa2_v_rmse(j,:),OUT_xa2_v_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa3_v_rmse(j,:),OUT_xa3_v_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa4_v_rmse(j,:),OUT_xa4_v_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa1_velc_rmse(j,:),OUT_xa1_velc_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa2_velc_rmse(j,:),OUT_xa2_velc_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa3_velc_rmse(j,:),OUT_xa3_velc_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa4_velc_rmse(j,:),OUT_xa4_velc_rmse(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa1_u_sprd(j,:),OUT_xa1_u_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa2_u_sprd(j,:),OUT_xa2_u_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa3_u_sprd(j,:),OUT_xa3_u_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa4_u_sprd(j,:),OUT_xa4_u_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa1_v_sprd(j,:),OUT_xa1_v_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa2_v_sprd(j,:),OUT_xa2_v_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa3_v_sprd(j,:),OUT_xa3_v_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa4_v_sprd(j,:),OUT_xa4_v_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa1_velc_sprd(j,:),OUT_xa1_velc_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa2_velc_sprd(j,:),OUT_xa2_velc_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa3_velc_sprd(j,:),OUT_xa3_velc_sprd(j))
      CALL com_mean(PARA_LOOP_NUM,mpi_xa4_velc_sprd(j,:),OUT_xa4_velc_sprd(j))
    ENDDO

OUT_FILE1 = trim(OUT_DIR)//'/'//'OUT_score.nc'
CALL kk_netcdf_define_dims   (trim(OUT_FILE1),2,'xdim',IM,'ydim',IM)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa1_u_rmse','xa1_u_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa1_u_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa2_u_rmse','xa2_u_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa2_u_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa3_u_rmse','xa3_u_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa3_u_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa4_u_rmse','xa4_u_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa4_u_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa1_v_rmse','xa1_v_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa1_v_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa2_v_rmse','xa2_v_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa2_v_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa3_v_rmse','xa3_v_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa3_v_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa4_v_rmse','xa4_v_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa4_v_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa1_velc_rmse','xa1_velc_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa1_velc_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa2_velc_rmse','xa2_velc_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa2_velc_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa3_velc_rmse','xa3_velc_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa3_velc_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa4_velc_rmse','xa4_velc_rmse','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa4_velc_rmse),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa1_u_sprd','xa1_u_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa1_u_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa2_u_sprd','xa2_u_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa2_u_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa3_u_sprd','xa3_u_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa3_u_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa4_u_sprd','xa4_u_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa4_u_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa1_v_sprd','xa1_v_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa1_v_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa2_v_sprd','xa2_v_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa2_v_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa3_v_sprd','xa3_v_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa3_v_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa4_v_sprd','xa4_v_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa4_v_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa1_velc_sprd','xa1_velc_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa1_velc_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa2_velc_sprd','xa2_velc_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa2_velc_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa3_velc_sprd','xa3_velc_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa3_velc_sprd),dummy_dim3,dummy_dim4)
CALL kk_netcdf_put_var_single(trim(OUT_FILE1),2,'xa4_velc_sprd','xa4_velc_sprd','xdim','ydim','UNDEF','UNDEF', &
                              IM,IM,9999,9999,dummy_dim1,sngl(OUT_xa4_velc_sprd),dummy_dim3,dummy_dim4)

  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi
STOP



END PROGRAM rankine_main
