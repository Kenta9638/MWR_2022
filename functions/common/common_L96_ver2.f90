MODULE common_L96

  !
  ! 2020.05.09 by KK
  ! ver.2
  !
!=======================================================================
  USE common
  USE common_L96DA

  CONTAINS
!=======================================================================
! [0] RUN Lorenz96
!=======================================================================
SUBROUTINE RUN_L96_NL(tmp_ini,tmp_ntimes,tmp_out,tmp_dt)
  IMPLICIT NONE

  REAL(r_size),INTENT(INOUT) :: tmp_ini(PARA_L96_J)
  INTEGER     ,INTENT(IN)    :: tmp_ntimes
  REAL(r_size),INTENT(OUT)   :: tmp_out(PARA_L96_J,tmp_ntimes+1)
  REAL(r_size),INTENT(IN),OPTIONAL   :: tmp_dt

  !--- ONLY for the first step in SPINUP ---
  IF ((maxval(tmp_ini)==minval(tmp_ini)) .AND. (minval(tmp_ini)==0.0d0)) THEN
    tmp_ini = tmp_ini + PARA_L96_F
    tmp_ini(PARA_L96_J/2) = tmp_ini(PARA_L96_J/2) + PARA_L96_P
  ENDIF
  
  IF (PRESENT(tmp_dt)) THEN ! only for 4DVAR
    CALL JP_FUNC_RUN_L96_core(tmp_ntimes,tmp_ini,tmp_out,tmp_dt)
  ELSE
    CALL JP_FUNC_RUN_L96_core(tmp_ntimes,tmp_ini,tmp_out)
  ENDIF

  RETURN
END SUBROUTINE RUN_L96_NL

!=======================================================================
! [2] core part of Lorenz96
!=======================================================================
!--[2.1] NL ------------------------------------------------------------
SUBROUTINE JP_FUNC_RUN_L96_core(ntimes,xin,xout,in_dt)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ntimes
  REAL(r_size),INTENT(IN)  :: xin  (1:PARA_L96_J)
  REAL(r_size),INTENT(OUT) :: xout (1:PARA_L96_J,1:ntimes+1)
  REAL(r_size),INTENT(IN),OPTIONAL :: in_dt
  REAL(r_size) :: x (PARA_L96_J)
  REAL(r_size) :: x0(PARA_L96_J+3)
  REAL(r_size),DIMENSION(PARA_L96_J):: k1,k2,k3,k4,y1,y2,y3,y4
  REAL(r_size) :: dt,F
  INTEGER      :: Nx,Nx_new ! <-- CAUTION !!
  INTEGER      :: i,j,t
  
  IF (PRESENT(in_dt)) THEN
    dt = in_dt
  ELSE
    dt  = PARA_L96_dt
  ENDIF
  F   = PARA_L96_F
  Nx  = PARA_L96_J
  xout(:,1) = xin
  Nx_new = Nx+3

  DO t = 1,ntimes
    x = xout(:,t)

    ! Place variables in vectors
    CALL dxt(x     ,k1); k1 = dt*k1
    CALL dxt(x+k1/2,k2); k2 = dt*k2
    CALL dxt(x+k2/2,k3); k3 = dt*k3
    CALL dxt(x+k3  ,k4); k4 = dt*k4

    ! Update variables
    x0(1)              = x(Nx-1)
    x0(2)              = x(Nx  )
    x0(3:Nx_new-1)     = x(:   )
    x0(Nx_new)         = x(1   )
    y1 = x0(4:Nx_new)  
    y2 = x0(3:Nx_new-1)  
    y3 = x0(2:Nx_new-2)  
    y4 = x0(1:Nx_new-3)  
    x  = y2 + k1/6 + k2/3 + k3/3 + k4/6
    xout(:,t+1) = x
  ENDDO

END SUBROUTINE JP_FUNC_RUN_L96_core

SUBROUTINE dxt(xin,xout)
  IMPLICIT NONE

  REAL(r_size),INTENT(IN)  :: xin  (1:PARA_L96_J)
  REAL(r_size),INTENT(OUT) :: xout (1:PARA_L96_J)
  REAL(r_size) :: x0(PARA_L96_J+3)
  REAL(r_size),DIMENSION(PARA_L96_J):: y1,y2,y3,y4
  REAL(r_size) :: F
  INTEGER      :: Nx,Nx_new
  F   = PARA_L96_F
  Nx  = PARA_L96_J
  Nx_new = Nx+3

  ! Create buffer zones on x for periodic domain
  x0(1)              = xin(Nx-1)
  x0(2)              = xin(Nx  )
  x0(3:Nx_new-1)     = xin(:   )
  x0(Nx_new)         = xin(1   )

  ! Place variables in vectors
  y1 = x0(4:Nx_new)  
  y2 = x0(3:Nx_new-1)  
  y3 = x0(2:Nx_new-2)  
  y4 = x0(1:Nx_new-3)  

  xout = (y1 - y4)*y3 - y2 + F

END SUBROUTINE dxt

!--[9.1] TL ------------------------------------------------------------
SUBROUTINE TL_lorenz96_core2(ntimes,xtraj,xin,in_dt,xout)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ntimes
  REAL(r_size),INTENT(IN)  :: xtraj(1:PARA_L96_J,1:ntimes+1)
  REAL(r_size),INTENT(IN)  :: xin  (1:PARA_L96_J)
  REAL(r_size),INTENT(IN)  :: in_dt
  REAL(r_size),INTENT(OUT) :: xout (1:PARA_L96_J,1:ntimes+1)
  REAL(r_size) :: x (PARA_L96_J+3,ntimes+1)
  REAL(r_size) :: xb(PARA_L96_J+3,ntimes+1)
  REAL(r_size),DIMENSION(PARA_L96_J):: x1,x2,x3,y1,y2,y3,y4
  REAL(r_size) :: dt,F
  INTEGER      :: Nx,Nx_new ! <-- CAUTION !!
  INTEGER      :: i,j,t
  
  dt  = in_dt
  F   = PARA_L96_F
  Nx  = PARA_L96_J
  Nx_new = PARA_L96_J+3 ! <-- CAUTION !!

  ! Create buffer zones on x for periodic domain
  x(1,1)      = xin(Nx-1)
  x(2,1)      = xin(Nx)
  x(3:Nx+2,1) = xin(:)
  x(Nx+3,1)   = xin(1)

  xb(1,:)      = xtraj(Nx-1,:)
  xb(2,:)      = xtraj(Nx,:)
  xb(3:Nx+2,:) = xtraj(:,:)
  xb(Nx+3,:)   = xtraj(1,:)

  DO t = 1,ntimes
    ! Place variables in vectors
    y1 = x(4:Nx_new  ,t); x1 = xb(4:Nx_new  ,t);
    y2 = x(3:Nx_new-1,t); x2 = xb(2:Nx_new-2,t);
    y3 = x(2:Nx_new-2,t); x3 = xb(1:Nx_new-3,t);
    y4 = x(1:Nx_new-3,t);
    ! Update variables
    x(3:Nx_new-1,t+1) = y2+((x1-x3)*y3+x2*y1-y2-x2*y4)*dt*(6-3*dt+dt*dt-dt*dt*dt/4)/6

    ! Update buffer zones
    x(1:2    ,t+1) = x(Nx_new-2:Nx_new-1,t+1)
    x(Nx_new ,t+1) = x(3,t+1)
  ENDDO

  xout(:,:) = x(3:Nx+2,:)


END SUBROUTINE TL_lorenz96_core2

!--[9.2] AD ------------------------------------------------------------
SUBROUTINE AD_lorenz96_core2(ntimes,xtraj,xin,in_dt,xout)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ntimes
  REAL(r_size),INTENT(IN)  :: xtraj(1:PARA_L96_J,1:ntimes+1)
  REAL(r_size),INTENT(IN)  :: xin  (1:PARA_L96_J)
  REAL(r_size),INTENT(IN)  :: in_dt 
  REAL(r_size),INTENT(OUT) :: xout (1:PARA_L96_J,1:ntimes+1)
  REAL(r_size) :: tmp_x(PARA_L96_J+3,1)
  REAL(r_size) :: x    (PARA_L96_J+3,ntimes+1)
  REAL(r_size) :: xb   (PARA_L96_J+4,ntimes+1)
  REAL(r_size),DIMENSION(PARA_L96_J):: x1,x2,x3,x4,y1,y2,y3,y4
  REAL(r_size) :: dt,F
  INTEGER      :: Nx,Nx_new1,Nx_new2 ! <-- CAUTION !!
  INTEGER      :: i,j,t
  
  dt  = in_dt
  F   = PARA_L96_F
  Nx  = PARA_L96_J
  Nx_new1 = PARA_L96_J+3 ! <-- CAUTION !!
  Nx_new2 = PARA_L96_J+4 ! <-- CAUTION !!

  ! Create buffer zones on x for periodic domain
  tmp_x(1,1)      = xin(Nx)
  tmp_x(2:Nx+1,1) = xin(:)
  tmp_x(Nx+2,1)   = xin(1)
  tmp_x(Nx+3,1)   = xin(2)

  xb(1,:)      = xtraj(Nx-1,:)
  xb(2,:)      = xtraj(Nx,:)
  xb(3:Nx+2,:) = xtraj(:,:)
  xb(Nx+3,:)   = xtraj(1,:)
  xb(Nx+4,:)   = xtraj(2,:)

  ! Fill x with zeros from t = 1:T
  x = 0.0d0
  x(:,ntimes+1) = tmp_x(:,1)


  DO t = ntimes+1,2,-1
    ! Place variables and background states in vectors
    y1 = x(4:Nx_new1  ,t); x1 = xb(5:Nx_new2  ,t);
    y2 = x(3:Nx_new1-1,t); x2 = xb(4:Nx_new2-1,t);
    y3 = x(2:Nx_new1-2,t); x3 = xb(2:Nx_new2-3,t);
    y4 = x(1:Nx_new1-3,t); x4 = xb(1:Nx_new2-4,t);

    ! Update variables
    x(2:Nx_new1-2,t-1) = y3+(x4*y4-y3+(x1-x3)*y2-x2*y1)*dt*(6-3*dt+dt*dt-dt*dt*dt/4)/6

    ! Update buffer zones
    x(1,t-1) = x(Nx_new1-2,t-1)
    x(Nx_new1-1:Nx_new1,t-1) = x(2:3,t-1)
  ENDDO

  xout(:,:) = x(2:Nx+1,:)


END SUBROUTINE AD_lorenz96_core2

END MODULE common_L96
