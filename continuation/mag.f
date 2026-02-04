C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
C  mag  : Non-Dimensional Formulation of a Cosserat rod under 
C         uniform dimensionless load parameter in a magnetic field
C         from a given homoclinic solution using Euler Parameters.
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
C 
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
C     ---------- ---- 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION U(NDIM),PAR(*),F(NDIM)
C
	F1  = U(1)
	F2  = U(2)
	F3  = U(3)
C
	XM1 = U(4)
	XM2 = U(5)
	XM3 = U(6)
C
	E31 = U(7)
	E32 = U(8)
	E33 = U(9)
C
	Q1  = U(10)
	Q2  = U(11)
	Q3  = U(12)
	Q4  = U(13)
C
	X   = U(14)
	Y   = U(15)
	Z   = U(16)
C
	xnu      = PAR(1) 
	xrho     = PAR(2)
	xm       = PAR(3)
	xlambda  = PAR(4)
	xepsilon = PAR(5)
	xT       = PAR(6)
	xD       = PAR(7)
	xR       = PAR(8)
C
C	     The Governing equations
C
C i) Force Balance
C
	F(1)  = (1.d0+xnu)*F2*XM3 - F3*XM2 +
     & xlambda*xlambda*(1.d0+xepsilon*F3)*E32
	F(2)  = (1.d0+xrho)*F3*XM1 - (1.d0+xnu)*F1*XM3 -
     & xlambda*(1.d0+xepsilon*F3)*E31
	F(3)  = F1*XM2 - (1.d0+xrho)*F2*XM1 
C
C ii) Moment Balance 
C
	F(4)  = xnu*XM2*XM3 + F2/(xm**2.d0)
	F(5)  = (1.d0+xrho)*XM1*XM3 - (1.d0+xnu)*XM1*XM3 - F1/(xm**2.d0)
	F(6)  = -xrho*XM1*XM2
C
C iii) Magnetic field
C
	F(7)  = (1.d0+xnu)*E32*XM3 - E33*XM2
	F(8)  = (1.d0+xrho)*E33*XM1 - (1.d0+xnu)*E31*XM3
	F(9)  = E31*XM2 - (1.d0+xrho)*E32*XM1
C
C iv) EulerParameter
C
	F(7)  = (1.d0/2.d0)*( +(1.d0+xrho)*XM1*Q4 - XM2*Q3 +
     & (1.d0+xnu)*XM3*Q2 )
	F(8)  = (1.d0/2.d0)*( +(1.d0+xrho)*XM1*Q3 + XM2*Q4 -
     & (1.d0+xnu)*XM3*Q1 )
	F(9)  = (1.d0/2.d0)*( -(1.d0+xrho)*XM1*Q2 + XM2*Q1 +
     & (1.d0+xnu)*XM3*Q4 )
	F(10) = (1.d0/2.d0)*( -(1.d0+xrho)*XM1*Q1 - XM2*Q2 -
     & (1.d0+xnu)*XM3*Q3 )
C
C v) Centre line equations
C
	F(11) = 2.d0*(Q1*Q3 + Q2*Q4)
	F(12) = 2.d0*(Q2*Q3 - Q1*Q4)
	F(13) = Q3**2.d0 + Q4**2.d0 - Q1**2.d0 - Q2**2.d0
C
C vi) Realign the solution to the initial data
C
	do i = 1,NDIM
		F(i) = xT*F(i)
	end do
C
      RETURN 
      END 
C 
      SUBROUTINE STPNT(NDIM,U,PAR,S) 
C     ---------- ---- 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION U(NDIM),PAR(*) 
C
C Load initial values
C
	open(unit=1,file='parameters.dat')
	read(1,'(8e20.12)') (PAR(i), i=1,8)
	close(1)
C
      RETURN
      END
C
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
C     ---------- ----
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)       
      PARAMETER (n=9, nrdens=9, lrpar=8, lipar=8)                           
      PARAMETER (lwork=11*n+8*nrdens+20, liwork=nrdens+20)      
      DIMENSION PAR(*),ICP(*),U0(NDIM),U1(NDIM),FB(NBC)  
      DIMENSION y(n),work(lwork),iwork(liwork)
      DIMENSION rpar(lrpar),ipar(lipar)
      DIMENSION z(n,n),vr(n,n),vi(n,n),vrdum(n,n),vidum(n,n)
      DIMENSION wrdum(n),widum(n),wr(n),wi(n),fv1(n)
      DIMENSION V(n,n), xProjection(n,n)
      DIMENSION xnorm(n)      
      DIMENSION static(n),v1(n),v2(n),vecdum(n)
      INTEGER, DIMENSION(n)   :: iv1
C
C define constants
C
	matz = 1 !not zero as only eigenvalues then
	zero = 0.00000001
	ierr = 0
	pi   = 4.d0*datan(1.d0) 
C
C clear vectors
C
	do i = 1,n
		do j = 1,n
			vr(i,j)          = 0.d0
			vi(i,j)          = 0.d0
			vrdum(i,j)       = 0.d0
			vidum(i,j)       = 0.d0
			xProjection(i,j) = 0.d0
			z(i,j)           = 0.d0
			V(i,j)           = 0.d0
		end do
		fv1(i)        = 0.d0
		iv1(i)        = 0
		wr(i)         = 0.d0
		wi(i)         = 0.d0
		xnorm(i)      = 0.d0
		wrdum(i)      = 0.d0
		widum(i)      = 0.d0
		static(i)     = 0.d0
	end do
C
C load parameters
C
	xnu      = PAR(1)
	xrho     = PAR(2)
	xm       = PAR(3)
	xlambda  = PAR(4)
	xepsilon = PAR(5)
	xT       = PAR(6)
	xD       = PAR(7)
	xR       = PAR(8)
C
C Create Projection matrix - linearise about the trivial solution
C
	y1 = 0.d0
	y2 = 0.d0
	y3 = 1.d0
C
	y4 = 0.d0
	y5 = 0.d0
	y6 = 1.d0
C
	y7 = 0.d0
	y8 = 0.d0
	y9 = 1.d0
C
	vrdum(1,2) = (1.d0+xnu)*y6 
	vrdum(1,3) =-y5 
	vrdum(1,5) =-1.d0*y3
	vrdum(1,6) = (1.d0+xnu)*y2
	vrdum(1,8) = xlambda*(1.d0+xepsilon)
C
	vrdum(2,1) = -(1.d0+xnu)*y6 
	vrdum(2,3) = +(1.d0+xrho)*y4 
	vrdum(2,4) = +(1.d0+xrho)*y3 
	vrdum(2,6) = -(1.d0+xnu)*y1
	vrdum(2,7) = -xlambda*(1.d0+xepsilon)
C
	vrdum(3,1) = y5
	vrdum(3,2) =-(1.d0+xrho)*y4
	vrdum(3,4) =-(1.d0+xrho)*y2 
	vrdum(3,5) = y1
C
	vrdum(4,2) = 1.d0/xm**2.d0 
	vrdum(4,5) = xnu*y6
	vrdum(4,6) = xnu*y5
C
	vrdum(5,1) =-1.d0/xm**2.d0 
	vrdum(5,4) = (xrho - xnu)*y6
	vrdum(5,6) = (xrho - xnu)*y4
C
	vrdum(6,4) = -xrho*y5
	vrdum(6,5) = -xrho*y4
C
	vrdum(7,5) = (1.d0+xnu)
	vrdum(7,8) =-1.d0
C
	vrdum(8,4) = (1.d0+xrho)
	vrdum(8,7) = -(1.d0+xnu) 
C
C Take transpose of linearised matrix
C
	do i = 1,n
		do j= 1,n
			xProjection(j,i) = vrdum(i,j)
		end do
	end do
C
C Find eigenvalues and eigenvectors of transpose of linearised matrix
C
	call rg(n,n,xProjection,wr,wi,matz,z,iv1,fv1,ierr)
	if ( ierr.ne.0 ) then
		print*, "failure in eigenvalue, eigenvector solver"
	end if
	write(99,*) ( wr(i), i = 1,9 )
	write(99,*) ( wi(i), i = 1,9 )
C
C Create matrices vr, vi from output matrix z
C
	iflag = 0
	do i = 1,n
		if ( wi(i+iflag).eq.0.d0 ) then
			do j = 1,n
				vr(j,i+iflag) = z(j,i+iflag) 
				vi(j,i+iflag) = 0.d0
			end do
		else if ( abs(wi(i+iflag)).eq.abs(wi(i+1+iflag)) ) then
			do j = 1,n
				vr(j,i+iflag)   = z(j,i+iflag)
				vr(j,i+1+iflag) = z(j,i+iflag)
				vi(j,i+iflag)   = z(j,i+1+iflag)
				vi(j,i+1+iflag) =-z(j,i+1+iflag)
			end do
			iflag = 1+iflag
		end if
		if ( i+iflag.eq.n ) goto 10
	end do
C
C Normalize vectors
C
10	do i = 1,n
		do j = 1,n
			xnorm(i) = xnorm(i) + vr(j,i)**2.d0 + vi(j,i)**2.d0 
		end do
	end do
	do i = 1,n
		do j = 1,n
			vr(j,i) = vr(j,i)/dsqrt(xnorm(i))
			vi(j,i) = vi(j,i)/dsqrt(xnorm(i))
		end do
	end do
C
C Order the eigenvectors/values according to the size of the real part of the eigenvalue. 
C
	do i = 1,n-1
		do j = i+1,n
			if ( wr(i).gt.wr(j) ) then 
				wrdum(i) = wr(i)
				widum(i) = wi(i)
				wr(i)    = wr(j)
				wi(i)    = wi(j)
				wr(j)    = wrdum(i)
				wi(j)    = widum(i)
				do k = 1,n
					vrdum(k,i) = vr(k,i)
					vidum(k,i) = vi(k,i)
					vr(k,i)    = vr(k,j)
					vi(k,i)    = vi(k,j)
					vr(k,j)    = vrdum(k,i)
					vi(k,j)    = vidum(k,i)
				end do
			end if
		end do
	end do
C
C Use real and imaginary parts of complex eigenvectors
C
	do i = 1,n-1
		jflag = 0
		do j = 1,n
			if ( abs(vr(j,i)-vr(j,i+1)).lt.zero ) then
				jflag = jflag + 1
				if ( jflag.eq.n ) then 
					do k = 1,n
						vr(k,i+1) = vi(k,i)
					end do
				end if
			end if
		end do
	end do
C
C Let V contain the stable as well as the unstable 9D eigenvectors
C
	do i = 1,n
		do j = 1,n
			V(i,j) = vr(j,i)
		end do
	end do
C
C define trivial solution
C
	static(3) = 1.d0
	static(6) = 1.d0
	static(9) = 1.d0
C
C Code the boundary conditions
C
	do i = 1,NBC
		FB(i) = 0.d0
	end do
C
C Four Projection conditions
C
	do i = 1,7
		do j = 1,n
			FB(i) = FB(i) + ( U0(j) - static(j) )*V(i,j)
		end do
	end do
C
C symmetric section
	FB(8)  = U1(1)
	FB(9)  = U1(4)
	FB(10)  = U1(7)
! C
! C casimirs
! 	FB(5)  = U0(3)-1.d0
! 	FB(6)  = U0(6)-1.d0
! 	FB(7)  = U0(9)-1.d0
C
C euler parameters
	FB(11) = U0(10)
	FB(12) = U0(11)
	FB(13) = U0(12) 
	FB(14) = U0(13)-1.d0
C
C centreline
	FB(15) = U0(14)
	FB(16) = U0(15)
	FB(17) = U0(16) 
C
! 	write(99,*) FB
! C
! C Two additional parameters in order to free parameters containing the end displacements: the end-shortening and the end rotation.
! C
! 	short  = xT - U1(16)
! 	FB(18) = xD - short
! C
! 	D11  = U1(10)**2.d0 - U1(11)**2.d0 - U1(12)**2.d0 + U1(13)**2.d0
! 	D12  = 2.d0*( U1(10)*U1(11) + U1(12)*U1(13) )
! 	twopi = 2.d0*pi
! 	if ( D12.lt.-1.d0 ) D12 = -1.d0
! 	if ( D12.gt.+1.d0 ) D12 = +1.d0
! 	xxnorm = DSQRT( D11**2.d0 + D12**2.d0)
! 	D1_11 = D11/xxnorm
! 	D1_12 = D12/xxnorm
! 	if ( D1_12.gt.0.d0 ) then
! 		rot = ( dasin(D1_12) - (1.d0 + xnu)*xT )/twopi
! 		nsign = +1
! 	else
! 		rot = ( twopi - dasin(D1_12) - (1.d0 + xnu)*xT )/twopi
! 		nsign = -1
! 	end if
! 	rotate = rot - DINT(rot)
! 	if ( rotate.lt.0.d0 ) rotate = rotate + 1.d0
! 	FB(19) = xR - rotate
C
C
100	format(9e15.6)
C
      RETURN 
      END 
C 
      SUBROUTINE ICND 
      RETURN 
      END 
C 
      SUBROUTINE FOPT 
      RETURN 
      END 
C
      SUBROUTINE PVLS 
      RETURN 
      END
C