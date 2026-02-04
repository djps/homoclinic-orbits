C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
C  ext :  Non-Dimensional Formulation of a Cosserat rod under 
C         uniform dimensionless load parameter from a given homoclinic 
C         solution using Euler Parameters.
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
	Q1  = U(7)
	Q2  = U(8)
	Q3  = U(9)
	Q4  = U(10)
C
	X   = U(11)
	Y   = U(12)
	Z   = U(13)
C
	xnu      = PAR(1) 
	xrho     = PAR(2)
	xm       = PAR(3)
	xgamma   = PAR(4)
	xsigma   = PAR(5)
	xepsilon = PAR(6)
	xT       = PAR(7)
	xD       = PAR(8)
	xR       = PAR(9)
C 
	F(1)  = (1.d0+xnu)*F2*XM3 - F3*XM2
	F(2)  = (1.d0+xrho)*F3*XM1 - (1.d0+xnu)*F1*XM3
	F(3)  = F1*XM2 - (1.d0+xrho)*F2*XM1
C
	F(4)  = xnu*XM2*XM3 + F2*(1.d0+xepsilon*F3)/(xm**2.d0) 
	F(5)  = (xrho-xnu)*XM1*XM3 - F1*(1.d0+xepsilon*F3)/(xm**2.d0)  
	F(6)  =-xrho*XM1*XM2
C
	F(7)  = 0.5*(+(1.d0+xrho)*XM1*Q4-XM2*Q3+(1.d0+xnu)*XM3*Q2)
	F(8)  = 0.5*(+(1.d0+xrho)*XM1*Q3+XM2*Q4-(1.d0+xnu)*XM3*Q1)
	F(9)  = 0.5*(-(1.d0+xrho)*XM1*Q2+XM2*Q1+(1.d0+xnu)*XM3*Q4)
	F(10) = 0.5*(-(1.d0+xrho)*XM1*Q1-XM2*Q2-(1.d0+xnu)*XM3*Q3)
C
	F(11) = (1.d0+xepsilon*F3)*2.d0*(Q1*Q3 + Q2*Q4)
C
	F(12) = (1.d0+xepsilon*F3)*2.d0*(Q2*Q3 - Q1*Q4)
C
	F(13) = (1.d0+xepsilon*F3)*(Q3**2.d0+Q4**2.d0-Q1**2.d0-Q2**2.d0)
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
	open(unit=1,file='parameters.dat')
	read(1,*) (PAR(i), i=1,9)
	close(1)
C
      
      RETURN
      END
C
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
C     ---------- ----
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)       
      PARAMETER (n=6)                                
      DIMENSION PAR(*),ICP(*),U0(NDIM),U1(NDIM),FB(NBC)  
      DIMENSION z(n,n),vr(n,n),vi(n,n),vrdum(n,n),vidum(n,n)
      DIMENSION wrdum(n),widum(n),wr(n),wi(n),fv1(n),iv1(n)
      DIMENSION xnorm(n), V(n,n), c(n,n), a(n,n), xI(4), xI1(4)      
C
	xnu      = PAR(1) 
	xrho     = PAR(2)
	xm       = PAR(3)
	xgamma   = PAR(4)
	xsigma   = PAR(5)
	xepsilon = PAR(6)
	xT       = PAR(7)
	xD       = PAR(8)
	xR       = PAR(9)
C
	matz = 1
	zero = 0.00000001
	pi   = 4.d0*datan(1.d0)
	do i = 1,n
		wrdum(i) = 0.d0
		widum(i) = 0.d0
		xnorm(i) = 0.d0
		do j = 1,n
			c(i,j)     = 0.d0
			a(i,j)     = 0.d0
			vr(i,j)    = 0.d0
			vi(i,j)    = 0.d0
			vrdum(i,j) = 0.d0
			vidum(i,j) = 0.d0
		end do
	end do
C
	c(1,2) = (1.d0+xnu)
	c(1,5) =-1.d0
	c(2,1) =-(1.d0+xnu)
	c(2,4) = (1.d0+xrho)
	c(4,2) = (1.d0+xepsilon)/(xm**2.d0)
	c(4,5) = xnu
	c(5,1) =-(1.d0+xepsilon)/(xm**2.d0)
	c(5,4) = xrho - xnu
C
	do i = 1,n
		do j = 1,n
			a(i,j) = c(j,i)
		end do
	end do
C
	call rg(n,n,a,wr,wi,matz,z,iv1,fv1,ierr)
	if ( ierr.ne.0 ) then
		print*, "failure in eigenvalue, eigenvector solver"
	end if
C
C Convert output z into matrices of real and imaginary parts of eigenvectors.
C
	iflag = 0
	do i = 1,n
		if ( wi(i+iflag).eq.0.d0 ) then
			do j = 1,n
				vr(j,i+iflag) = z(j,i+iflag) 
				vi(j,i+iflag) = 0.d0
			end do
		else if ( abs(wi(i+iflag)).eq.abs(wi(i+iflag+1)) ) then
			do j = 1,n
				vr(j,i+iflag)   = z(j,i+iflag)
				vr(j,i+iflag+1) = z(j,i+iflag)
				vi(j,i+iflag)   = z(j,i+iflag+1)
				vi(j,i+iflag+1) =-z(j,i+iflag+1)
			end do
			iflag = 1 + iflag
		end if
		if ( i+iflag.eq.n ) goto 18
	end do
C
C Normalize vectors
C
18	do i = 1,n
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
C Order the eigenvectors/values according to the size of the real part of the
C eigenvalue. 
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
C Let V contain the stable as well as the unstable 6D eigenvectors
C
	do i = 1,n
		do j = 1,n
			V(i,j) = 0.d0
		end do 
	end do
	do i = 1,n
		do j = 1,n
			V(i,j) = vr(j,i)
		end do
	end do
C
C Code the boundary conditions
C
	do i = 1,NBC
		FB(i) = 0.d0
	end do
C
C Two Projection conditions
C
	do i = 1,2
		do j = 1,n
			FB(i) = FB(i) + U0(j)*V(i,j)
		end do
	end do
C
C Integral fixing conditions
C
	FB(3) = U0(3) - 1.d0
	FB(4) = U0(6) - 1.d0
C
C Two end point conditions in the symmetry plane, corresponding to specific
C reversibility
C
	FB(5) = U1(1)
	FB(6) = U1(4)
C
C The Postprocessed Euler parameters
C
	FB(7)  = U0(7)
	FB(8)  = U0(8)
	FB(9)  = U0(9)
	FB(10) = U0(10) - 1.d0
C
C The Centreline Equation
C
	FB(11) = U0(11)
	FB(12) = U0(12)
	FB(13) = U0(13)
C
C Two additional parameters in order to free parameters containing the end
C displacements: the end-shortening and the end rotation.
C
	short  = 2.d0*( xT*(1.d0+xepsilon) - U1(13) )
	FB(14) = xD - short
C
	D11  = U1(7)**2.d0 - U1(8)**2.d0 - U1(9)**2.d0 + U1(10)**2.d0
	D12  = 2.d0*( U1(7)*U1(8) + U1(9)*U1(10) )
	twopi = 2.d0*pi
	if ( D11.lt.-1.d0 ) D11 = -1.d0
	if ( D11.gt.+1.d0 ) D11 = +1.d0
	xxnorm = DSQRT( D11**2.d0 + D12**2.d0)
	D1_11 = D11/xxnorm
	D1_12 = D12/xxnorm
	if ( D1_12.gt.0.d0 ) then
		rot = ( dacos(D1_11) - (1.d0 + xnu)*xT )/twopi
		nsign = +1
	else
		rot = ( twopi - dacos(D1_11) - (1.d0 + xnu)*xT )/twopi
		nsign = -1
	end if
	rotate = rot - DINT(rot)
	if ( rotate.lt.0.d0 ) rotate = rotate + 1.d0
	FB(15) = xR - 2.d0*rotate
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