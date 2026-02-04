      SUBROUTINE MONODROMY(vector1,vector2,idid)
C     ------- --------
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0 
      EXTERNAL USRFUN_M, SOLOUT_M
      PARAMETER (n=6, nrdens=6, lwork=11*n+8*nrdens+20)
      PARAMETER (liwork=nrdens+20, lrpar=1, lipar=1)
	  PARAMETER (m=10)
      DIMENSION z(n,n),vector1(m),vector2(m)
      DIMENSION vr(n,n),vi(n,n),vrdum(n,n),vidum(n,n)
      DIMENSION xMonodromy(n,n),xnorm(n),wrdum(n),widum(n)
      DIMENSION fv1(n),y(n)
      DIMENSION wr(n),wi(n),xInitial(n,n),iv1(n)
      DIMENSION rpar(lrpar),ipar(lipar),work(lwork),iwork(liwork)

      DIMENSION xMonodromy_copy(n,n), work_ev(100)
      INTEGER lwork_ev, info_ev
      CHARACTER jobvl, jobvr

C
	x    = 0.d0
	xend = 4.d0 * pi / (1.d0 + xnu)
C
	itol   = 0
C
	iout   = 0
C
	do i = 1,lwork
		work(i) = 0.d0
	end do
C
	do i = 1,liwork
		iwork(i) = 0
	end do
C
	irtrn = 0
C
	iwork(3) =-1
	iwork(4) =+1
	iwork(5) = nrdens
C
	idid     =  0
C
	matz = 1
	zero = 0.00000001
	do i = 1,n
		do j = 1,n
			vr(i,j)         = 0.d0
			vi(i,j)         = 0.d0
			vrdum(i,j)      = 0.d0
			vidum(i,j)      = 0.d0
			xInitial(i,j)   = 0.d0
			xMonodromy(i,j) = 0.d0
		end do
		xInitial(i,i) = 1.d0
		xnorm(i)      = 0.d0
		wrdum(i)      = 0.d0
		widum(i)      = 0.d0
	end do
C
C Create Monodromy matrix
C		
	do k = 1,n
		do i = 1,n
			y(i) = xInitial(i,k)
		end do
		call dop853(n,usrfun_m,x,y,xend,tolx,tolf,itol,solout_m,iout,work,lwork,iwork,liwork,rpar,ipar,idid)
		x    = 0.d0
		xend = 4.d0*pi/(1.d0+xnu)
		do i = 1,n
			xMonodromy(i,k) = y(i)
		end do
	end do
C
	if (idid.ne.1) then
		print*, "argh", idid
		return
	else
	end if
C
! 	print*, "Monodromy Matrix"
! 	do k = 1,n 
! 		write(6,100) (xMonodromy(k,i), i=1,n)
! 	end do

C
C Find floquet multipliers using LAPACK DGEEV
C        
C --- Copy matrix (DGEEV destroys input)
	do k = 1,n
		do i = 1,n
		xMonodromy_copy(i,k) = xMonodromy(i,k)
		end do
	end do
	
C --- Setup DGEEV parameters
	jobvl = 'N'
	jobvr = 'V'
	lwork_ev = 100
	
C --- Call LAPACK DGEEV instead of EISPACK RG
	call DGEEV(jobvl, jobvr, n, xMonodromy_copy, n, wr, wi, z, n, z, n, work_ev, lwork_ev, info_ev)
	
	ierr = info_ev
	if (ierr.ne.0) then
		write(6,*) " "
		write(6,*) "DGEEV failed with INFO = ", ierr
		write(6,*) " "
		return
	else
		idid = 1
	end if


!	call rg(n,n,xMonodromy,wr,wi,matz,z,iv1,fv1,ierr)
!	if (ierr.ne.0) then
!		print*, "argh"
!		return
!	else
!		idid = 1
!	end if
C
!	if (ierr.ne.0) then
!		write(6,*) " "
!		write(6,*) "failure in eigenvalue, eigenvector solver"
!!	else
! 		write(6,*) " "
! 		write(6,*) "real part of eigenvalues"
! 		do i = 1,n 
! 			write(6,'(1e18.10)') wr(i)
! 		end do
! 		write(6,*) " "
! 		write(6,*) "imaginary part of eigenvalues"
! 		do i = 1,n 
! 			write(6,'(1e18.10)') wi(i)
! 		end do	
! 		write(6,*) " "
! 		write(6,*) "magnitude of eigenvalues"
! 		do i = 1,n 
! 			write(6,'(1e18.10)') dsqrt(wi(i)**2.d0 + wr(i)**2.d0)
! 		end do	
! 		write(6,*) " "
!	end if


C
C Create matrices vr, vi from output matrix z
C
	iflag = 0
	do i = 1,n
		if ( wi(i+iflag).eq.0.00 ) then
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
! 	write(6,*) "normalised REAL matrix"
! 	do i = 1,n
! 		write(6,100) ( vr(i,j), j = 1,n )
! 	end do
! 	write(6,*) " "
! 	write(6,*) "normalised IMAGINARY matrix"
! 	do i = 1,n
! 		write(6,100) ( vi(i,j), j = 1,n )
! 	end do
! 	write(6,*) " "
C
C Order vectors according to magnitude of eigenvalues
C
	do i = 1,n-1
		do j = i+1,n
			if ( (wi(i)**2.d0 + wr(i)**2.d0).gt.(wi(j)**2.d0 + wr(j)**2.d0) ) then 
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
! 	write(6,*) "final REAL matrix"
! 	do i = 1,n
! 		write(6,100) ( vr(i,j), j = 1,n )
! 	end do
! 	write(6,*) " "
! 	write(6,*) "final IMAGINARY matrix"
! 	do i = 1,n
! 		write(6,100) ( vi(i,j), j = 1,n )
! 	end do
! 	write(6,*) " "
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
! 	write(6,*) "final matrix"
! 	do i = 1,n
! 		write(6,100) ( vr(i,j), j = 1,n )
! 	end do
! 	write(6,*) " "
! C
! 	write(6,*) "required vectors"
! 	do i = n-1,n
! 		write(6,100) ( vr(j,i), j = 1,n  )
! 	end do
! 	write(6,*) " "
C
	do i =1,m
		vector1(i) = 0.d0
		vector2(i) = 0.d0
	end do
	do i = 1,2
		vector1(i)   = vr(i+0,n-1)
		vector2(i)   = vr(i+0,n)
		vector1(i+3) = vr(i+2,n-1)
		vector2(i+3) = vr(i+2,n)
		vector1(i+6) = vr(i+4,n-1)
		vector2(i+6) = vr(i+4,n)
	end do
C
	p = abs( (pi*xend)/( 2.d0*( datan( wi(n)/wr(n) ) - 2.d0*pi ) ) )

C
	write(6,'(" imaginary part of floquet multipler = ", 1e19.10)') wi(n)
	write(6,'(" real part of floquet multipler      = ", 1e19.10)') wr(n)
	write(6,'(" absolute value of floquet multipler = ", 1e19.10)') dsqrt(wi(n)**2.d0 + wr(n)**2.d0)
	write(6,'(" argument of floquet multipler       = ", 1e19.10)') datan(wi(n)/wr(n))
	write(6,'(" real part of eigenvalue             = ", 1e19.10)') (1.d0/(xend*2.d0))*dlog( wi(n)**2.d0 + wr(n)**2.d0 )
	write(6,'(" imaginary part of eigenvalue        = ", 1e19.10)') ( datan(wi(n)/wr(n)) - 2.d0*pi )/xend
	write(6,'(" period of multimodals               = ", 1e19.10)') pi - nint(p/pi)*p
	write(6,*) " "
C
	write(6,*) "required vectors"
	write(6,'(10e12.4)') ( vector1(i), i = 1,m )
	write(6,'(10e12.4)') ( vector2(i), i = 1,m )
	write(6,*) " "
C
100	format(6e12.4)
C
      RETURN
      END
C
C
C
      SUBROUTINE USRFUN_M(N,X,Y,F,RPAR,IPAR)
C     ---------- ------ 
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0
      PARAMETER (lrpar=1, lipar=1) 
      DIMENSION Y(N), F(N), rpar(lrpar),ipar(lipar)
C
	F1  = Y(1)
	F2  = Y(2)
	XM1 = Y(3)
	XM2 = Y(4)
	Q1  = Y(5)
	Q2  = Y(6)
C
	Q3_star = dsin( X*(1.d0+xnu)/2.d0 )
	Q4_star = dcos( X*(1.d0+xnu)/2.d0 )
C
	F(1) = (1.d0+xnu)*F2 - XM2 + 2.d0*xlambda*(Q4_star*Q1+Q3_star*Q2) +
     + xepsilon*xlambda*(1.d0+xgamma)*2.*(Q4_star*Q1+Q3_star*Q2) -
     + xepsilon*xlambda*F2 
C
	F(2) =-(1.d0+xnu)*F1 + (1.d0+xrho)*XM1 - 2.d0*xlambda*(Q3_star*Q1-Q4_star*Q2) +
     + xepsilon*xlambda*((1.d0+xsigma)*F1) -
     + xepsilon*xlambda*(1.d0+xgamma)*2.*(Q3_star*Q1-Q4_star*Q2)
C
	F(3) = F2*(1.d0+xepsilon*xsigma)/(xm**2.d0) + xnu*XM2
C
	F(4) =-F1*(1.d0+xepsilon*(xsigma-xgamma))/(xm**2.d0) + (xrho-xnu)*XM1
C
	F(5) = (1.d0+xrho)*Q4_star*XM1/2.d0 + (1.d0+xnu)*Q2/2.d0 - Q3_star*XM2/2.d0
C
	F(6) = (1.d0+xrho)*Q3_star*XM1/2.d0 - (1.d0+xnu)*Q1/2.d0 + Q4_star*XM2/2.d0 
C
      RETURN
      END      
C
C
C
      SUBROUTINE SOLOUT_M(NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN) 
C     ---------- ------ 
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /intern/xout      
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0
      PARAMETER(np=10, xplot=2d-3, lrpar=1, lipar=1)
      DIMENSION y(n),z(np),con(8*nd),icomp(nd)
      DIMENSION rpar(lrpar), ipar(lipar)
C
      RETURN
      END
