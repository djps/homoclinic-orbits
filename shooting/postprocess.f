      SUBROUTINE POSTPROCESS(vector1,vector2)
C     ---------- -----------
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0
      PARAMETER (n=13, nrdens=13, lrpar=1, lipar=1)                           
      PARAMETER (lwork=11*n+8*nrdens+20, liwork=nrdens+20)
      DIMENSION y(n),work(lwork),iwork(liwork),w(n)
      DIMENSION rpar(lrpar),ipar(lipar)
      DIMENSION vector1(10),vector2(10)
      EXTERNAL USRFUN_P, SOLOUT_P
C
	n_monodromy = 10
C
	x      = 0.d0
	xend   = 1.d0
C
	do i = 1,n_monodromy
		y(i) = epsilon*( dcos(xdelta)*vector1(i) + dsin(xdelta)*vector2(i) )
	end do
	do i = n_monodromy+1,n
		y(i) = 0.d0
	end do
	y(3)  = 1.d0 + y(3)
	y(6)  = 1.d0 + y(6)
	y(10) = 1.d0 + y(10)
C
	itol   = 0
C
	iout   = 2
C
	do i = 1,lwork
		work(i) = 0.d0
	end do
C
	do i = 1,liwork
		iwork(i) = 0
	end do
	iwork(3) =-1
	iwork(4) = 1
	iwork(5) = nrdens
C
	idid   =  0
C
	call dop853(n,usrfun_p,x,y,xend,tolx,tolf,itol,solout_p,iout,work,lwork,iwork,liwork,rpar,ipar,idid)
	print*, y
C       
	print*, " "
	if (idid.eq.+1) then 
		print*, "Computation successful"
		print*, " Number of function evaluations =", iwork(17)
		print*, " Number of computed steps       =", iwork(18)
		print*, " Number of accepted steps       =", iwork(19)
		print*, " Number of rejected steps       =", iwork(20)
	else
	end if
	if (idid.eq.+2) print*, "Computation successful (interupted by solout)"
	if (idid.eq.-1) print*, "Input is not consistant"
	if (idid.eq.-2) print*, "Larger nmax is needed"
	if (idid.eq.-3) print*, "Stepsize becomes too small"
	if (idid.eq.-4) print*, "Problem is probably stiff (interupted)"
	print*, " "        
C
	xD = xT - y(13)
C
	D11  = y(7)**2.d0 - y(8)**2.d0 - y(9)**2.d0 + y(10)**2.d0
	D12  = 2.d0*( y(7)*y(8) + y(9)*y(10) )
	if ( D11.lt.-1.d0 ) D11 = -1.d0
	if ( D11.gt.+1.d0 ) D11 = +1.d0
	xxnorm = DSQRT( D11**2.d0 + D12**2.d0)
	D1_11 = D11/xxnorm
	D1_12 = D12/xxnorm
	if ( D1_12.gt.0.d0 ) then
		rot = ( dacos(D1_11) - (1.d0 + xnu)*T )/(2.d0*pi)
		nsign = +1
	else
		rot = ( 2.d0*pi - dacos(D1_11) - (1.d0 + xnu)*T )/(2.d0*pi)
		nsign = -1
	end if
	xR = rot - dint(rot)
	if ( xR.lt.0.d0 ) xR = xR + 1.d0
C
	if (mod(xdelta,2.d0*pi).lt.0.d0) then
		xdelta = mod(xdelta,2.d0*pi) + 2.d0*pi
	else
		xdelta = mod(xdelta,2.d0*pi)
	end if
C
	write(6,'(" delta = ", 1e19.12)') xdelta
	write(6,'(" T     = ", 1e19.12)') xT
	write(6,'(" R     = ", 1e19.12)') xR
	write(6,'(" D     = ", 1e19.12)') xD
	print*, " "
C
	open(unit=2,file='/home/dave/auto/97/code/Magnetic/parameters.dat')
	write(2,'(8e20.12)') xnu,xrho,xm,xlambda,xT,xdelta,xD,xR
	write(2,*) " "
	write(2,'("  xnu      = ", 1e19.12)') xnu
	write(2,'("  xrho     = ", 1e19.12)') xrho
	write(2,'("  xm       = ", 1e19.12)') xm
	write(2,'("  xgamma   = ", 1e19.12)') xgamma
	write(2,'("  xsigma   = ", 1e19.12)') xsigma
	write(2,'("  xepsilon = ", 1e19.12)') xepsilon
	write(2,'("  xlambda  = ", 1e19.12)') xlambda
	write(2,'("  xT       = ", 1e19.12)') xT
	write(2,'("  xR       = ", 1e19.12)') xD
	write(2,'("  xD       = ", 1e19.12)') xR
	write(2,'("  epsilon  = ", 1e19.12)') epsilon
	write(2,'("  delta    = ", 1e19.12)') xdelta 
	close(2)
C
      RETURN
      END           
C
C
C
      SUBROUTINE USRFUN_P(N,X,Y,F,RPAR,IPAR)      
C     ---------- ------      
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0
      DIMENSION Y(N),F(N),RPAR(IPAR)
C
	F1  = Y(1)
	F2  = Y(2)
	F3  = Y(3)
	XM1 = Y(4)
	XM2 = Y(5)
	XM3 = Y(6)
	Q1  = Y(7)
	Q2  = Y(8)
	Q3  = Y(9)
	Q4  = Y(10)
	xX  = Y(11)
	xY  = Y(12)
	xZ  = Y(13)
C       
	F(1)  = (1.d0+xnu)*F2*XM3 - F3*XM2 + xlambda*2.d0*(Q1*Q4 + Q2*Q3)+
     + xlambda*xepsilon*F3*(1.d0+xsigma)*2.d0*(Q1*Q4+Q2*Q3)-
     + xlambda*xepsilon*F2*(Q3**2.d0+Q4**2.d0-Q1**2.d0-Q2**2.d0)
C
	F(2)  = (1.d0+xrho)*F3*XM1 - (1.d0+xnu)*F1*XM3 - xlambda*2.*(Q1*Q3 - Q2*Q4)+
     + xlambda*xepsilon*F1*(1.+xgamma)*2.*(Q3**2.+Q4**2.-Q1**2.-Q2**2.)-
     + xlambda*xepsilon*F3*(1.d0+xsigma)*2.d0*(Q1*Q3-Q2*Q4)
C
	F(3)  = F1*XM2 - (1.d0+xrho)*F2*XM1 +
     + xlambda*xepsilon*2.d0*F2*(Q1*Q3 - Q2*Q4) -
     + xlambda*xepsilon*(1.d0+xgamma)*F1*2.d0*(Q1*Q4 + Q2*Q3)
C
	F(4)  = xnu*XM2*XM3 + F2*(1.d0+xepsilon*xgamma*F3)/(xm**2.d0) 
	F(5)  = (xrho-xnu)*XM1*XM3 - F1*(1.d0+xepsilon*F3*(xgamma-xsigma))/(xm**2.d0)  
	F(6)  =-xrho*XM1*XM2 - xepsilon*xsigma*F1*F2/xm**2.d0 
C
	F(7)  = (1.d0/2.d0)*( +(1.d0+xrho)*XM1*Q4 - XM2*Q3 + (1.d0+xnu)*XM3*Q2 )
	F(8)  = (1.d0/2.d0)*( +(1.d0+xrho)*XM1*Q3 + XM2*Q4 - (1.d0+xnu)*XM3*Q1 )
	F(9)  = (1.d0/2.d0)*( -(1.d0+xrho)*XM1*Q2 + XM2*Q1 + (1.d0+xnu)*XM3*Q4 )
	F(10) = (1.d0/2.d0)*( -(1.d0+xrho)*XM1*Q1 - XM2*Q2 - (1.d0+xnu)*XM3*Q3 )
C
	F(11) = (1.d0+xepsilon*F3*(1.d0+xgamma))*2.d0*(Q1*Q3 + Q2*Q4)+
     + xepsilon*F2*2.d0*(Q1*Q2-Q3*Q4)+
     + xepsilon*F1*(1.d0+xsigma)*(Q1**2.d0-Q2**2.d0-Q3**2.d0+Q4**2.d0)
C
	F(12) = (1.d0+xepsilon*F3*(1.d0+xgamma))*2.d0*(Q2*Q3 - Q1*Q4)+
     + xepsilon*F2*(Q2**2.d0+Q4**2.d0-Q1**2.d0-Q3**2.d0)+
     + xepsilon*F1*(1.d0+xsigma)*2.d0*(Q1*Q2+Q3*Q4)
C
	F(13) = (1.d0+xepsilon*F3*(1.d0+xgamma))*(Q3**2.d0+Q4**2.d0-Q1**2.d0-Q2**2.d0)+
     + xepsilon*F2*2.d0*(Q1*Q4+Q2*Q3)+
     + xepsilon*F1*(1.d0+xsigma)*2.d0*(Q1*Q3-Q2*Q4)
C
	do i=1,n
		F(i) = xT*F(i)
	end do
C
      RETURN
      END      
C
C
C
      SUBROUTINE SOLOUT_P(NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)      
C     ---------- ------      
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /intern/xout      
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0
      PARAMETER(np=13, xplot=5d-4, lrpar=1, lipar=1)
      DIMENSION y(n),z(np),con(8*nd),icomp(nd),w(np)
      DIMENSION rpar(lrpar), ipar(lipar)
C
	if (nr.eq.1) then
		open(unit=2,file='/home/dave/auto/97/code/Magnetic/mag.dat')
		write(2,100) x, (y(i),i=1,n)
		xout=x+xplot
	else
10		continue
		if (x.ge.xout) then
			do i=1,n
				z(i) = contd8(i,xout,con,icomp,nd)
			end do
			open(unit=2,file='/home/dave/auto/97/code/Magnetic/mag.dat')
			write(2,100) xout, (z(i),i=1,n)
			xout = xout+xplot
			goto 10
		end if
	end if
C
100	format(14e20.10)
C
      RETURN
      END      