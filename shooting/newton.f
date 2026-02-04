      SUBROUTINE NEWTON(R,vector1,vector2,idid)
C     ---------- ------
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0    
      PARAMETER (n=30, nrdens=30, lrpar=1, lipar=1)                           
      PARAMETER (lwork=11*n+8*nrdens+20, liwork=nrdens+20)      
      DIMENSION s(2),G(2),DG(2,2),DG_INV(2,2)      
      DIMENSION T_residue(nmax+1),delta_residue(nmax+1) 
      DIMENSION vector1(10),vector2(10)
      DIMENSION y(n),work(lwork),iwork(liwork),rpar(lrpar),ipar(lipar)
      EXTERNAL USRFUN_N, SOLOUT_N  
C
	x      = 0.d0
	xend   = 1.d0
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
C
	irtrn = 0
C
	iwork(3) =-1
	iwork(4) =+1
	iwork(5) = nrdens
C
	idid     =  0
C
	s(1) = xdelta
	s(2) = xT
C
	n_monodromy  = 10
	n_iterations = 1
C
	delta_residue(1) = 1.d0
	T_residue(1)     = 1.d0
C
	if (R.eq.1) then
		n_1 = 1
		n_2 = 4
	else 
		n_1 = 2
		n_2 = 5
	end if
C
	print*, "n    delta              T                  residual delta     residual T"
	write(6,'(1i2,4f20.14)') n_iterations-1, s, dabs(delta_residue(n_iterations)), dabs(T_residue(n_iterations)) 
C
	call cpu_time(t1)
C
100	if ( n_iterations.le.nmax ) then
C  
200		if ( (dabs(delta_residue(n_iterations))+dabs(T_residue(n_iterations)).ge.tolx).or.dabs(y(n_1))+dabs(y(n_2)).ge.tolf) then
C
			do i = 1,n_monodromy
				y(i)               = epsilon*( dcos(xdelta)*vector1(i) + dsin(xdelta)*vector2(i) )
				y(i+1*n_monodromy) = epsilon*( dcos(xdelta)*vector2(i) - dsin(xdelta)*vector1(i) )
				y(i+2*n_monodromy) = 0.d0
			end do
			y(10) = 1.d0 + y(10)
			y(3)  = 1.d0 + y(3)
			y(6)  = 1.d0 + y(6)
C
			do i = 1,lwork
				work(i) = 0.d0
			end do
			do i = 1,liwork
				iwork(i) = 0
			end do
			iwork(3) =-1
			iwork(4) =+1
			iwork(5) = nrdens 
C
			call dop853(n,usrfun_n,x,y,xend,tolx,tolf,itol,solout_n,iout,work,lwork,iwork,liwork,rpar,ipar,idid)
C
			G(1)    = y(n_1)
			G(2)    = y(n_2)
C
			DG(1,1) = y(n_1+1*n_monodromy)
			DG(1,2) = y(n_1+2*n_monodromy)
			DG(2,1) = y(n_2+1*n_monodromy)
			DG(2,2) = y(n_2+2*n_monodromy)
C
			det = DG(1,1)*DG(2,2) - DG(1,2)*DG(2,1)
C
			DG_INV(1,1) = +DG(2,2)/det 
			DG_INV(1,2) = -DG(1,2)/det
			DG_INV(2,1) = -DG(2,1)/det
			DG_INV(2,2) = +DG(1,1)/det
C
			s(1) = s(1) - ( DG_INV(1,1)*G(1) + DG_INV(1,2)*G(2) )
			s(2) = s(2) - ( DG_INV(2,1)*G(1) + DG_INV(2,2)*G(2) )
C
			n_iterations                = n_iterations + 1
			delta_residue(n_iterations) = s(1) - xdelta
			T_residue(n_iterations)     = s(2) - xT
			x = 0.d0 
C
			write(6,'(1i2,4f20.14)') n_iterations-1, s, dabs(delta_residue(n_iterations)), dabs(T_residue(n_iterations)) 
C
			xdelta = s(1)
			xT     = s(2)
C
			goto 100
		else 
			print*, " "
			idid = 1
			print*, "program complete - tolerances reached"
			write(6,'("  delta tolerance = ", 1e19.12)') dabs(T_residue(n_iterations))
			write(6,'("  T tolerance     = ", 1e19.12)') dabs(delta_residue(n_iterations))
			if (R.eq.1) then
				write(6,'("  F1(T)           = ", 1e19.12)') y(n_1)
				write(6,'("  M1(T)           = ", 1e19.12)') y(n_2)
			else
				write(6,'("  F2(T)           = ", 1e19.12)') y(n_1)
				write(6,'("  M2(T)           = ", 1e19.12)') y(n_2)
			end if
			call cpu_time(t2)
			print*, " "
			print*, "time taken=", t2-t1, " seconds"
			goto 400
		end if
		goto 100
	else
		print*, " "
		idid = 2  
		print*, "program finished - maximum iterations reached:", n_iterations-1
		write(6,'("  delta tolerance = ", 1e19.12)') dabs(T_residue(n_iterations))
		write(6,'("  T tolerance     = ", 1e19.12)') dabs(delta_residue(n_iterations))
		if (R.eq.1) then
			write(6,'("  F1(T)           = ", 1e19.12)') y(n_1)
			write(6,'("  M1(T)           = ", 1e19.12)') y(n_2)
		else
			write(6,'("  F2(T)           = ", 1e19.12)') y(n_1)
			write(6,'("  M2(T)           = ", 1e19.12)') y(n_2)
		end if
		call cpu_time(t2)
		print*, " "
		print*, "time taken=", t2-t1, " seconds"
		goto 400
400	end if
C
      RETURN
      END
C
C
C
      SUBROUTINE USRFUN_N(N,X,Y,F,RPAR,IPAR)      
C     ---------- ------      
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0      
      PARAMETER (lrpar=1, lipar=1) 
      DIMENSION Y(N), F(N), RPAR(LRPAR), IPAR(LIPAR)      
C
	x1  = Y(1)
	x2  = Y(2)
	x3  = Y(3)
	x4  = Y(4)
	x5  = Y(5)
	x6  = Y(6)
	q1  = Y(7)
	q2  = Y(8)
	q3  = Y(9)
	q4  = Y(10)
C
	w1  = Y(11)
	w2  = Y(12)
	w3  = Y(13)
	w4  = Y(14)
	w5  = Y(15)
	w6  = Y(16)
	qw1 = Y(17)
	qw2 = Y(18)
	qw3 = Y(19)
	qw4 = Y(20)
C
	v1  = Y(21)
	v2  = Y(22)
	v3  = Y(23)
	v4  = Y(24)
	v5  = Y(25)
	v6  = Y(26)
	qv1 = Y(27)
	qv2 = Y(28)
	qv3 = Y(29)
	qv4 = Y(30)  
C
	xep = xepsilon
	xl  = xlambda
C
	F(1)  = xT*( (1.d0+xnu)*x2*x6 - x3*x5 + xlambda*2.d0*(q1*q4 + q2*q3) ) +  
     + xT*( xlambda*xepsilon*x3*(1.d0+xsigma)*2.d0*(q1*q4+q2*q3) ) -
     + xT*( xlambda*xepsilon*x2*(q3**2.d0+q4**2.d0-q1**2.d0-q2**2.d0) )
C
	F(2)  = xT*( (1.d0+xrho)*x3*x4-(1.d0+xnu)*x1*x6-xlambda*2.d0*(q1*q3-q2*q4) ) +
     + xT*(xl*xep*x1*(1.+xgamma)*(q3**2.+q4**2.-q1**2.-q2**2.)) -
     + xT*(xl*xep*x3*(1.+xsigma)*2.d0*(q1*q3-q2*q4))   
C
	F(3)  = xT*( x1*x5 - (1.d0+xrho)*x2*x4 ) +    
     + xT*( xlambda*xepsilon*2.d0*x2*(q1*q3 - q2*q4) ) -
     + xT*( xlambda*xepsilon*(1.d0+xgamma)*x1*2.d0*(q1*q4 + q2*q3) )
C
	F(4)  = xT*( xnu*x5*x6 + x2*(1.d0 + xepsilon*xgamma*x3)/(xm**2.d0) )
	F(5)  = xT*( (xrho-xnu)*x4*x6 - x1*(1.d0 + xepsilon*(xgamma-xsigma)*x3)/(xm**2.d0) )
	F(6)  = xT*( -xrho*x4*x5 - xepsilon*xsigma*x1*x2/xm**2.d0 )
	F(7)  = (xT/2.d0)*( +(xrho+1.d0)*x4*q4 - x5*q3 + x6*(xnu+1.d0)*q2 )
	F(8)  = (xT/2.d0)*( +(xrho+1.d0)*x4*q3 + x5*q4 - x6*(xnu+1.d0)*q1 )
	F(9)  = (xT/2.d0)*( -(xrho+1.d0)*x4*q2 + x5*q1 + x6*(xnu+1.d0)*q4 )
	F(10) = (xT/2.d0)*( -(xrho+1.d0)*x4*q1 - x5*q2 - x6*(xnu+1.d0)*q3 )  
C    
	F(11) = xT*( (1.d0+xnu)*( w2*x6 + x2*w6 ) - ( w5*x3 + x5*w3 )) + 
     + xT*(xlambda*2.d0*(qw1*q4 + q1*qw4 + qw2*q3 + q2*qw3)) +
     + xT*( xl*xep*w3*(1.d0+xsigma)*2.d0*(q1*q4+q2*q3) ) -
     + xT*( xl*xep*w2*(q3**2.d0+q4**2.d0-q1**2.d0-q2**2.d0) ) +
     + xT*( xl*xep*x3*(1.d0+xsigma)*2.d0*(q1*qw4+q2*qw3) ) -
     + xT*( xl*xep*x2*2.*(q3*qw3+q4*qw4-q1*qw1-q2*qw2) ) -
     + xT*( xl*xep*x3*(1.d0+xsigma)*2.d0*(qw1*q4+qw2*q3) )
C
	F(12) = xT*( (1.d0+xrho)*( w4*x3 + x4*w3 ) - (1.d0+xnu)*( w1*x6 + x1*w6 ) ) - 
     + xT*(xlambda*2.d0*(qw1*q3 + q1*qw3 - qw2*q4 - q2*qw4)) +
     + xT*(xl*xep*w1*(1.+xgamma)*(q3**2.+q4**2.-q1**2.-q2**2.) ) -
     + xT*(xl*xep*w3*(1.+xsigma)*2.*(q1*q3-q2*q4) ) +
     + xT*(xl*xep*x1*(1.+xgamma)*2.*(qw3*q3+qw4*q4-qw1*q1-qw2*q2) ) -
     + xT*(xl*xep*x3*(1.+xsigma)*2.*(qw1*q3-qw2*q4) ) -
     + xT*(xl*xep*x3*(1.+xsigma)*2.*(q1*qw3-q2*qw4) )
C
	F(13) = xT*( w1*x5 + x1*w5 - (1.d0+xrho)*(w2*x4 + x2*w4) ) +
     + xT*(xl*xep*2.*(w2*(q1*q3-q2*q4)-x2*(qw1*q3-qw2*q4)) )-     
     + xT*(xl*xep*2.*x2*(q1*qw3-q2*qw4)) +
     + xT*(xl*xep*(1.+xgamma)*2.*(w1*(q1*q4+q2*q3)+x1*(qw1*q4+qw2*q3)))+
     + xT*(xl*xep*(1.+xgamma)*2.*x1*(q1*qw4+q2*qw3))
C
	F(14) = xT*( xnu*( w5*x6 + x5*w6 ) + (w2*(1.d0+xepsilon*xgamma*x3)+x2*xepsilon*xgamma*w3 )/(xm**2.d0) )
	F(15) = xT*((xrho-xnu)*(w4*x6+x4*w6) - (w1*(1.+xepsilon*(xgamma-xsigma)*x3)+x1*xepsilon*(xgamma-xsigma)*w3)/(xm**2.d0) )
	F(16) = xT*( -xrho*(w4*x5 + x4*w5) - xepsilon*xsigma*(w1*x2+x1*w2)/xm**2.d0 )
	F(17) = (xT/2.d0)*( +(xrho+1.d0)*(w4*q4+x4*qw4) - (w5*q3+x5*qw3) + (xnu+1.d0)*(qw2*x6+q2*w6))
	F(18) = (xT/2.d0)*( +(xrho+1.d0)*(w4*q3+x4*qw3) + (w5*q4+x5*qw4) - (xnu+1.d0)*(qw1*x6+q1*w6))
	F(19) = (xT/2.d0)*( -(xrho+1.d0)*(w4*q2+x4*qw2) + (w5*q1+x5*qw1) + (xnu+1.d0)*(qw4*x6+q4*w6))
	F(20) = (xT/2.d0)*( -(xrho+1.d0)*(w4*q1+x4*qw1) - (w5*q2+x5*qw2) - (xnu+1.d0)*(qw3*x6+q3*w6)) 
C    
	F(21) = xT*( (1.d0+xnu)*( v2*x6 + x2*v6) - ( v5*x3 + x5*v3 ) ) + 
     + xT*(xlambda*2.d0*(qv1*q4 + q1*qv4 + qv2*q3 + q2*qv3)) +
     + (1.d0+xnu)*x2*x6-x3*x5+xlambda*2.d0*(q1*q4+q2*q3) +
     + xl*xep*x3*(1.d0+xsigma)*2.d0*(q1*q4+q2*q3) -
     + xl*xep*x2*(q3**2.d0+q4**2.d0-q1**2.d0-q2**2.d0) +
     + xT*( xl*xep*v3*(1.d0+xsigma)*2.d0*(q1*q4+q2*q3) ) -
     + xT*( xl*xep*v2*(q3**2.d0+q4**2.d0-q1**2.d0-q2**2.d0) ) +
     + xT*( xl*xep*x3*(1.d0+xsigma)*2.d0*(q1*qv4+q2*qv3) ) -
     + xT*( xl*xep*x2*2.*(q3*qv3+q4*qv4-q1*qv1-q2*qv2) ) -
     + xT*( xl*xep*x3*(1.d0+xsigma)*2.d0*(qv1*q4+qv2*q3) )
C
	F(22) = xT*( (1.d0+xrho)*(v4*x3+x4*v3)-(1.d0+xnu)*(v1*x6+x1*v6) ) - 
     + xT*(xlambda*2.d0*(qv1*q3 + q1*qv3 - qv2*q4 - q2*qv4)) + 
     + (1.+xrho)*x3*x4-(1.+xnu)*x1*x6-xlambda*2.*(q1*q3-q2*q4) + 
     + xlambda*xep*x1*(1.+xgamma)*2.*(q3**2.+q4**2.-q1**2.-q2**2.) -
     + xlambda*xep*x3*(1.+xsigma)*2.*(q1*q3-q2*q4) +
     + xT*(xl*xep*v1*(1.+xgamma)*2.*(q3**2.+q4**2.-q1**2.-q2**2.) ) -
     + xT*(xl*xep*v3*(1.+xsigma)*2.*(q1*q3-q2*q4) ) +
     + xT*(xl*xep*x1*(1.+xgamma)*4.*(qv3*q3+qv4*q4-qv1*q1-qv2*q2) ) -
     + xT*(xl*xep*x3*(1.+xsigma)*2.*(qv1*q3-qv2*q4) ) -
     + xT*(xl*xep*x3*(1.+xsigma)*2.*(q1*qv3-q2*qv4) )
C
	F(23) = xT*( v1*x5 + x1*v5 - (1.d0+xrho)*( v2*x4 + x2*v4 ) ) + 
     + xT*(xl*xep*2.*(v2*(q1*q3-q2*q4)-x2*(qv1*q3-qv2*q4))) -     
     + xT*(xl*xep*2.*x2*(q1*qv3-q2*qv4)) +
     + xT*(xl*xep*(1.+xgamma)*2.*(v1*(q1*q4+q2*q3)+x1*(qv1*q4+qv2*q3)))+
     + xT*(xl*xep*(1.+xgamma)*2.*x1*(q1*qv4+q2*qv3))+
     + x1*x5 - (1.d0+xrho)*x2*x4 +
     + xl*xep*2.*x2*(q1*q3 - q2*q4) -
     + xl*xep*(1.+xgamma)*x1*2.*(q1*q4 + q2*q3)
C
	F(24) = xT*( xnu*( v5*x6 + x5*v6 ) + (v2*(1.d0+xepsilon*xgamma*x3)+x2*xepsilon*xgamma*v3 )/(xm**2.d0) )+ 
     + xnu*x5*x6 + x2*(1.d0 + xepsilon*xgamma*x3)/(xm**2.d0)
C
	F(25) = xT*((xrho-xnu)*(v4*x6+x4*v6)-(v1*(1.+xepsilon*(xgamma-xsigma)*x3)+x1*xepsilon*(xgamma-xsigma)*v3)/(xm**2.d0) ) +
     + (xrho-xnu)*x4*x6-x1*(1.+xepsilon*(xgamma-xsigma)*x3)/(xm**2.d0)
C
	F(26) = xT*( -xrho*(v4*x5 + x4*v5) - xepsilon*xsigma*(v1*x2+x1*v2)/xm**2.d0 ) - 
     + xrho*x4*x5 - xepsilon*xsigma*x1*x2/xm**2.d0
C
	F(27) = (xT/2.d0)*( +(xrho+1.d0)*(v4*q4+x4*qv4) - (v5*q3+x5*qv3) + (xnu+1.d0)*(qv2*x6+q2*v6)) + 
     + (1.d0/2.d0)*(+(xrho+1.d0)*x4*q4-x5*q3+x6*(xnu+1.d0)*q2)
C
	F(28) = (xT/2.d0)*( +(xrho+1.d0)*(v4*q3+x4*qv3) + (v5*q4+x5*qv4) - (xnu+1.d0)*(qv1*x6+q1*v6)) +
     + (1.d0/2.d0)*(+(xrho+1.d0)*x4*q3+x5*q4-x6*(xnu+1.d0)*q1)
C
	F(29) = (xT/2.d0)*( -(xrho+1.d0)*(v4*q2+x4*qv2) + (v5*q1+x5*qv1) + (xnu+1.d0)*(qv4*x6+q4*v6)) +
     + (1.d0/2.d0)*(-(xrho+1.d0)*x4*q2+x5*q1+x6*(xnu+1.d0)*q4)
C
	F(30) = (xT/2.d0)*( -(xrho+1.d0)*(v4*q1+x4*qv1) - (v5*q2+x5*qv2) - (xnu+1.d0)*(qv3*x6+q3*v6)) +
     + (1.d0/2.d0)*(-(xrho+1.d0)*x4*q1-x5*q2-x6*(xnu+1.d0)*q3)    
C
      RETURN
      END      
C
C
C
      SUBROUTINE SOLOUT_N(NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)      
C     ---------- ------      
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /intern/xout      
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0
      PARAMETER(np=30, xplot=2d-3, lrpar=1, lipar=1)
      DIMENSION y(n),z(np),con(8*nd),icomp(nd)
      DIMENSION RPAR(LRPAR), IPAR(LIPAR)
C
	if (NR0.eq.1) then NR=1
C
	if (NR.eq.1) then
		write(9,100) x,(y(i),i=1,n)
		xout=x+xplot
	else
10		continue
		if (x.ge.xout) then
			do i=1,n
				z(i) = contd8(i,xout,con,icomp,nd)
			end do
			write(9,100) xout,(z(i),i=1,n)
			xout = xout+xplot
			goto 10
		end if
	end if
C
	NR0 = 2
C
100	format(31e14.4)
C
      RETURN
      END
