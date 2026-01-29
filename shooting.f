      PROGRAM SHOOTING
C     ------- --------
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PAR/xnu,xrho,xm,xgamma,xsigma,xepsilon,xlambda,xT,xdelta
      COMMON /CONST/epsilon,tolx,tolf,pi,nmax,NNR,NR0 
      DIMENSION vector1(10), vector2(10)
C
	xnu      = 1.d0/3.d0
	xrho     = 0.d0/2.d0 
	xm       = 17.d0/10.d0
	xgamma   = 0.d0/100.d0
	xsigma   = 0.d0/100.d0
	xepsilon = 0.d0/100.d0
	xlambda  = 0.d0/100.d0
C
	xT       = 0.d0
	xdelta   = 0.d0
C
	epsilon = 10.d0**(-5.d0)
	tolx    = 10.d0**(-12.d0)
	tolf    = 10.d0**(-12.d0)
	pi      = 4.d0*datan(1.d0)
	nmax    = 15
C
	R       = 1
C
	write(6,*), " "
	write(6, '(" nu       = ", 1f10.8)') xnu
	write(6, '(" xrho     = ", 1f10.8)') xrho
	write(6, '(" xm       = ", 1f10.8)') xm
	write(6, '(" xgamma   = ", 1f10.8)') xgamma
	write(6, '(" xsigma   = ", 1f10.8)') xsigma
	write(6, '(" xepsilon = ", 1f10.8)') xepsilon
	write(6, '(" xlambda  = ", 1f10.8)') xlambda
	write(6,*), " "
C
	call monodromy(vector1,vector2,idid)
C
	if (idid.eq.1) then
		print*, "Enter initial delta, T"
		read*, xdelta, xT
		print*, " "
		call newton(R,vector1,vector2,ierr)
	else
	end if
C
	if (ierr.eq.1) call postprocess(vector1,vector2)
C
      RETURN
      END