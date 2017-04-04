	subroutine Visibility(image,b,angle,lam0,V,phase,Diam)
	use Parameters
	IMPLICIT NONE
	type(RPhiImage) image
	real*8 b,lam0,V,k,phi,theta,phi1,phi2,angle,phase ! Gijsexp
	real*8 Im1,Im2,Diam,Rsigma,Im1_norm,Im2_norm
	complex*32 cI,Int1,Int2,r1,r2,cV,Ftot,a1,b1,c1,expa1r1,expc1r1,expa1r2,expc1r2
	integer i,j

	if(abs(angle).lt.(pi/360d0)) then
		theta=pi/2d0-pi/360d0+D%PA*pi/180d0
	else
		theta=pi/2d0-angle+D%PA*pi/180d0
	endif		
	cI=(0d0,1d0)
	k=2d0*pi*b/(lam0*D%distance/parsec)
	k=k*1d6*pi/(180d0*3600d0)
	
	write(*,'("Computing visibility")')
	write(*,'("Baseline:",f10.3," meter")') b
	write(*,'("Angle:   ",f10.3," degrees")') angle*180d0/pi
	write(9,'("Computing visibility")')
	write(9,'("Baseline:",f10.3," meter")') b
	write(9,'("Angle:   ",f10.3," degrees")') angle*180d0/pi
	cV=0d0
	Ftot=0d0

	if(Diam.gt.0d0) then
		Rsigma=0.42*lam0*1d-6/Diam
		Rsigma=Rsigma*180d0/pi
		Rsigma=Rsigma*3600d0
		Rsigma=Rsigma*D%distance/parsec
	else
		Rsigma=1d200
	endif

	do i=1,image%nr-1
	call tellertje(i,image%nr-1)
	do j=1,image%nphi
	r1=image%R(i)
	r2=image%R(i+1)

	Im1=image%image(i,j)*exp(-r1**2/(2d0*Rsigma**2))
	Im2=image%image(i+1,j)*exp(-r2**2/(2d0*Rsigma**2))

	if(vis_norm_fov) then
		Im1_norm=Im1
		Im2_norm=Im2
	else
		Im1_norm=image%image(i,j)
		Im2_norm=image%image(i+1,j)
	endif

	if(r1.eq.r2) goto 1

	phi1=image%phi(j)-0.5d0*pi/real(image%nphi)
	phi2=image%phi(j)+0.5d0*pi/real(image%nphi)

	c1=cos(theta+phi1)
	b1=sin((phi1+phi2)/2d0+theta)
	a1=c1+2d0*sin((phi1-phi2)/2d0)*b1

	if(abs(k*r2).gt.0.1d0.or.abs(k*r1).gt.0.1d0) then
	expa1r1=cqexp(-cI*k*r1*a1)
	expa1r2=cqexp(-cI*k*r2*a1)
	expc1r1=cqexp(-cI*k*r1*c1)
	expc1r2=cqexp(-cI*k*r2*c1)
	Int1=(-1d0/(k**2*b1)) * ((expa1r2-expa1r1)/a1-
     &						 (expc1r2-expc1r1)/c1)
	Int2=(-1d0/(cI*k**3*b1)) * 
     &	((1d0+cI*k*a1*r2)*expa1r2/a1**2-
     &	 (1d0+cI*k*c1*r2)*expc1r2/c1**2-
     &	 (1d0+cI*k*a1*r1)*expa1r1/a1**2+
     &	 (1d0+cI*k*c1*r1)*expc1r1/c1**2)
	else
	Int1=(r2**2-r1**2)*(a1-c1)/(2d0*b1)
     & -(cI*(r2**3-r1**3)*(a1**2-c1**2)/(6d0*b1))*k
     & -((r2**4-r1**4)*(a1**3-c1**3)/(24d0*b1))*k**2
     & +(cI*(r2**5-r1**5)*(a1**4-c1**4)/(120d0*b1))*k**3
     & +((r2**6-r1**6)*(a1**5-c1**5)/(720d0*b1))*k**4
     & -(cI*(r2**7-r1**7)*(a1**6-c1**6)/(5040d0*b1))*k**5
     & -((r2**8-r1**8)*(a1**7-c1**7)/(40320d0*b1))*k**6
	Int2=(r2**3-r1**3)*(a1-c1)/(3d0*b1)
     & -(cI*(r2**4-r1**4)*(a1**2-c1**2)/(8d0*b1))*k
     & -((r2**5-r1**5)*(a1**3-c1**3)/(30d0*b1))*k**2
     & +(cI*(r2**6-r1**6)*(a1**4-c1**4)/(144d0*b1))*k**3
     & +((r2**7-r1**7)*(a1**5-c1**5)/(840d0*b1))*k**4
     & -(cI*(r2**8-r1**8)*(a1**6-c1**6)/(5760d0*b1))*k**5
     & -((r2**9-r1**9)*(a1**7-c1**7)/(45360d0*b1))*k**6
	endif


	cV=cV+Im1*((1d0+r1/(r2-r1))*Int1-Int2/((r2-r1)))
	cV=cV+Im2*(Int2-r1*Int1)/((r2-r1))

	Int1=(r2**2-r1**2)*(a1-c1)/(2d0*b1)
	Int2=(r2**3-r1**3)*(a1-c1)/(3d0*b1)

	Ftot=Ftot-Im1_norm*((1d0+r1/(r2-r1))*Int1-Int2/((r2-r1)))
	Ftot=Ftot-Im2_norm*(Int2-r1*Int1)/((r2-r1))

c and the other half

	phi1=-image%phi(j)-0.5d0*pi/real(image%nphi)
	phi2=-image%phi(j)+0.5d0*pi/real(image%nphi)

	c1=cos(theta+phi1)
	b1=sin((phi1+phi2)/2d0+theta)
	a1=c1+2d0*sin((phi1-phi2)/2d0)*b1

	if(abs(k*r2).gt.0.1d0.or.abs(k*r1).gt.0.1d0) then
	expa1r1=cqexp(-cI*k*r1*a1)
	expa1r2=cqexp(-cI*k*r2*a1)
	expc1r1=cqexp(-cI*k*r1*c1)
	expc1r2=cqexp(-cI*k*r2*c1)
	Int1=(-1d0/(k**2*b1)) * ((expa1r2-expa1r1)/a1-
     &						 (expc1r2-expc1r1)/c1)
	Int2=(-1d0/(cI*k**3*b1)) * 
     &	((1d0+cI*k*a1*r2)*expa1r2/a1**2-
     &	 (1d0+cI*k*c1*r2)*expc1r2/c1**2-
     &	 (1d0+cI*k*a1*r1)*expa1r1/a1**2+
     &	 (1d0+cI*k*c1*r1)*expc1r1/c1**2)
	else
	Int1=(r2**2-r1**2)*(a1-c1)/(2d0*b1)
     & -(cI*(r2**3-r1**3)*(a1**2-c1**2)/(6d0*b1))*k
     & -((r2**4-r1**4)*(a1**3-c1**3)/(24d0*b1))*k**2
     & +(cI*(r2**5-r1**5)*(a1**4-c1**4)/(120d0*b1))*k**3
     & +((r2**6-r1**6)*(a1**5-c1**5)/(720d0*b1))*k**4
     & -(cI*(r2**7-r1**7)*(a1**6-c1**6)/(5040d0*b1))*k**5
     & -((r2**8-r1**8)*(a1**7-c1**7)/(40320d0*b1))*k**6
	Int2=(r2**3-r1**3)*(a1-c1)/(3d0*b1)
     & -(cI*(r2**4-r1**4)*(a1**2-c1**2)/(8d0*b1))*k
     & -((r2**5-r1**5)*(a1**3-c1**3)/(30d0*b1))*k**2
     & +(cI*(r2**6-r1**6)*(a1**4-c1**4)/(144d0*b1))*k**3
     & +((r2**7-r1**7)*(a1**5-c1**5)/(840d0*b1))*k**4
     & -(cI*(r2**8-r1**8)*(a1**6-c1**6)/(5760d0*b1))*k**5
     & -((r2**9-r1**9)*(a1**7-c1**7)/(45360d0*b1))*k**6
	endif


	cV=cV+Im1*((1d0+r1/(r2-r1))*Int1-Int2/((r2-r1)))
	cV=cV+Im2*(Int2-r1*Int1)/((r2-r1))

	Int1=(r2**2-r1**2)*(a1-c1)/(2d0*b1)
	Int2=(r2**3-r1**3)*(a1-c1)/(3d0*b1)

	Ftot=Ftot-Im1_norm*((1d0+r1/(r2-r1))*Int1-Int2/((r2-r1)))
	Ftot=Ftot-Im2_norm*(Int2-r1*Int1)/((r2-r1))

1	continue

	enddo
	enddo

	V=cqabs(cV/Ftot)
	phase = ATAN2 ( AIMAG(cV/Ftot), REAL(cV/Ftot) ) ! Gijsexp: complex phase
	
	return
	end
	
