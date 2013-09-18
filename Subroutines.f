c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	real*8 function Planck(T,lam)
	IMPLICIT NONE
	real*8 T,k,c,h,nu,lam,x
	k=1.3807d-16
	c=2.9979d10
	h=6.6261d-27
	nu=c/(lam*1d-4)
	x=h*nu/(k*T)
	Planck=(2d0*h*nu**3/c**2)/(exp(x)-1d0)
c	Planck=Planck*1e23

	return
	end

	real*8 function Luminosity(T,R)
	IMPLICIT NONE
	real*8 T,k,c,h,R,pi
	k=1.3807d-16
	c=2.9979d10
	h=6.6261d-27
	pi=3.1415926536

	Luminosity=(2d0*pi**4*k**4/(15d0*h**3*c**2))*T**4*pi*R**2
	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine regrid(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),y(n),x0,y0,x1,y1
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,RECL=100)
	i=1
1	read(20,*,end=102,err=1) x0,y0
103	if(x0.ge.grid(i)) then
		y(i)=y0
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y1
101	if(grid(i).le.x1.and.grid(i).ge.x0) then
		y(i)=y1+(grid(i)-x1)*(y0-y1)/(x0-x1)
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	do j=i,n
		y(j)=y(i-1)*grid(i-1)/grid(j)
	enddo
	close(unit=20)
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine regridlog(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),y(n),x0,y0,x1,y1
	real*8 lx0,ly0,lx1,ly1,lx,ly
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,RECL=100)
	i=1
1	read(20,*,end=102,err=1) x0,y0
103	if(x0.ge.grid(i)) then
		y(i)=y0
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y1
	if(y1.le.0d0) y1=y0*1d-50
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		y(i)=10d0**(log10(y1)+(log10(grid(i)/x1))*(log10(y0/y1))/(log10(x0/x1)))
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	do j=i,n
		y(j)=y(i-1)*grid(i-1)/grid(j)
	enddo
	close(unit=20)
	do i=1,n
		if(y(i).le.0d0) y(i)=y0*1d-60
	enddo
	return
	end

c-----------------------------------------------------------------------
c The new readstar subroutine uses a boxcar filtering to read in 
c high resolution spectra.
c-----------------------------------------------------------------------

	subroutine readstar(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n,nls
	real*8 grid(n),y(n),x0,y0,xedge(n+1)
	real*8 grid2(n),y2(n),tot(n)
	real*8,allocatable :: ls(:),Fs(:),dls(:)
	character*500 input
	logical truefalse,done(n)

	do i=1,n-1
		xedge(i+1)=sqrt(grid(i)*grid(i+1))
	enddo
	xedge(1)=grid(1)**2/xedge(2)
	xedge(n)=grid(n)**2/xedge(n-1)

	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,RECL=100)
	j=1
1	read(20,*,end=2,err=1) x0,y0
	j=j+1
	goto 1
2	continue
	nls=j-1
	close(unit=20)
	allocate(ls(nls))
	allocate(dls(nls))
	allocate(Fs(nls))

	open(unit=20,file=input,RECL=100)
	do i=1,nls
3		read(20,*,end=4,err=3) ls(i),Fs(i)
	enddo
4	close(unit=20)
	do i=2,nls-1
		dls(i)=(1d0/ls(i-1)-1d0/ls(i+1))/2d0
	enddo
	dls(1)=dls(2)
	dls(nls)=dls(nls-1)
	
	j=1
	tot=0d0
	y=0d0
	done=.false.
	do i=1,nls
5		continue
		if(ls(i).lt.xedge(j)) goto 6
		if(ls(i).gt.xedge(j+1)) then
			j=j+1
			if(j.gt.n) goto 7
			goto 5
		endif
		y(j)=y(j)+Fs(i)*dls(i)
		tot(j)=tot(j)+dls(i)
		done(j)=.true.
6		continue
	enddo
7	continue
	j=0
	do i=1,n
		if(done(i)) then
			y(i)=y(i)/tot(i)
		else
			j=j+1
			grid2(j)=grid(i)
		endif
	enddo
	if(j.ne.0) then
		call readstar_interpol(input,grid2,y2,j)
		j=0
		do i=1,n
			if(.not.done(i)) then
				j=j+1
				y(i)=y2(j)
			endif
		enddo
	endif
	
	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine readstar_interpol(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),y(n),x0,y0,x1,y1
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,RECL=100)
	i=1
1	read(20,*,end=102,err=1) x0,y0
103	if(x0.ge.grid(i)) then
		y(i)=y0
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y1
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		if(y1.gt.1d-60.and.y0.gt.1d-60) then
			y(i)=10d0**(log10(y1)+(log10(grid(i))-log10(x1))*(log10(y0)-log10(y1))/(log10(x0)-log10(x1)))
		else
			y(i)=y1+(grid(i)-x1)*(y0-y1)/(x0-x1)
		endif
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	do j=i,n
		y(j)=y1*x1**2/grid(j)**2
	enddo
	close(unit=20)
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine regrid2(input,grid,y,yerr,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),y(n),yerr(n),x0,y0,x1,y1,yerr0,yerr1
	character*500 input
	character*500 line1
	logical truefalse,sqrtflux
	inquire(file=input,exist=truefalse)
	sqrtflux=.false.
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,RECL=1000)
	i=1
	read(20,'(a500)') line1(1:500)
	read(line1,*,end=300) x0,y0,yerr0
	goto 103
300	print*,'No errors, using sqrt(flux)'
	sqrtflux=.true.
	read(line1,*) x0,y0
	yerr0=sqrt(abs(y0))
103	if(x0.ge.grid(i)) then
		y(i)=y0
		yerr(i)=yerr0
		i=i+1
		goto 103
	endif
100	if(.not.sqrtflux) then
		read(20,*,end=102) x1,y1,yerr1
	else
		read(20,*,end=102) x1,y1
		yerr1=sqrt(abs(y1))
	endif	
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		y(i)=y1+(grid(i)-x1)*(y0-y1)/(x0-x1)
		yerr(i)=yerr1+(grid(i)-x1)*(yerr0-yerr1)/(x0-x1)
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	yerr0=yerr1
	goto 100
102	continue
	do j=i,n
		y(j)=y(i-1)*grid(i-1)/grid(j)
		yerr(j)=yerr(i-1)*grid(i-1)/grid(j)
	enddo
	close(unit=20)
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine readrefind(input,grid,e1,e2,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),e1(n),e2(n),x0,y01,y02,x1,y11,y12,wp,gamma
	character*500 input
	open(unit=20,file=input,RECL=1000)
	i=1
1	read(20,*,end=102,err=1) x0,y01,y02
	wp=(1d0-y01)/x0**2
	gamma=y02/x0**3
103	if(x0.ge.grid(i)) then
		e1(i)=1d0-wp*grid(i)**2
		e2(i)=gamma*grid(i)**3
c		e1(i)=y01
c		e2(i)=y02
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y11,y12
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		e1(i)=y11+(grid(i)-x1)*(y01-y11)/(x0-x1)
		e2(i)=y12+(grid(i)-x1)*(y02-y12)/(x0-x1)
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y01=y11
	y02=y12
	goto 100
102	continue
	do j=i,n
		e1(j)=e1(i-1)
		e2(j)=e2(i-1)*grid(i-1)/grid(j)
	enddo
	close(unit=20)
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine q_cde(e1,e2,lam,rad,cab)
	IMPLICIT NONE
	real*8 e1,e2,lam,rad,cab,V,k,pi
	complex*16 m,alpha
	pi=3.1415926536
	m=dcmplx(e1,e2)
	alpha=2d0*m*m*log(m*m)/(m*m-1d0)-2d0
	V=4d0*pi*rad**3/3d0
	k=2d0*pi/lam
	cab=k*V*dimag(alpha)
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine ignorestar(un)
	IMPLICIT NONE
	integer un
	character c
1	read(un,fmt=3,end=2) c
	if(c.eq.'*') goto 1
	backspace(unit=un)
2	continue
3	format(a1)
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine dosmooth(lam,Cabs,nlam,rsmooth)
	IMPLICIT  NONE
	integer nlam
	real*8 rsmooth,Cabs(nlam),lam(nlam)
	real*8 newCabs(nlam),profile,dx,tot,prof
	integer i,j

	do j=1,nlam
	call tellertje(j,nlam)
	prof=profile(lam(1),lam(j),rsmooth)
	dx=(lam(2)-lam(1))/2.
	newCabs(j)=Cabs(1)*prof*dx
	tot=prof*dx
	do i=2,nlam-1
		prof=profile(lam(i),lam(j),rsmooth)
		dx=(lam(i+1)-lam(i-1))/2.
		newCabs(j)=newCabs(j)+Cabs(i)*prof*dx
		tot=tot+prof*dx
	enddo
	prof=profile(lam(nlam),lam(j),rsmooth)
	dx=(lam(nlam)-lam(nlam-1))/2.
	newCabs(j)=newCabs(j)+Cabs(nlam)*prof*dx
	tot=tot+prof*dx
	newCabs(j)=newCabs(j)/tot
	enddo

	do i=1,nlam
		Cabs(i)=newCabs(i)
	enddo

	return
	end

	real*8 function profile(x,x0,R)
	IMPLICIT NONE
	real*8 x,x0,R

c	Gamma distribution profile
	profile=exp(-((x-x0)*R*2d0/x0)**2)

c	Block profile
c	profile=0d0
c	if(abs(x-x0).lt.(x0/R/2d0)) profile=1d0

	return
	end

	subroutine ConvOld(im,IMDIM,lam,D,fov0,width0)
	IMPLICIT NONE
	integer n,IMDIM
	complex*16,allocatable :: iml(:,:),psf(:,:),seeing(:,:),imh(:,:)
	real*8 bessj1,x,y,r,pi,rscale,im(IMDIM,IMDIM),tot,tot2,max,min,max0
	integer nn(2),ndim,i,j,ix,iy,nsplit,isplit,isign,ii,jj
	parameter(pi=3.1415926536)
	real*8 lam,D,fov0,fov,width,width0,phi

	max0=1d7
	
	n=1
1	n=n*2
	if(n.lt.IMDIM*4) goto 1

	allocate(psf(n,n))
	allocate(seeing(n,n))
	allocate(iml(n,n))

	min=1d200
	do i=1,IMDIM
	do j=1,IMDIM
		if(im(i,j).gt.0d0.and.im(i,j).lt.min) min=im(i,j)
	enddo
	enddo
	max=1d0
	do i=1,IMDIM
	do j=1,IMDIM
		if(im(i,j)/min.gt.max) max=im(i,j)/min
	enddo
	enddo

	iml=0d0
	if(max.gt.max0) then
		allocate(imh(n,n))
		imh=0d0
	endif
	do i=1,IMDIM
	do j=1,IMDIM
		if(im(i,j)/min.le.max0) then
			iml(i,j)=im(i,j)
		else
			imh(i,j)=im(i,j)/max0
		endif
	enddo
	enddo

	fov=fov0*real(n)/real(IMDIM)
	width=width0*real(n)/real(IMDIM)

	rscale=real(n)*lam/(pi*D*1d6)*(3600d0*180d0/pi)/fov

	psf=0d0
	seeing=0d0

	tot=0d0
	tot2=0d0
	do i=1,IMDIM
	call tellertje(i,IMDIM)
	do j=1,IMDIM
		if(i.eq.1.and.j.eq.1) then
			psf(i,j)=1d0
		else
		x=real(i-1)
		y=real(j-1)
		r=sqrt(x**2+y**2)/rscale
		if(r.lt.50d0) then
			nsplit=int((5d0/rscale))+1
		else if(r.lt.500d0) then
			nsplit=int((1d0/rscale))+1
		else
			nsplit=1
		endif

		do ii=1,nsplit
		do jj=1,nsplit

		x=real(i-1)+real(ii-1)/real(nsplit)
		y=real(j-1)+real(jj-1)/real(nsplit)
		r=sqrt(x**2+y**2)/rscale

		psf(i,j)=psf(i,j)+(2d0*bessj1(r)/r)**2/real(nsplit**2)
c		psf(i,j)=psf(i,j)+(exp(-r**2)/2d0+3d0/(r**3+6d0))/real(nsplit**2)
		
		enddo
		enddo
		endif

		psf(n+1-i,j)=psf(i,j)
		psf(n+1-i,n+1-j)=psf(i,j)
		psf(i,n+1-j)=psf(i,j)

		tot=tot+psf(i,j)+psf(n+1-i,j)+psf(n+1-i,n+1-j)+psf(i,n+1-j)

		x=real(i-1)
		y=real(j-1)
		r=sqrt(x**2+y**2)*fov/real(n)
		seeing(i,j)=exp(-(r/width)**2)
		seeing(n+1-i,j)=seeing(i,j)
		seeing(n+1-i,n+1-j)=seeing(i,j)
		seeing(i,n+1-j)=seeing(i,j)
		tot2=tot2+seeing(i,j)*4d0
	enddo
	enddo

	psf=psf/tot
	seeing=seeing/tot2

	ndim=2
	nn(1:2)=n
	isign=1
	call fourn(iml,nn,ndim,isign)
	call fourn(psf,nn,ndim,isign)
	call fourn(seeing,nn,ndim,isign)

	do i=1,n
	do j=1,n
		iml(i,j)=iml(i,j)*psf(i,j)*seeing(i,j)
	enddo
	enddo

	ndim=2
	nn(1:2)=n
	isign=-1
	call fourn(iml,nn,ndim,isign)

	do i=1,IMDIM
	do j=1,IMDIM
		im(i,j)=real(iml(i,j))
	enddo
	enddo

	if(max.gt.max0) then
		min=1d200
		do i=1,IMDIM
		do j=1,IMDIM
			if(real(imh(i,j)).gt.0d0.and.real(imh(i,j)).lt.min) min=real(imh(i,j))
		enddo
		enddo
		max=1d0
		do i=1,IMDIM
		do j=1,IMDIM
			if(real(imh(i,j))/min.gt.max) max=real(imh(i,j))/min
		enddo
		enddo
		iml=imh
		imh=0d0
		do i=1,IMDIM
		do j=1,IMDIM
			if(real(iml(i,j))/min.gt.max0) then
				imh(i,j)=iml(i,j)/max0**2
				iml(i,j)=0d0
			endif
		enddo
		enddo
		isign=1
		call fourn(iml,nn,ndim,isign)
		if(max.gt.max0) call fourn(imh,nn,ndim,isign)
		do i=1,n
		do j=1,n
			iml(i,j)=iml(i,j)*psf(i,j)*seeing(i,j)
			if(max.gt.max0) imh(i,j)=imh(i,j)*psf(i,j)*seeing(i,j)
		enddo
		enddo
		isign=-1
		call fourn(iml,nn,ndim,isign)
		if(max.gt.max0) call fourn(imh,nn,ndim,isign)
		do i=1,IMDIM
		do j=1,IMDIM
			im(i,j)=im(i,j)+real(iml(i,j))*max0
			if(max.gt.max0) im(i,j)=im(i,j)+real(imh(i,j))*max0**3
		enddo
		enddo
	endif

	im=im/(real(n*n))

	deallocate(psf)
	deallocate(seeing)
	deallocate(iml)
	if(allocated(imh)) deallocate(imh)

	return
	end



	subroutine Convolution(im,imQ,imU,imV,IMDIM,lam0,Diam,Diam2,SpW,mask,masksize,owa,strehl,fov0,width0,snoise)
	use Parameters
	IMPLICIT NONE
	integer n,IMDIM
	complex*16,allocatable :: psf(:,:),image(:,:),seeing(:,:),peakseeing(:,:)
	complex*16,allocatable :: imstar(:,:),psfcor(:,:)
	complex*16 ic
	parameter(ic=(0d0,1d0))
	real*8 bessj1,x,y,r,rscale,im(IMDIM,IMDIM),tot,tot2,masksize,strehl,iwa,owa
	real*8 imQ(IMDIM,IMDIM),imU(IMDIM,IMDIM),imV(IMDIM,IMDIM),fcoro,tot3,width_ao
	integer nn(2),ndim,i,j,ix,iy,nsplit,isplit,isign,ii,jj,nsplitpix,n2
	real*8 lam0,Diam,fov0,fov,width,width0,phi,rw,Diam2,SpW,dx,mask,psf2
	integer*8 plan
	integer FFTW_ESTIMATE
	parameter(FFTW_ESTIMATE=64)
	real*8 gasdev,ran2,snoise
	character*500 filename
	complex*16 ctot1,ctot2
	complex*16,allocatable :: CimI(:,:),CimQ(:,:),CimU(:,:)

	nsplitpix=0
	fov=fov0

1	nsplitpix=nsplitpix+1
	n=IMDIM*nsplitpix
	rscale=real(n)*lam0/(pi*1d6)*(3600d0*180d0/pi)/fov
	x=real(n/4)*rscale*2d0*pi/real(n)
	if(x.lt.Diam.and.IMDIM*(nsplitpix+1)*4.lt.NMAX_CONVOLUTION) goto 1
	if(mask.lt.1d0.and.(masksize*real(n)/fov).lt.10d0.and.IMDIM*(nsplitpix+1)*4.lt.NMAX_CONVOLUTION) goto 1
	if(rscale.lt.(Diam/pi).and.rscale.lt.((3600d0*180d0/pi**2)*lam0/1d6/width0).and.IMDIM*(nsplitpix+1)*4.lt.NMAX_CONVOLUTION) goto 1


	n=1
2	n=n*2
	fov=fov0*real(n)/real(IMDIM)/real(nsplitpix)
	rscale=real(n)*lam0/(pi*1d6)*(3600d0*180d0/pi)/fov
	if(n.lt.IMDIM*nsplitpix*4.and.n.lt.NMAX_CONVOLUTION) goto 2
	if((Diam/(rscale*2d0*pi/real(n))).lt.20d0.and.n.lt.NMAX_CONVOLUTION.and.Diam.gt.0d0) goto 2
	if((Diam2/(rscale*2d0*pi/real(n))).lt.2d0.and.n.lt.NMAX_CONVOLUTION.and.Diam2.gt.0d0) goto 2

	fov=fov0*real(n)/real(IMDIM)/real(nsplitpix)
	rscale=real(n)*lam0/(pi*1d6)*(3600d0*180d0/pi)/fov
	x=real(n/4)*rscale*2d0*pi/real(n)

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Adjusting pixel scale factor ",i3)') nsplitpix
	write(9,'("Adjusting pixel scale factor ",i3)') nsplitpix
	if(n.ge.10000) then
		write(*,'("FFT on grid of ",i5,"x",i5)') n,n
		write(9,'("FFT on grid of ",i5,"x",i5)') n,n
	else if(n.ge.1000) then
		write(*,'("FFT on grid of ",i4,"x",i4)') n,n
		write(9,'("FFT on grid of ",i4,"x",i4)') n,n
	else if(n.ge.100) then
		write(*,'("FFT on grid of ",i3,"x",i3)') n,n
		write(9,'("FFT on grid of ",i3,"x",i3)') n,n
	else
		write(*,'("FFT on grid of ",i2,"x",i2)') n,n
		write(9,'("FFT on grid of ",i2,"x",i2)') n,n
	endif
	write(*,'("Telescope aperture sampled with ",i5," pixels")') int(Diam/(rscale*2d0*pi/real(n)))
	write(9,'("Telescope aperture sampled with ",i5," pixels")') int(Diam/(rscale*2d0*pi/real(n)))

	if(n.lt.IMDIM*nsplitpix*4) then
		write(*,'("Warning!! cannot do the convolution properly!!")')
		write(9,'("Warning!! cannot do the convolution properly!!")')
	endif

	allocate(seeing(n,n))
	allocate(psf(n,n))
	allocate(image(n,n))

	if(mask.lt.1d0) then
		allocate(psfcor(n,n))
		allocate(imstar(n,n))
	endif
	width=width0

	fcoro=1d0
	if(width.gt.0d0) then
	fcoro=0d0
	seeing=0d0
	tot2=0d0
	if(strehl.gt.0d0) then
		allocate(peakseeing(n,n))
		peakseeing=0d0
		tot3=0d0
		iwa=(3600d0*180d0*(1.22d0*lam0*1d-6/Diam)/pi)
		width_ao=iwa*sqrt(1d0/(-log(1d0-strehl)))
	endif
	n2=IMDIM*nsplitpix
	if(n2.gt.(n/2)) n2=n/2
	do i=1,n2
	do j=1,n2
		if(i.eq.1.and.j.eq.1) then
			seeing(i,j)=1d0
		else
		x=real(i-1)
		y=real(j-1)
		r=sqrt(x**2+y**2)/rscale
		nsplit=int((5d0/rscale))+10

		do ii=1,nsplit
		do jj=1,nsplit

		x=real(i-1)+real(ii-1)/real(nsplit)
		y=real(j-1)+real(jj-1)/real(nsplit)
		r=sqrt(x**2+y**2)
		rw=r*fov/real(n)

		if(strehl.gt.0d0) then
			if(rw.lt.iwa) then
				peakseeing(i,j)=peakseeing(i,j)+exp(-(rw/width_ao)**2)/real(nsplit*nsplit)
				if(rw.le.masksize) then
					fcoro=fcoro+exp(-(rw/width_ao)**2)/real(nsplit*nsplit)
				endif
			else
				if(rw.gt.owa) then
					seeing(i,j)=seeing(i,j)+((1d0+((rw-owa)/width)**2)**(-2.5))/real(nsplit*nsplit)
				else
					seeing(i,j)=seeing(i,j)+((1d0+(0.5d0*(rw-owa)/width)**2)**(-2.5))/real(nsplit*nsplit)
				endif
			endif
		else
			seeing(i,j)=seeing(i,j)+((1d0+(rw/width)**2)**(-2.5))/real(nsplit*nsplit)
			if(rw.le.masksize) then
				fcoro=fcoro+((1d0+(rw/width)**2)**(-2.5))/real(nsplit*nsplit)
			endif
		endif
		
		enddo
		enddo
		endif

		tot2=tot2+seeing(i,j)
		seeing(n+1-i,j)=seeing(i,j)
		tot2=tot2+seeing(n+1-i,j)
		seeing(n+1-i,n+1-j)=seeing(i,j)
		tot2=tot2+seeing(n+1-i,n+1-j)
		seeing(i,n+1-j)=seeing(i,j)
		tot2=tot2+seeing(i,n+1-j)

		if(strehl.gt.0d0) then
			tot3=tot3+peakseeing(i,j)
			peakseeing(n+1-i,j)=peakseeing(i,j)
			tot3=tot3+peakseeing(n+1-i,j)
			peakseeing(n+1-i,n+1-j)=peakseeing(i,j)
			tot3=tot3+peakseeing(n+1-i,n+1-j)
			peakseeing(i,n+1-j)=peakseeing(i,j)
			tot3=tot3+peakseeing(i,n+1-j)
		endif
	enddo
	enddo


	if(strehl.gt.0d0) then
		seeing=peakseeing*strehl/tot3+(1d0-strehl)*seeing/tot2

		tot=0d0
		tot3=0d0
		do i=1,n/2
		do j=1,n/2
			x=real(i-1)
			y=real(j-1)
			r=sqrt(x**2+y**2)
			rw=r*fov/real(n)
			tot=tot+seeing(i,j)
			if(rw.lt.iwa) tot3=tot3+seeing(i,j)
		enddo
		enddo
		print*,'Strehl:',tot3/tot

		deallocate(peakseeing)
		fcoro=strehl*4d0*(fcoro+1d0)/tot3
		if(fcoro.gt.1d0) fcoro=1d0
	else
		seeing=seeing/tot2
		fcoro=4d0*(fcoro+1d0)/tot2
		if(fcoro.gt.1d0) fcoro=1d0
	endif

	if(snoise.gt.0d0) then
		do i=1,n
		do j=1,n
			psf(i,j)=gasdev(idum)
		enddo
		enddo

		isign=1
		call dfftw_plan_dft_2d(plan, n,n, psf,psf,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, psf, psf)
		call dfftw_destroy_plan(plan)

		psf=psf/real(n*n)

		ctot1=0d0
		do i=1,n
		do j=1,n
			ctot1=ctot1+psf(i,j)
		enddo
		enddo

		ctot2=snoise
c		psf=ctot1+(psf-ctot1)*ctot2
		psf=psf**ctot2
		psf=psf/cdabs(psf)

		seeing=seeing*psf
	
		tot=0d0
		do i=1,n
		do j=1,n
			tot=tot+seeing(i,j)
		enddo
		enddo
		seeing=seeing/tot
	endif
	
	isign=1
	call dfftw_plan_dft_2d(plan, n,n, seeing,seeing,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, seeing, seeing)
	call dfftw_destroy_plan(plan)
	
c	write(filename,'("seeing.fits")')
c	call writefitsfile(filename,seeing(1:IMDIM,1:IMDIM),IMDIM)

	else
	seeing=1d0
	endif

	if(Diam.gt.0d0) then
	psf=0d0
	tot=0d0
	do i=1,n/2
	do j=1,n/2
		if(i.eq.1.and.j.eq.1.and.Diam2.le.0d0) then
			psf(i,j)=1d0
		else
		x=real(i-1)
		y=real(j-1)
		r=sqrt(x**2+y**2)/rscale
		nsplit=int((5d0/rscale))+10

		dx=rscale*2d0*pi/real(n)/real(nsplit)
		psf2=0d0
		do ii=1,nsplit
		do jj=1,nsplit

		x=real(i-1)+real(ii-1)/real(nsplit)
		y=real(j-1)+real(jj-1)/real(nsplit)
		r=sqrt(x**2+y**2)
		rw=r*fov/real(n)
		r=r*rscale*2d0*pi/real(n)
		x=x*rscale*2d0*pi/real(n)
		y=y*rscale*2d0*pi/real(n)

		if(r.le.Diam.and.r.ge.Diam2) then
			if(abs(x).ge.(SpW).and.abs(y).ge.(SpW)) then
				psf(i,j)=psf(i,j)+1d0/real(nsplit**2)
			endif
			if(abs(x).lt.(SpW).and.(abs(x)+dx).ge.(SpW).and.abs(y).ge.(SpW)) then
				psf(i,j)=psf(i,j)+(1d0/real(nsplit**2))*(abs(x)+dx-SpW)/dx
			endif
			if(abs(y).lt.(SpW).and.(abs(y)+dx).ge.(SpW).and.abs(x).ge.(SpW)) then
				psf(i,j)=psf(i,j)+(1d0/real(nsplit**2))*(abs(y)+dx-SpW)/dx
			endif
			if(cdabs(psf(i,j)).gt.1d0) psf(i,j)=1d0
		endif
		
		enddo
		enddo
		endif

		tot=tot+psf(i,j)
		psf(n+1-i,j)=psf(i,j)
		tot=tot+psf(n+1-i,j)
		psf(n+1-i,n+1-j)=psf(i,j)
		tot=tot+psf(n+1-i,n+1-j)
		psf(i,n+1-j)=psf(i,j)
		tot=tot+psf(i,n+1-j)
	enddo
	enddo

	if(mask.lt.1d0) then
		psfcor=psf/tot
		isign=-1
		call dfftw_plan_dft_2d(plan, n,n, psfcor,psfcor,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, psfcor, psfcor)
		call dfftw_destroy_plan(plan)
		image=1d0
		do i=1,n
		do j=1,n
			nsplit=10
			do ii=1,nsplit
			do jj=1,nsplit
				if(i.lt.n/2) then
					x=real(i-1)+real(ii-1)/real(nsplit)
				else
					x=real(i-n)-real(ii-1)/real(nsplit)
				endif
				if(j.lt.n/2) then
					y=real(j-1)+real(jj-1)/real(nsplit)
				else
					y=real(j-n)-real(ii-1)/real(nsplit)
				endif
				x=x*fov/real(n)
				y=y*fov/real(n)
				r=sqrt(x**2+y**2)
				if(r.le.masksize) then
					image(i,j)=image(i,j)-(1d0-mask)/real(nsplit*nsplit)
				endif
			enddo
			enddo
		enddo
		enddo

		psfcor=psfcor*image

		isign=1
		call dfftw_plan_dft_2d(plan, n,n, psfcor,psfcor,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, psfcor, psfcor)
		call dfftw_destroy_plan(plan)
		psfcor=psfcor*psf
		isign=-1
		call dfftw_plan_dft_2d(plan, n,n, psfcor,psfcor,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, psfcor, psfcor)
		call dfftw_destroy_plan(plan)
		do i=1,n
		do j=1,n
			psfcor(i,j)=cdabs(psfcor(i,j))**2
		enddo
		enddo
		isign=1
		call dfftw_plan_dft_2d(plan, n,n, psfcor,psfcor,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, psfcor, psfcor)
		call dfftw_destroy_plan(plan)
		psfcor=psfcor/real(n*n)/real(n*n)/real(n*n)
	endif

	psf=psf/tot
	ndim=2
	nn(1:2)=n

	isign=-1
	call dfftw_plan_dft_2d(plan, n,n, psf,psf,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, psf, psf)
	call dfftw_destroy_plan(plan)

	tot=0d0
	do i=1,n
	do j=1,n
		psf(i,j)=cdabs(psf(i,j))**2
		if(i.gt.IMDIM*nsplitpix.and.i.lt.(n-IMDIM*nsplitpix)) psf(i,j)=0d0
		if(j.gt.IMDIM*nsplitpix.and.j.lt.(n-IMDIM*nsplitpix)) psf(i,j)=0d0
		if(mask.lt.1d0) then
			if(i.gt.IMDIM*nsplitpix.and.i.lt.(n-IMDIM*nsplitpix)) psfcor(i,j)=0d0
			if(j.gt.IMDIM*nsplitpix.and.j.lt.(n-IMDIM*nsplitpix)) psfcor(i,j)=0d0
		endif
		tot=tot+psf(i,j)
	enddo
	enddo
	psf=psf/tot
	if(mask.lt.1d0) psfcor=psfcor/tot

	isign=1

	call dfftw_plan_dft_2d(plan, n,n, psf,psf,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, psf, psf)
	call dfftw_destroy_plan(plan)

	psf=psf/real(n*n)

	else
	psf=1d0/real(n*n)
	endif

	image=0d0
	if(mask.lt.1d0) imstar=0d0
	do i=1,IMDIM
	do j=1,IMDIM
		x=real(i)-real(IMDIM+1)/2d0
		y=real(j)-real(IMDIM+1)/2d0
		if(IMDIM.eq.(2*(IMDIM/2))) then
			x=abs(x)-0.5d0
			y=abs(y)-0.5d0
		endif
c		x=(x*fov0/real(IMDIM))*D%distance*AU/parsec
c		y=(y*fov0/real(IMDIM))*D%distance*AU/parsec
		r=sqrt(x**2+y**2)*fov/real(n)
		if(r.gt.iwa.or.mask.ge.1d0) then
			do ii=1,nsplitpix
			do jj=1,nsplitpix
				image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj)=im(i,j)/real(nsplitpix*nsplitpix)
			enddo
			enddo
		else
			do ii=1,nsplitpix
			do jj=1,nsplitpix
				image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj)=(1d0-fcoro)*im(i,j)/real(nsplitpix*nsplitpix)
				imstar((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj)=fcoro*im(i,j)/real(nsplitpix*nsplitpix)
			enddo
			enddo
		endif
	enddo
	enddo
	isign=1
	call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, image, image)
	call dfftw_destroy_plan(plan)
	image=image*psf*seeing
	isign=-1
	call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, image, image)
	call dfftw_destroy_plan(plan)
	if(mask.lt.1d0) then
		isign=1
		call dfftw_plan_dft_2d(plan, n,n, imstar,imstar,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, imstar, imstar)
		call dfftw_destroy_plan(plan)
		imstar=imstar*psfcor*seeing
		isign=-1
		call dfftw_plan_dft_2d(plan, n,n, imstar,imstar,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, imstar, imstar)
		call dfftw_destroy_plan(plan)
	endif

	allocate(CimI(IMDIM,IMDIM))

	do i=1,IMDIM
	do j=1,IMDIM
		CimI(i,j)=0d0
		do ii=1,nsplitpix
		do jj=1,nsplitpix
			CimI(i,j)=CimI(i,j)+(image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj))
			if(mask.lt.1d0) then
				CimI(i,j)=CimI(i,j)+(imstar((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj))
			endif
		enddo
		enddo
	enddo
	enddo

	im=cdabs(CimI)

	if(scat_how.eq.2) then
		allocate(CimQ(IMDIM,IMDIM))
		allocate(CimU(IMDIM,IMDIM))
		image=0d0
		do i=1,IMDIM
		do j=1,IMDIM
			do ii=1,nsplitpix
			do jj=1,nsplitpix
				image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj)=imQ(i,j)/real(nsplitpix*nsplitpix)
			enddo
			enddo
		enddo
		enddo
		isign=1
c		call fourn(image,nn,ndim,isign)
		call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, image, image)
		call dfftw_destroy_plan(plan)
		image=image*psf*seeing
		isign=-1
c		call fourn(image,nn,ndim,isign)
		call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, image, image)
		call dfftw_destroy_plan(plan)
		do i=1,IMDIM
		do j=1,IMDIM
			CimQ(i,j)=0d0
			do ii=1,nsplitpix
			do jj=1,nsplitpix
				CimQ(i,j)=CimQ(i,j)+(image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj))
			enddo
			enddo
		enddo
		enddo
		
		imQ=(cdabs(CimI+CimQ)-cdabs(CimI-CimQ))/2d0

		image=0d0
		do i=1,IMDIM
		do j=1,IMDIM
			do ii=1,nsplitpix
			do jj=1,nsplitpix
				image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj)=imU(i,j)/real(nsplitpix*nsplitpix)
			enddo
			enddo
		enddo
		enddo
		isign=1
c		call fourn(image,nn,ndim,isign)
		call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, image, image)
		call dfftw_destroy_plan(plan)
		image=image*psf*seeing
		isign=-1
c		call fourn(image,nn,ndim,isign)
		call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, image, image)
		call dfftw_destroy_plan(plan)
		do i=1,IMDIM
		do j=1,IMDIM
			CimU(i,j)=0d0
			do ii=1,nsplitpix
			do jj=1,nsplitpix
				CimU(i,j)=CimU(i,j)+real(image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj))
			enddo
			enddo
		enddo
		enddo

		imU=(cdabs(CimI+CimU)-cdabs(CimI-CimU))/2d0

		image=0d0
		do i=1,IMDIM
		do j=1,IMDIM
			do ii=1,nsplitpix
			do jj=1,nsplitpix
				image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj)=imV(i,j)/real(nsplitpix*nsplitpix)
			enddo
			enddo
		enddo
		enddo
		isign=1
c		call fourn(image,nn,ndim,isign)
		call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, image, image)
		call dfftw_destroy_plan(plan)
		image=image*psf*seeing
		isign=-1
c		call fourn(image,nn,ndim,isign)
		call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, image, image)
		call dfftw_destroy_plan(plan)
		do i=1,IMDIM
		do j=1,IMDIM
			imV(i,j)=0d0
			do ii=1,nsplitpix
			do jj=1,nsplitpix
				imV(i,j)=imV(i,j)+real(image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj))
			enddo
			enddo
		enddo
		enddo
	endif


	deallocate(psf)
	deallocate(image)
	deallocate(seeing)
	if(mask.lt.1d0) then
		deallocate(psfcor)
		deallocate(imstar)
	endif
	deallocate(CimI)
	if(scat_how.eq.2) then
		deallocate(CimQ)
		deallocate(CimU)
	endif
	
	return
	end


	subroutine ReadFits(filename,im,n)
	IMPLICIT NONE
	character*500 filename
	integer n
	real*8 im(n,n)
!
!*******************************************************************************
!
!! READ_IMAGE reads a FITS image and determines the pixel range.
!
	integer status,unit,readwrite,blocksize,naxes(2),nfound
	real*8 nullval
	logical anynull
	integer group,firstpix,nbuffer,npixels,i
!
	status=0
!
!  Get an unused Logical Unit Number to use to open the FITS file.
!
	call ftgiou(unit,status)
!
!  open the FITS file previously created by WRITE_IMAGE.
!
	readwrite=0
	call ftopen(unit,filename,readwrite,blocksize,status)
!
!  determine the size of the image
!
	print*,filename
	call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
!
!  check that it found both NAXIS1 and NAXIS2 keywords
!
	print*,naxes(1),naxes(2)
	if (nfound /= 2)then
		print *,'READ_IMAGE failed to read the NAXISn keywords.'
c		return
	end if
!
!  initialize variables.
!
	npixels=naxes(1)*naxes(2)
	group=1
	firstpix=1
	nullval=-999

!
	call ftgpvd(unit,group,firstpix,npixels,nullval,im(1:naxes(1),1:naxes(2)),anynull,status)
!
!  close the file and free the unit number
!
	call ftclos(unit, status)
	call ftfiou(unit, status)
!
!  Check for any error, and if so print out error messages
!

	return
	end

	subroutine ConvolutionFILE(im,imQ,imU,imV,IMDIM,psffile)
	use Parameters
	IMPLICIT NONE
	integer n,IMDIM
	complex*16,allocatable :: image(:,:),psf(:,:)
	real*8 psftemp(IMDIM,IMDIM)
	real*8 bessj1,x,y,r,rscale,im(IMDIM,IMDIM),tot,tot2
	real*8 imQ(IMDIM,IMDIM),imU(IMDIM,IMDIM),imV(IMDIM,IMDIM)
	integer nn(2),ndim,i,j,ix,iy,nsplit,isplit,isign,ii,jj
	real*8 lam0,Diam,fov0,fov,width,width0,phi
	character*500 psffile

	n=1
1	n=n*2
	if(n.lt.IMDIM*2.and.n.lt.5000) goto 1

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	if(n.ge.10000) then
		write(*,'("FFT on grid of ",i5,"x",i5)') n,n
		write(9,'("FFT on grid of ",i5,"x",i5)') n,n
	else if(n.ge.1000) then
		write(*,'("FFT on grid of ",i4,"x",i4)') n,n
		write(9,'("FFT on grid of ",i4,"x",i4)') n,n
	else if(n.ge.100) then
		write(*,'("FFT on grid of ",i3,"x",i3)') n,n
		write(9,'("FFT on grid of ",i3,"x",i3)') n,n
	else
		write(*,'("FFT on grid of ",i2,"x",i2)') n,n
		write(9,'("FFT on grid of ",i2,"x",i2)') n,n
	endif

	if(n.lt.IMDIM*2) then
		write(*,'("Warning!! cannot do the convolution properly!!")')
		write(9,'("Warning!! cannot do the convolution properly!!")')
	endif

	allocate(psf(n,n))
	allocate(image(n,n))

	call ReadFITS(psffile,psftemp,IMDIM)

	tot=0d0
	psf=0d0
	do i=1,IMDIM/2
	do j=1,IMDIM/2
		psf(j,i)=psftemp(IMDIM/2+i,IMDIM/2+j)
		tot=tot+psf(j,i)
		psf(n+1-j,i)=psftemp(IMDIM/2+i,IMDIM/2+1-j)
		tot=tot+psf(n+1-j,i)
		psf(j,n+1-i)=psftemp(IMDIM/2+1-i,IMDIM/2+j)
		tot=tot+psf(j,n+1-i)
		psf(n+1-j,n+1-i)=psftemp(IMDIM/2+1-i,IMDIM/2+1-j)
		tot=tot+psf(n+1-j,n+1-i)
	enddo
	enddo

	psf=psf/tot

	ndim=2
	nn(1:2)=n
	isign=1
	call fourn(psf,nn,ndim,isign)

	image=0d0
	do i=1,IMDIM
	do j=1,IMDIM
		image(i,j)=im(i,j)
	enddo
	enddo
	call fourn(image,nn,ndim,isign)
	image=image*psf
	isign=-1
	call fourn(image,nn,ndim,isign)
	do i=1,IMDIM
	do j=1,IMDIM
		im(i,j)=real(image(i,j))/real(n*n)
	enddo
	enddo

	if(scat_how.eq.2) then
		image=0d0
		do i=1,IMDIM
		do j=1,IMDIM
			image(i,j)=imQ(i,j)
		enddo
		enddo
		isign=1
		call fourn(image,nn,ndim,isign)
		image=image*psf
		isign=-1
		call fourn(image,nn,ndim,isign)
		do i=1,IMDIM
		do j=1,IMDIM
			imQ(i,j)=real(image(i,j))/real(n*n)
		enddo
		enddo

		image=0d0
		do i=1,IMDIM
		do j=1,IMDIM
			image(i,j)=imU(i,j)
		enddo
		enddo
		isign=1
		call fourn(image,nn,ndim,isign)
		image=image*psf
		isign=-1
		call fourn(image,nn,ndim,isign)
		do i=1,IMDIM
		do j=1,IMDIM
			imU(i,j)=real(image(i,j))/real(n*n)
		enddo
		enddo

		image=0d0
		do i=1,IMDIM
		do j=1,IMDIM
			image(i,j)=imV(i,j)
		enddo
		enddo
		isign=1
		call fourn(image,nn,ndim,isign)
		image=image*psf
		isign=-1
		call fourn(image,nn,ndim,isign)
		do i=1,IMDIM
		do j=1,IMDIM
			imV(i,j)=real(image(i,j))/real(n*n)
		enddo
		enddo
	endif

	deallocate(psf)
	deallocate(image)

	return
	end




      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL*8 data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL*8 tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=(wr)*data(k2)-(wi)*data(k2+1)
                tempi=(wr)*data(k2+1)+(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END





	real*8 function psfdev(idum)
	IMPLICIT NONE
	integer idum
	real*8 bessj1,x,y,ran2,r

c	psfdev=2d0*sqrt(-log(ran2(idum)))
c	return

1	continue
	r=4d0*ran2(idum)
	if(r.gt.2d0) then
		r=(r-2d0)*2d0
		x=8d0/(4d0-r)
		y=ran2(idum)*4d0/x**2
	else
		x=r
		y=ran2(idum)
	endif
	if(y.gt.(2d0*bessj1(x))**2/x) goto 1

	psfdev=x

	return
	end



	real*8 function bessdev(idum)
	IMPLICIT NONE
	integer idum
	real*8 bessj1,x,y,ran2,r,bessapprox,gasdev

	bessdev=gasdev(idum)*sqrt(2d0)
	return

1	continue
	r=2d0*ran2(idum)
	if(r.gt.1d0) then
		r=(r-1d0)*4d0
		x=8d0/(4d0-r)
		y=ran2(idum)*4d0/x**2
	else
		x=r*2d0
		y=ran2(idum)
	endif

c	if(y.gt.(2d0*bessj1(x))**2/x) goto 1

	bessapprox=(exp(-x**2/3d0)/2d0+3d0/(x**3+6d0))*x
	if(y.gt.bessapprox) goto 1

	bessdev=x

	return
	end



      FUNCTION bessj1(x)
      REAL*8 bessj1,x
      REAL*8 ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     *s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     *99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.d0)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-2.356194491d0
        bessj1=sqrt(.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.d0,x)
      endif
      return
      END



      FUNCTION poidev(xm,idum)
      INTEGER idum
      REAL*8 poidev,xm,PI
      PARAMETER (PI=3.141592654)
CU    USES gammln,ran2
      REAL*8 alxm,em,g,oldm,sq,t,y,gammln,ran2
      SAVE alxm,g,oldm,sq
      DATA oldm /-1d0/
      if (xm.lt.12d0)then
        if (xm.ne.oldm) then
          oldm=xm
          g=exp(-xm)
        endif
        em=-1
        t=1d0
2       em=em+1d0
        t=t*ran2(idum)
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2d0*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1d0)
        endif
1       y=tan(PI*ran2(idum))
        em=sq*y+xm
        if (em.lt.0d0) goto 1
        em=int(em)
        t=0.9d0*(1d0+y**2)*exp(em*alxm-gammln(em+1d0)-g)
        if (ran2(idum).gt.t) goto 1
      endif
      poidev=em
      return
      END


      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END


      FUNCTION gasdev(idum)
      INTEGER idum
      REAL*8 gasdev
CU    USES ran1
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END


      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END

C      FUNCTION ran0(idum)
C      INTEGER idum,IA,IM,IQ,IR,MASK
C      REAL*8 ran0,AM
C      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
C     *MASK=123459876)
C      INTEGER k
C      idum=ieor(idum,MASK)
C      k=idum/IQ
C      idum=IA*(idum-k*IQ)-IR*k
C      if (idum.lt.0) idum=idum+IM
C      ran0=AM*idum
C      idum=ieor(idum,MASK)
C      return
C      END
C
C
C      FUNCTION ran3(idum)
C      INTEGER idum
C      INTEGER MBIG,MSEED,MZ
CC     REAL MBIG,MSEED,MZ
C      REAL*8 ran3,FAC
C      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
CC     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
C      INTEGER i,iff,ii,inext,inextp,k
C      INTEGER mj,mk,ma(55)
CC     REAL mj,mk,ma(55)
C      SAVE iff,inext,inextp,ma
C      DATA iff /0/
C      if(idum.lt.0.or.iff.eq.0)then
C        iff=1
C        mj=abs(MSEED-abs(idum))
C        mj=mod(mj,MBIG)
C        ma(55)=mj
C        mk=1
C        do 11 i=1,54
C          ii=mod(21*i,55)
C          ma(ii)=mk
C          mk=mj-mk
C          if(mk.lt.MZ)mk=mk+MBIG
C          mj=ma(ii)
C11      continue
C        do 13 k=1,4
C          do 12 i=1,55
C            ma(i)=ma(i)-ma(1+mod(i+30,55))
C            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
C12        continue
C13      continue
C        inext=0
C        inextp=31
C        idum=1
C      endif
C      inext=inext+1
C      if(inext.eq.56)inext=1
C      inextp=inextp+1
C      if(inextp.eq.56)inextp=1
C      mj=ma(inext)-ma(inextp)
C      if(mj.lt.MZ)mj=mj+MBIG
C      ma(inext)=mj
C      ran3=mj*FAC
C      return
C      END
C

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

	real*8 ran1
	ran2=ran1(idum)
	return

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END


	subroutine stati(dat,n,aver,var1,var2)
	IMPLICIT NONE
	integer n,i,n1,n2
	real*8 dat(n),aver,var1,var2
	aver=0d0
	var1=0d0
	var2=0d0
	n1=0
	n2=0
	do i=1,n
		aver=aver+dat(i)
	enddo
	aver=aver/real(n)
	do i=1,n
		if(dat(i).lt.aver) then
			var1=var1+(aver-dat(i))**2
			n1=n1+1
		else
			var2=var2+(aver-dat(i))**2
			n2=n2+1
		endif
	enddo
	var1=sqrt(var1/real(n1-1))
	var2=sqrt(var2/real(n2-1))
	if(n1.le.1) var1=0d0
	if(n2.le.1) var2=0d0
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotate(x,y,z,u,v,w,theta)
	IMPLICIT NONE
	real*8 x,y,z,u,v,w,yy(3),theta,inp
	real*8 cost,sint,u2,v2,w2
	cost=cos(theta)
	sint=sin(theta)
	u2=u*u
	v2=v*v
	w2=w*w

	inp=x*u+y*v+z*w
	yy(1)=u*inp
     & +(x*(v2+w2)-u*(v*y+w*z))*cost
     & +(v*z-w*y)*sint
	yy(2)=v*inp
     & +(y*(u2+w2)-v*(u*x+w*z))*cost
     & +(w*x-u*z)*sint
	yy(3)=w*inp
     & +(z*(u2+v2)-w*(u*x+v*y))*cost
     & +(u*y-v*x)*sint
	x=yy(1)
	y=yy(2)
	z=yy(3)
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotateZ(x,y,z,cost,sint)
	IMPLICIT NONE
	real*8 x,y,z,xx,yy,r
	real*8 cost,sint

	xx=x*cost-y*sint
	yy=y*cost+x*sint
	x=xx
	y=yy
	

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotatevec(x,u,theta)
	IMPLICIT NONE
	real*8 x(3),u(3),y(3),theta,inp,cost,sint
	cost=cos(theta)
	sint=sin(theta)
	
	inp=x(1)*u(1)+x(2)*u(2)+x(3)*u(3)
	y(1)=u(1)*inp
     & +(x(1)*(u(2)**2+u(3)**2)-u(1)*(u(2)*x(2)+u(3)*x(3)))*cost
     & +(u(2)*x(3)-u(3)*x(2))*sint
	y(2)=u(2)*inp
     & +(x(2)*(u(1)**2+u(3)**2)-u(2)*(u(1)*x(1)+u(3)*x(3)))*cost
     & +(u(3)*x(1)-u(1)*x(3))*sint
	y(3)=u(3)*inp
     & +(x(3)*(u(1)**2+u(2)**2)-u(3)*(u(1)*x(1)+u(2)*x(2)))*cost
     & +(u(1)*x(2)-u(2)*x(1))*sint
	x(1)=y(1)
	x(2)=y(2)
	x(3)=y(3)
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine tellertje(i,n)
	use NAG
	IMPLICIT NONE
	integer i,n,f
	
	if(i.eq.1) write(*,'("....................")')
	if(i.eq.1) write(9,'("....................")')
	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*real(i-1)/real(n).lt.real(f)
     &   .and.20d0*real(i+1)/real(n).gt.real(f)) then
		write(*,'(".",$)')
		call flush(6)
		write(9,'(".",$)')
		call flush(9)
	endif
	
	if(i.eq.n) write(*,*)
	if(i.eq.n) write(9,*)
c	if(i.ne.n) then
c		write(*,'(5a1,$)') (char(8),j=1,5)
c		write(*,'(5a1,$)') (' ',j=1,5)
c		write(*,'(5a1,$)') (char(8),j=1,5)
c	else
c		write(*,*)
c	endif


	return
	end

	
	subroutine icosahedron(x,y,z)
	IMPLICIT NONE
	real*8 phi
	parameter(phi=1.61803399)
	real*8 x(20),y(20),z(20)
	real*8 vx(12),vy(12),vz(12)
	real*8 d

	d=sqrt((2d0*phi+1d0)**2+phi**2)
	vx( 1)=-1d0
	vy( 1)=-phi
	vz( 1)= 0d0
	vx( 2)=-1d0
	vy( 2)= phi
	vz( 2)= 0d0
	vx( 3)= 0d0
	vy( 3)=-1d0
	vz( 3)=-phi
	vx( 4)= 0d0
	vy( 4)=-1d0
	vz( 4)= phi
	vx( 5)= 0d0
	vy( 5)= 1d0
	vz( 5)=-phi
	vx( 6)= 0d0
	vy( 6)= 1d0
	vz( 6)= phi
	vx( 7)= 1d0
	vy( 7)=-phi
	vz( 7)= 0d0
	vx( 8)= 1d0
	vy( 8)= phi
	vz( 8)= 0d0
	vx( 9)=-phi
	vy( 9)= 0d0
	vz( 9)=-1d0
	vx(10)=-phi
	vy(10)= 0d0
	vz(10)= 1d0
	vx(11)= phi
	vy(11)= 0d0
	vz(11)=-1d0
	vx(12)= phi
	vy(12)= 0d0
	vz(12)= 1d0

	x(1)=(vx(1)+vx(3)+vx(7))/d
	x(2)=(vx(1)+vx(3)+vx(9))/d
	x(3)=(vx(1)+vx(4)+vx(7))/d
	x(4)=(vx(1)+vx(4)+vx(10))/d
	x(5)=(vx(1)+vx(9)+vx(10))/d
	x(6)=(vx(2)+vx(5)+vx(8))/d
	x(7)=(vx(2)+vx(5)+vx(9))/d
	x(8)=(vx(2)+vx(6)+vx(8))/d
	x(9)=(vx(2)+vx(6)+vx(10))/d
	x(10)=(vx(2)+vx(9)+vx(10))/d
	x(11)=(vx(3)+vx(5)+vx(9))/d
	x(12)=(vx(3)+vx(5)+vx(11))/d
	x(13)=(vx(3)+vx(7)+vx(11))/d
	x(14)=(vx(4)+vx(6)+vx(10))/d
	x(15)=(vx(4)+vx(6)+vx(12))/d
	x(16)=(vx(4)+vx(7)+vx(12))/d
	x(17)=(vx(5)+vx(8)+vx(11))/d
	x(18)=(vx(6)+vx(8)+vx(12))/d
	x(19)=(vx(7)+vx(11)+vx(12))/d
	x(20)=(vx(8)+vx(11)+vx(12))/d

	y(1)=(vy(1)+vy(3)+vy(7))/d
	y(2)=(vy(1)+vy(3)+vy(9))/d
	y(3)=(vy(1)+vy(4)+vy(7))/d
	y(4)=(vy(1)+vy(4)+vy(10))/d
	y(5)=(vy(1)+vy(9)+vy(10))/d
	y(6)=(vy(2)+vy(5)+vy(8))/d
	y(7)=(vy(2)+vy(5)+vy(9))/d
	y(8)=(vy(2)+vy(6)+vy(8))/d
	y(9)=(vy(2)+vy(6)+vy(10))/d
	y(10)=(vy(2)+vy(9)+vy(10))/d
	y(11)=(vy(3)+vy(5)+vy(9))/d
	y(12)=(vy(3)+vy(5)+vy(11))/d
	y(13)=(vy(3)+vy(7)+vy(11))/d
	y(14)=(vy(4)+vy(6)+vy(10))/d
	y(15)=(vy(4)+vy(6)+vy(12))/d
	y(16)=(vy(4)+vy(7)+vy(12))/d
	y(17)=(vy(5)+vy(8)+vy(11))/d
	y(18)=(vy(6)+vy(8)+vy(12))/d
	y(19)=(vy(7)+vy(11)+vy(12))/d
	y(20)=(vy(8)+vy(11)+vy(12))/d
	
	z(1)=(vz(1)+vz(3)+vz(7))/d
	z(2)=(vz(1)+vz(3)+vz(9))/d
	z(3)=(vz(1)+vz(4)+vz(7))/d
	z(4)=(vz(1)+vz(4)+vz(10))/d
	z(5)=(vz(1)+vz(9)+vz(10))/d
	z(6)=(vz(2)+vz(5)+vz(8))/d
	z(7)=(vz(2)+vz(5)+vz(9))/d
	z(8)=(vz(2)+vz(6)+vz(8))/d
	z(9)=(vz(2)+vz(6)+vz(10))/d
	z(10)=(vz(2)+vz(9)+vz(10))/d
	z(11)=(vz(3)+vz(5)+vz(9))/d
	z(12)=(vz(3)+vz(5)+vz(11))/d
	z(13)=(vz(3)+vz(7)+vz(11))/d
	z(14)=(vz(4)+vz(6)+vz(10))/d
	z(15)=(vz(4)+vz(6)+vz(12))/d
	z(16)=(vz(4)+vz(7)+vz(12))/d
	z(17)=(vz(5)+vz(8)+vz(11))/d
	z(18)=(vz(6)+vz(8)+vz(12))/d
	z(19)=(vz(7)+vz(11)+vz(12))/d
	z(20)=(vz(8)+vz(11)+vz(12))/d

	return
	end
	
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine readradialabun(input,grid,y,n,m)
	IMPLICIT NONE
	integer i,j,n,m,input
	real*8 grid(n),y(n,m),x0,y0(m),x1,y1(m)
	i=1
	read(input,*,end=102) x0,y0(1:m)
103	if(x0.ge.grid(i)) then
		y(i,1:m)=y0(1:m)
		i=i+1
		goto 103
	endif
100	read(input,*,end=102) x1,y1(1:m)
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		y(i,1:m)=y1(1:m)+(grid(i)-x1)*(y0(1:m)-y1(1:m))/(x0-x1)
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	do j=i,n
		y(j,1:m)=y(i-1,1:m)*grid(i-1)/grid(j)
	enddo
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine sort(x,n)
	IMPLICIT NONE
	integer n,i,j,imin
	real*8 x(n),min
	
	do j=1,n-1
	min=x(j)
	imin=j
	do i=j,n
		if(x(i).lt.min) then
			min=x(i)
			imin=i
		endif
	enddo
	min=x(j)
	x(j)=x(imin)
	x(imin)=min
	enddo
	
	return
	end

c-----------------------------------------------------------------------
c Numerical recipes subroutine to solve an ordinary differential
c equation. Used for the disk structure.
c-----------------------------------------------------------------------
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,derivs)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      real*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.d-150)
      INTEGER i,kmax,kount,nstp
      real*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) return!pause
c     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
c      pause 'too many steps in odeint'
      return
      END


      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      real*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER i
      real*8 errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     *PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        xnew=x+h
c        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END



      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      real*8 h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs
      INTEGER i
      real*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END


      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END



      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
c     changed by Allard Jan van Marle
c      if (h.eq.0.) pause 'bad xa input in splint'
c
      if (h.eq.0.) then 
         print*, 'bad xa input in splint'
	 read(*,*)
      endif
c
c
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END


      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(n)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END


	subroutine intext(x,y,n,x0,y0)
	IMPLICIT NONE
	integer n,i
	real*8 x(n),y(n),x0,y0
	real*8 x1,x2,y1,y2
	
	if(x0.lt.x(1)) then
		y0=y(1)
		return
	endif
	if(x0.gt.x(n)) then
		y0=10d0**(log10(y(n-1))+log10(y(n)/y(n-1))*log10(x0/x(n-1))/log10(x(n)/x(n-1)))
		return
	endif
	do i=1,n-1
		if(x(i).lt.x0.and.x(i+1).ge.x0) then
			y0=y(i)+(y(i+1)-y(i))*(x0-x(i))/(x(i+1)-x(i))
			return
		endif
	enddo
	
	y0=y(n)

	return
	end


	subroutine SimpleSmooth(x,y,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 x(n),y(n),yy(n),ya
	
	do j=1,n
		yy(1)=y(1)
		do i=2,n-1
			ya=y(i-1)+(y(i+1)-y(i-1))*(x(i)-x(i-1))/(x(i+1)-x(i-1))
			yy(i)=(y(i)+ya)/2d0
		enddo
		yy(n)=y(n)
		y=yy
	enddo
	
	return
	end


	module GaussFit
		integer nx
		real*8,allocatable :: x(:),y(:)
	end module GaussFit

	subroutine FitGauss(x0,y0,nx0,sigma)
	use GaussFit
	IMPLICIT NONE
	integer ndim,iter,nx0,i,np,mp
	real*8 sigma,x0(nx0),y0(nx0)
	real*8 p(3,2),chi(3),ftol,Gauss
	external Gauss
	nx=nx0
	allocate(x(nx))
	allocate(y(nx))
	x(1:nx)=x0(1:nx0)
	y(1:nx)=y0(1:nx0)

	ndim=2
	mp=3
	np=2
	p(3,1)=-x(nx)/real(nx)
	p(1,2)=x(nx)/2d0
	p(2,1)=0d0
	p(2,2)=x(nx)/real(nx)
	p(3,1)=x(nx)/real(nx)
	p(3,2)=x(nx)/2d0

	do i=1,3
		chi(i)=Gauss(ndim,p(i,1:2))
	enddo
	ftol=1d-5

	call amoeba(p,chi,mp,np,ndim,ftol,Gauss,iter)

	sigma=p(1,2)
	
	deallocate(x)
	deallocate(y)
	return
	end
	


	real*8 function Gauss(ndim,var)
	use GaussFit
	IMPLICIT NONE
	integer ndim
	integer i
	real*8 G(nx),sum1,sum2,scale,var(ndim)
	
	do i=1,nx
		G(i)=exp(-0.5d0*((x(i)-var(1))/var(2))**2)
	enddo
	sum1=0d0
	sum2=0d0
	do i=1,nx
		sum1=sum1+G(i)*y(i)
		sum2=sum2+G(i)**2
	enddo
	scale=sum1/sum2
	Gauss=0d0
	do i=1,nx
		Gauss=Gauss+(y(i)-G(i)*scale)**2
	enddo

	return
	end
	
	
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      REAL*8 ftol,p(mp,np),y(mp),funk,TINY
      PARAMETER (NMAX=20,ITMAX=10000,TINY=1.e-10)
      EXTERNAL funk
CU    USES amotry,funk
      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL*8 rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      iter=0
1     do 12 n=1,ndim
        sum=0.
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
      if (iter.ge.ITMAX) then
      	print*,'ITMAX exceeded in amoeba'
		return
	  endif
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1d0)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2d0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(ndim,psum)
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END
      
      
      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      INTEGER ihi,mp,ndim,np,NMAX
      REAL*8 amotry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk
CU    USES funk
      INTEGER j
      REAL*8 fac1,fac2,ytry,ptry(NMAX)
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ndim,ptry)
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END



	subroutine writefitsfile(filename,im,n)
	IMPLICIT NONE
	character*500 filename
	integer n
	complex*16 im(n,n)
	real*8 array(n,n)

      integer status,unit,blocksize,bitpix,naxis,naxes(2)
      integer i,j,group,fpixel,nelements
      logical simple,extend,truefalse

	inquire(file=filename,exist=truefalse)
	if(truefalse) then
		write(*,'("FITS file already exists, overwriting")')
		write(9,'("FITS file already exists, overwriting")')
		open(unit=90,file=filename)
		close(unit=90,status='delete')
	endif

      status=0
C     Get an unused Logical Unit Number to use to create the FITS file
      call ftgiou(unit,status)
C     create the new empty FITS file
      blocksize=1
      call ftinit(unit,filename,blocksize,status)

C     initialize parameters about the FITS image (IMDIM x IMDIM 64-bit reals)
      simple=.true.
      bitpix=-64
      naxis=2
      naxes(1)=n
      naxes(2)=n
      extend=.true.

C     write the required header keywords
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

		array=real(im)

C     write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)
      call ftpprd(unit,group,fpixel,nelements,array,status)

C     close the file and free the unit number
      call ftclos(unit, status)
      call ftfiou(unit, status)


	return
	end


c-----------------------------------------------------------------------
c  Linear interpolation of value x between x1 and x2
c-----------------------------------------------------------------------

      real*8 function lininterpol(x,x1,x2,y1,y2)
      IMPLICIT NONE
      real*8 x1,x2,y1,y2,x,rc
      rc=(y2-y1)/(x2-x1)
      lininterpol=y1+(x-x1)*rc

      return
      end

c-----------------------------------------------------------------------
c  Logarithmic interpolation of value x between x1 and x2
c-----------------------------------------------------------------------

	real*8 function loginterpol(x,x1,x2,y1,y2)
	IMPLICIT NONE
	real*8 x1,x2,y1,y2,x,pow
	pow=log10(y2/y1)/log10(x2/x1)
	loginterpol=y1*(x/x1)**pow
	return
	end

c-----------------------------------------------------------------------
c  Linear/Log interpolation of value x in array arr of dimension n
c  log=0 is linear, log=1 is log
c-----------------------------------------------------------------------

	real*8 function arrinterpol(x,xarr,yarr,n,log)
	IMPLICIT NONE
	integer i,ix,n,log
	real*8 x,xarr(1:n),yarr(1:n)
	real*8 lininterpol,loginterpol
	
	! x smaller than array
	if (x.le.xarr(1)) then
	   arrinterpol=xarr(1)

	! x larger than array
	else if (x.ge.xarr(n)) then
	   arrinterpol=xarr(n)

	! x inside array
	else

	   do i=2,n
	      if (x.le.xarr(i)) ix=i 
	   enddo

	   if (log.eq.0) then
	      arrinterpol=lininterpol(x,xarr(ix-1),xarr(ix),yarr(ix-1),yarr(ix))
	   else
	      arrinterpol=loginterpol(x,xarr(ix-1),xarr(ix),yarr(ix-1),yarr(ix))
	   endif
	endif

	return
	end

c-----------------------------------------------------------------------
c This routine sets an isothermal density structure for a given
c temperature. Currently only works for power-law density structures.
c-----------------------------------------------------------------------

	subroutine SetIsoTdensity()
	use Parameters
	implicit none
	integer i,j,jj,njj
	real*8 sh(1:D%nR-1),Temp(1:D%nR-1,1:D%nTheta-1)
	real*8 theta,r,z,scale,f1,f2,sh1AU
	real*8 Mtot,Mtot0
	real*8 shgas,shdust ! test
	real*8 arrinterpol ! function

	call GetTemp(Temp(1:D%nR-1,1:D%nTheta-1))
	call IsoTsh(Temp(1:D%nR-1,D%nTheta-1),sh(1:D%nR-1))

	! sh at 1 AU	
c	print*,'Rin: ',D%Rin 
c	print*,'Mtot: ',D%Mtot / Msun
c	stop

	do i=1,D%nR-1
	   ! mass before
	   Mtot0=0d0
	   do j=1,D%nTheta-1
	      Mtot0=Mtot0+C(i,j)%dens*C(i,j)%V
	   enddo

	   do j=1,D%nTheta-1
	      njj=10  ! 10 calc per grid cell
	      C(i,j)%dens=0d0
	      do jj=1,njj
		 theta=D%thet(j)+(D%thet(j+1)-D%thet(j))*real(jj)/real(njj+1)
		 r=D%R_av(i)*sin(theta)
		 z=D%R_av(i)*cos(theta)
		 f1=D%Mtot/(2d0*pi*(D%Rout-D%Rin)*AU)* r**(-D%denspow)
		 f2=exp(-0.5*(z/sh(i))**2) / (sh(i)*sqrt(2d0*pi))
		 C(i,j)%dens=C(i,j)%dens+f1*f2/real(njj)
	      enddo
	   enddo		! j

	   ! mass after and rescale
	   Mtot=0d0
	   do j=1,D%nTheta-1
	      Mtot=Mtot+C(i,j)%dens*C(i,j)%V
	   enddo
	   do j=1,D%nTheta-1
	      C(i,j)%dens=C(i,j)%dens*Mtot0/Mtot
	      C(i,j)%gasdens=C(i,j)%dens ! for Deadzone
	   enddo
!	   print*,i,'scale:',Mtot0/Mtot

	   ! test if this went well
	   call calcscaleheight(i,shgas,shdust)
c	   print*,i,'sh correct to:',abs(1d0-sh(i)/AU / shgas) ! 7%

	enddo ! i

	end subroutine

c-----------------------------------------------------------------------
c This routine returns the scaleheight for a given dust temperature,
c in cgs (cm)
c-----------------------------------------------------------------------

	subroutine IsoTsh(Temp,sh)
	use Parameters
	implicit none
	integer i
	real*8 Temp(1:D%nR-1)
	real*8 sh(1:D%nR-1)
	doubleprecision, PARAMETER :: GG=6.6720000d-08
	doubleprecision, PARAMETER :: mu=2.3d0,amu=1.66053886d-24
	
	! H=c_s / Omega_k = sqrt(K T / mu m) / sqrt(G M_* / r^3)
	do i=1,D%nR-1
	   sh(i)= sqrt(kb*Temp(i)*D%R_av(i)**3d0 / (mu*amu*GG*D%Mstar))
	enddo

	end subroutine


