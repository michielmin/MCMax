	subroutine ReadParticleFits(input,p,Rayleigh,asym,asym2,wasym2,Pmax,iopac)
	use Parameters
	implicit none
	character*500 input
	character*80 comment,errmessage
	character*30 errtext
	integer*4 :: status,stat2,stat3,readwrite,unit,blocksize,nfound,group
	integer*4 :: firstpix,nbuffer,npixels,hdunum,hdutype,ix,iz,ilam
	integer*4 :: istat,stat4,tmp_int,stat5,stat6
	real*8  :: nullval,rho_gr,tmp
	logical*4 :: anynull
	integer*4, dimension(4) :: naxes

	real*8,allocatable :: array(:,:),matrix(:,:,:)

	type(particle) p,p0,p1
	integer i,j,ia,iopac,iread,nl_read
	logical truefalse,Rayleigh
	real*8 l0,l1,tot,tot2,theta,asym,Pmax,HG,asym2,wasym2

	p%qhp=.false.
	p%gascoupled=.false.

	allocate(p0%Kabs(1,nlam))
	allocate(p0%Ksca(1,nlam))
	allocate(p0%Kext(1,nlam))
	allocate(p0%F(1,nlam))
	allocate(p1%Kabs(1,nlam))
	allocate(p1%Ksca(1,nlam))
	allocate(p1%Kext(1,nlam))
	allocate(p1%F(1,nlam))

	! Get an unused Logical Unit Number to use to open the FITS file.
	status=0

	call ftgiou (unit,status)
	! Open file
	readwrite=0
	call ftopen(unit,input,readwrite,blocksize,status)
	if (status /= 0) then
	   write(*,'("Error reading particle file: ",a)') input(1:len_trim(input))
	   write(9,'("Error reading particle file: ",a)') input(1:len_trim(input))
	   write(*,'("--------------------------------------------------------")')
	   write(9,'("--------------------------------------------------------")')
	   stop
	endif
	group=1
	firstpix=1
	nullval=-999

	write(*,'("Reading particle file: ",a)') input(1:len_trim(input))
	write(9,'("Reading particle file: ",a)') input(1:len_trim(input))

	!------------------------------------------------------------------------
	! HDU0 : opacities
	!------------------------------------------------------------------------
	! Check dimensions
	call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

	nl_read = naxes(1)

	npixels=naxes(1)*naxes(2)

	! Read model info

	call ftgkye(unit,'a1',p%dust_moment1,comment,status)
	call ftgkye(unit,'a2',p%dust_moment2,comment,status)
	p%rv=sqrt(p%dust_moment2)*1d-4
	call ftgkye(unit,'a3',p%dust_moment3,comment,status)

c	call ftgkyj(unit,'mcfost2prodimo',mcfost(1)%mcfost2ProDiMo,comment,stat4)
 
	! read_image
	allocate(array(nl_read,4))

	call ftgpvd(unit,group,firstpix,npixels,nullval,array,anynull,status)


	!------------------------------------------------------------------------
	! HDU 1: matrix
	!------------------------------------------------------------------------
	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)	

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

	npixels=naxes(1)*naxes(2)*naxes(3)

	! read_image
	allocate(matrix(nl_read,6,180))

	call ftgpve(unit,group,firstpix,npixels,nullval,matrix,anynull,status)
   
				 
	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	!  Get the text string which describes the error
	if (status > 0) then
	   call ftgerr(status,errtext)
	   print *,'FITSIO Error Status =',status,': ',errtext

	   !  Read and print out all the error messages on the FITSIO stack
	   call ftgmsg(errmessage)
	   do while (errmessage .ne. ' ')
		  print *,errmessage
		  call ftgmsg(errmessage)
	   end do
	endif


	i=1
	iread=1
	l0=array(i,1)
	p0%Kext(1,1)=array(iread,2)
	p0%Kabs(1,1)=array(iread,3)
	p0%Ksca(1,1)=array(iread,4)
	do j=1,180
		p0%F(1,1)%F11(j)=matrix(iread,1,j)
		p0%F(1,1)%F12(j)=matrix(iread,2,j)
		p0%F(1,1)%F22(j)=matrix(iread,3,j)
		p0%F(1,1)%F33(j)=matrix(iread,4,j)
		p0%F(1,1)%F34(j)=matrix(iread,5,j)
		p0%F(1,1)%F44(j)=matrix(iread,6,j)
	enddo
	if(Rayleigh) then
	do ia=1,180
		theta=(real(ia)-0.5d0)*pi/180d0
		if(asym.lt.1d0.and.asym.gt.-1d0) then
			p0%F(1,1)%F11(ia)=HG(asym,theta)
			if(wasym2.ne.0d0) then
				p0%F(1,1)%F11(ia)=p0%F(1,1)%F11(ia)+wasym2*HG(asym2,theta)
			endif
			p0%F(1,1)%F12(ia)=-Pmax*p0%F(1,1)%F11(ia)*(1d0-cos(theta)**2)/(1d0+cos(theta)**2)
		else
			p0%F(1,1)%F11(ia)=(1d0+cos(theta)**2)/2d0
			p0%F(1,1)%F12(ia)=-(1d0-cos(theta)**2)/2d0
		endif
		p0%F(1,1)%F22(ia)=p0%F(1,1)%F11(ia)
		p0%F(1,1)%F33(ia)=cos(theta)*p0%F(1,1)%F11(ia)/(1d0+cos(theta)**2)
		p0%F(1,1)%F34(ia)=0d0
		p0%F(1,1)%F44(ia)=cos(theta)*p0%F(1,1)%F11(ia)/(1d0+cos(theta)**2)
	enddo
	endif
103	if(l0.ge.lam(i)) then
		p%Kext(iopac,i)=p0%Kext(1,1)
		p%Ksca(iopac,i)=p0%Ksca(1,1)
		p%Kabs(iopac,i)=p0%Kabs(1,1)
		p%F(iopac,i)=p0%F(1,1)
		call tellertje(i,nlam)
		i=i+1
		goto 103
	endif
100	iread=iread+1
	if(iread.gt.nl_read) goto 102
	l1=array(iread,1)
	p1%Kext(1,1)=array(iread,2)
	p1%Kabs(1,1)=array(iread,3)
	p1%Ksca(1,1)=array(iread,4)
	do j=1,180
		p1%F(1,1)%F11(j)=matrix(iread,1,j)
		p1%F(1,1)%F12(j)=matrix(iread,2,j)
		p1%F(1,1)%F22(j)=matrix(iread,3,j)
		p1%F(1,1)%F33(j)=matrix(iread,4,j)
		p1%F(1,1)%F34(j)=matrix(iread,5,j)
		p1%F(1,1)%F44(j)=matrix(iread,6,j)
	enddo
	if(Rayleigh) then
	do ia=1,180
		theta=(real(ia)-0.5d0)*pi/180d0
		if(asym.lt.1d0.and.asym.gt.-1d0) then
			p1%F(1,1)%F11(ia)=HG(asym,theta)
			if(wasym2.ne.0d0) then
				p1%F(1,1)%F11(ia)=p1%F(1,1)%F11(ia)+wasym2*HG(asym2,theta)
			endif
			p1%F(1,1)%F12(ia)=-Pmax*p1%F(1,1)%F11(ia)*(1d0-cos(theta)**2)/(1d0+cos(theta)**2)
		else
			p1%F(1,1)%F11(ia)=(1d0+cos(theta)**2)/2d0
			p1%F(1,1)%F12(ia)=-(1d0-cos(theta)**2)/2d0
		endif
		p1%F(1,1)%F22(ia)=p1%F(1,1)%F11(ia)
		p1%F(1,1)%F33(ia)=cos(theta)*p1%F(1,1)%F11(ia)/(1d0+cos(theta)**2)
		p1%F(1,1)%F34(ia)=0d0
		p1%F(1,1)%F44(ia)=cos(theta)*p1%F(1,1)%F11(ia)/(1d0+cos(theta)**2)
	enddo
	endif
101	if(lam(i).le.l1.and.lam(i).ge.l0) then
		p%Kext(iopac,i)=p1%Kext(1,1)+(lam(i)-l1)*(p0%Kext(1,1)-p1%Kext(1,1))/(l0-l1)
		p%Ksca(iopac,i)=p1%Ksca(1,1)+(lam(i)-l1)*(p0%Ksca(1,1)-p1%Ksca(1,1))/(l0-l1)
		p%Kabs(iopac,i)=p1%Kabs(1,1)+(lam(i)-l1)*(p0%Kabs(1,1)-p1%Kabs(1,1))/(l0-l1)
		p%F(iopac,i)%F11(1:180)=p1%F(1,1)%F11(1:180)+(lam(i)-l1)*(p0%F(1,1)%F11(1:180)-p1%F(1,1)%F11(1:180))/(l0-l1)
		p%F(iopac,i)%F12(1:180)=p1%F(1,1)%F12(1:180)+(lam(i)-l1)*(p0%F(1,1)%F12(1:180)-p1%F(1,1)%F12(1:180))/(l0-l1)
		p%F(iopac,i)%F22(1:180)=p1%F(1,1)%F22(1:180)+(lam(i)-l1)*(p0%F(1,1)%F22(1:180)-p1%F(1,1)%F22(1:180))/(l0-l1)
		p%F(iopac,i)%F33(1:180)=p1%F(1,1)%F33(1:180)+(lam(i)-l1)*(p0%F(1,1)%F33(1:180)-p1%F(1,1)%F33(1:180))/(l0-l1)
		p%F(iopac,i)%F34(1:180)=p1%F(1,1)%F34(1:180)+(lam(i)-l1)*(p0%F(1,1)%F34(1:180)-p1%F(1,1)%F34(1:180))/(l0-l1)
		p%F(iopac,i)%F44(1:180)=p1%F(1,1)%F44(1:180)+(lam(i)-l1)*(p0%F(1,1)%F44(1:180)-p1%F(1,1)%F44(1:180))/(l0-l1)
		call tellertje(i,nlam)
		i=i+1
		if(i.gt.nlam) goto 102
		goto 101
	endif
	l0=l1
	p0%Kext(1,1)=p1%Kext(1,1)
	p0%Ksca(1,1)=p1%Ksca(1,1)
	p0%Kabs(1,1)=p1%Kabs(1,1)
	p0%F(1,1)=p1%F(1,1)
	goto 100
102	continue
	do j=i,nlam
		call tellertje(j,nlam)
		p%Ksca(iopac,j)=p%Ksca(iopac,i-1)*(lam(i-1)/lam(j))**4
		p%Kabs(iopac,j)=p%Kabs(iopac,i-1)*(lam(i-1)/lam(j))**2
		p%Kext(iopac,j)=p%Kabs(iopac,j)+p%Ksca(iopac,j)
		p%F(iopac,j)=p%F(iopac,i-1)
	enddo


	if(nspike.gt.0.and.nspike.lt.180) then
c the nspike parameter removes the n degree spike in the forward direction.
		write(*,'("Making first ",i2," degrees isotropic")') nspike
		write(9,'("Making first ",i2," degrees isotropic")') nspike
	endif
	if(wasym2.ne.0d0.and..not.Rayleigh) then
		write(*,'("Adding HG fasefunction")')
		write(9,'("Adding HG fasefunction")')
	endif

	do j=1,nlam
	tot=0d0
	tot2=0d0
	do i=1,180
		tot=tot+p%F(iopac,j)%F11(i)*sin(pi*(real(i)-0.5)/180d0)
		tot2=tot2+sin(pi*(real(i)-0.5)/180d0)
	enddo
	do i=1,180
		p%F(iopac,j)%F11(i)=tot2*p%F(iopac,j)%F11(i)/tot
		p%F(iopac,j)%F12(i)=tot2*p%F(iopac,j)%F12(i)/tot
		p%F(iopac,j)%F22(i)=tot2*p%F(iopac,j)%F22(i)/tot
		p%F(iopac,j)%F33(i)=tot2*p%F(iopac,j)%F33(i)/tot
		p%F(iopac,j)%F34(i)=tot2*p%F(iopac,j)%F34(i)/tot
		p%F(iopac,j)%F44(i)=tot2*p%F(iopac,j)%F44(i)/tot
	enddo

	if(nspike.gt.0.and.nspike.lt.180) then
c the nspike parameter removes the n degree spike in the forward direction.
		do i=1,nspike
			p%F(iopac,j)%F12(i)=p%F(iopac,j)%F12(i)*p%F(iopac,j)%F11(nspike+1)/p%F(iopac,j)%F11(i)
			p%F(iopac,j)%F22(i)=p%F(iopac,j)%F22(i)*p%F(iopac,j)%F11(nspike+1)/p%F(iopac,j)%F11(i)
			p%F(iopac,j)%F33(i)=p%F(iopac,j)%F33(i)*p%F(iopac,j)%F11(nspike+1)/p%F(iopac,j)%F11(i)
			p%F(iopac,j)%F34(i)=p%F(iopac,j)%F34(i)*p%F(iopac,j)%F11(nspike+1)/p%F(iopac,j)%F11(i)
			p%F(iopac,j)%F44(i)=p%F(iopac,j)%F44(i)*p%F(iopac,j)%F11(nspike+1)/p%F(iopac,j)%F11(i)
			p%F(iopac,j)%F11(i)=p%F(iopac,j)%F11(nspike+1)
		enddo

		tot=0d0
		tot2=0d0
		do i=1,180
			tot=tot+p%F(iopac,j)%F11(i)*sin(pi*(real(i)-0.5)/180d0)
			tot2=tot2+sin(pi*(real(i)-0.5)/180d0)
		enddo
		print*,lam(j),tot/tot2
		p%Ksca(iopac,j)=p%Ksca(iopac,j)*tot/tot2
		p%Kext(iopac,j)=p%Kabs(iopac,j)+p%Ksca(iopac,j)
		do i=1,180
			p%F(iopac,j)%F11(i)=tot2*p%F(iopac,j)%F11(i)/tot
			p%F(iopac,j)%F12(i)=tot2*p%F(iopac,j)%F12(i)/tot
			p%F(iopac,j)%F22(i)=tot2*p%F(iopac,j)%F22(i)/tot
			p%F(iopac,j)%F33(i)=tot2*p%F(iopac,j)%F33(i)/tot
			p%F(iopac,j)%F34(i)=tot2*p%F(iopac,j)%F34(i)/tot
			p%F(iopac,j)%F44(i)=tot2*p%F(iopac,j)%F44(i)/tot
		enddo
	endif

	if(wasym2.ne.0d0.and..not.Rayleigh) then
		tot=0d0
		tot2=0d0
		do ia=1,180
			theta=(real(ia)-0.5d0)*pi/180d0
			p1%F(1,1)%F11(ia)=HG(asym2,theta)
			tot=tot+p1%F(1,1)%F11(i)*sin(theta)
			tot2=tot2+sin(theta)
		enddo
		do ia=1,180
			p%F(iopac,j)%F11(ia)=p%F(iopac,j)%F11(ia)*(1d0-wasym2)+wasym2*p1%F(1,1)%F11(ia)
		enddo
	endif

	enddo

	deallocate(p0%Kabs)
	deallocate(p0%Ksca)
	deallocate(p0%Kext)
	deallocate(p0%F)
	deallocate(p1%Kabs)
	deallocate(p1%Ksca)
	deallocate(p1%Kext)
	deallocate(p1%F)

	deallocate(array)
	deallocate(matrix)

	return

	end


