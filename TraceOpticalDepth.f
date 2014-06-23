	subroutine TraceOpticalDepth(image,lam0)
	use Parameters
	type(RPhiImage) image
	real*8 lam0,wl1,wl2
	integer i,j,ilam1,ilam2,ip,jp

	image%lam=lam0
	if(lam0.le.lam(1)) then
		ilam1=1
		ilam2=2
		wl1=1d0
		wl2=0d0
	endif
	do i=1,nlam-1
		if(lam0.ge.lam(i).and.lam0.le.lam(i+1)) then
			ilam1=i
			ilam2=i+1
			wl1=(lam(i+1)-lam0)/(lam(i+1)-lam(i))
			wl2=(lam0-lam(i))/(lam(i+1)-lam(i))
		endif
	enddo
	if(lam0.ge.lam(nlam)) then
		ilam1=nlam-1
		ilam2=nlam
		wl1=0d0
		wl2=1d0
	endif

	do i=0,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%Kext=0d0
			do ii=1,ngrains
				do iopac=1,Grain(ii)%nopac
					C(i,j)%Kext=C(i,j)%Kext+(wl1*Grain(ii)%Kext(iopac,ilam1)
     &			+wl2*Grain(ii)%Kext(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
				enddo
			enddo
		enddo
	enddo

	do i=1,image%nr
	call tellertje(i,image%nr)
	do j=1,image%nphi
		image%image(i,j)=0d0
		do k=1,image%p(i,j)%n
			ip=image%p(i,j)%i(k)
			jp=image%p(i,j)%j(k)
			image%image(i,j)=image%image(i,j)+image%p(i,j)%v(k)*C(ip,jp)%dens*C(ip,jp)%Kext*AU
		enddo
	enddo
	enddo
	
	return
	end
	
