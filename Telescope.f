	subroutine Observe(tel)
	use NAG
	use Parameters
	IMPLICIT NONE
	type(Telescope) tel
	type(RPhiImage) image
	real*8 spec(nlam),scatspec(nlam),flux,scatflux,R0,angle,FWHM1,FWHM2
	real*8 arrinterpol,fstar1,fstar2
	real*8 V(100),starttime,stoptime,tottime,fluxQ,specQ(nlam),clight,Resolution
	real*8 basegrid(100),V2(100,100),phase(100),phase2(100,100)	!Gijsexp : nbase=100
	integer nbase		!Gijsexp
	integer i,j,k,i_up,i_low
	character*500 specfile
	real*8,allocatable :: velo(:),velo_flux(:),velo_flux_R(:)
	character*10 poptype
	parameter(clight=2.9979d5) !km/s
	real*8 Reddening,compute_dlam,ExtISM,temp,tot
	real*8 radtau,radialtau,tau1,Av
	integer ntau1
	integer nlam_obs
	real*8,allocatable :: lam_obs(:)
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	call cpu_time(starttime)

	angle=pi*tel%angle/180d0

	if(.not.tracescat) then
		tel%Nphot=0
		tel%NphotAngle=0
	endif

	if(tel%nlam_obs.lt.0) then
		nlam_obs=nlam
		allocate(lam_obs(nlam_obs))
		lam_obs=lam
	else
		nlam_obs=tel%nlam_obs
		if(nlam_obs.gt.nlam) nlam_obs=nlam
		allocate(lam_obs(nlam_obs))
		do i=1,nlam_obs
			lam_obs(i)=10d0**(log10(tel%lam1)+log10(tel%lam2/tel%lam1)*real(i-1)/real(nlam_obs-1))
		enddo
	endif

	outfluxcontr=tel%fluxcontr

	fastobs=tel%fastobs
	if(scat_how.eq.2) fastobs=.false.

	do i=1,ngrains
		Grain(i)%trace=tel%trace(i)
	enddo
	tracestar=tel%tracestar
	traceemis=tel%traceemis
	tracescat=tel%tracescat
	tracegas=tel%tracegas

	makeangledependence=.false.

	Av=D%Av
	if(adjustAv) then
		radtau=radialtau(0.55d0,tau1,ntau1,1)
		Av=Av+2.5*log10(exp(-radtau))
		if(Av.lt.0d0) Av=0d0
		write(*,'("Av adjusted to: ",f5.2)') Av
	endif

	if(scattering.and.(tel%Nphot+tel%NphotAngle).ne.0.and..not.storescatt) then
		do i=0,D%nR
		do j=1,D%nTheta-1
			if(scat_how.ne.2) then
				C(i,j)%nrg=1
			else
				if(int((D%thet(j+1)-D%thet(j))*270d0/pi).gt.1) then
					C(i,j)%nrg=int((D%thet(j+1)-D%thet(j))*270d0/pi)
				else
					C(i,j)%nrg=1
				endif
			endif
			allocate(C(i,j)%thetrg(C(i,j)%nrg+1))
			do k=1,C(i,j)%nrg+1
				C(i,j)%thetrg(k)=cos(D%thet(j)+real(k-1)*(D%thet(j+1)-D%thet(j))/real(C(i,j)%nrg))**2
			enddo
			allocate(C(i,j)%scattfield(C(i,j)%nrg,0:NPHISCATT/2,2))
			if(scat_how.eq.2) then
				allocate(C(i,j)%scattQ(C(i,j)%nrg,0:NPHISCATT/2,2))
				allocate(C(i,j)%scattU(C(i,j)%nrg,0:NPHISCATT/2,2))
				allocate(C(i,j)%scattV(C(i,j)%nrg,0:NPHISCATT/2,2))
			endif
		enddo
		enddo
	else
		do i=0,D%nR
		do j=1,D%nTheta-1
			C(i,j)%nrg=1
			allocate(C(i,j)%thetrg(C(i,j)%nrg+1))
			do k=1,C(i,j)%nrg+1
				C(i,j)%thetrg(k)=cos(D%thet(j)+real(k-1)*(D%thet(j+1)-D%thet(j))/real(C(i,j)%nrg))**2
			enddo
		enddo
		enddo
	endif
	
	image%scaletype=tel%scaletype
	if(tel%scaletype.eq.1) then
		tel%width=tel%width/(2d0*sqrt(-2d0*log(0.5d0)))*parsec/D%Distance
	else if(tel%scaletype.eq.2) then
		tel%width=tel%width/(2d0*sqrt(-2d0*log(0.5d0)))
	endif
	
	useobspol=tel%usepol
	nexits=tel%nexits
	if(tel%kind(1:8).eq.'SPECTRUM') then
		readmcscat=tel%readmcscat
		call TracePath(image,angle,tel%nphi,tel%nr,0.55d0)
		write(specfile,'(a,"spectrum",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		open(unit=30,file=specfile,RECL=6000)
		call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		ExtISM=Reddening(tel%lam1,compute_dlam(tel%lam1),Av)
		fstar1=arrinterpol(tel%lam1,lam(nlam),D%Fstar,nlam,1)
		fstar2=arrinterpol(tel%lam2,lam(nlam),D%Fstar,nlam,1)
		if(scat_how.ne.2) then
		   write(30,*) tel%lam1,1d23*flux*ExtISM/D%distance**2,
     &                         1d23*scatflux*ExtISM/D%distance**2,1d23*fstar1/D%distance**2
		else
		   write(30,*) tel%lam1,1d23*flux*ExtISM/D%distance**2,
     &                         1d23*scatflux*ExtISM/D%distance**2,1d23*fstar1/D%distance**2,
     &                         1d23*fluxQ*ExtISM/D%distance**2
		endif
		do j=1,nlam
		   if(lam(j).gt.tel%lam1.and.lam(j).lt.tel%lam2) then
		      call TraceFlux(image,lam(j),spec(j),scatspec(j),specQ(j),tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		      ExtISM=Reddening(lam(j),compute_dlam(lam(j)),Av)
		      if(scat_how.ne.2) then
			 write(30,*) lam(j),1d23*spec(j)*ExtISM/D%distance**2,
     &                               1d23*scatspec(j)*ExtISM/D%distance**2,1d23*D%Fstar(j)/D%distance**2
		      else
			 write(30,*) lam(j),1d23*spec(j)*ExtISM/D%distance**2,
     &                               1d23*scatspec(j)*ExtISM/D%distance**2,1d23*D%Fstar(j)/D%distance**2,
     &                               1d23*specQ(j)*ExtISM/D%distance**2
		      endif
		      call flush(30)
		   endif
		enddo
		call TraceFlux(image,tel%lam2,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		ExtISM=Reddening(tel%lam2,compute_dlam(tel%lam2),Av)
		if(scat_how.ne.2) then
		   write(30,*) tel%lam2,1d23*flux*ExtISM/D%distance**2,
     &          	       1d23*scatflux*ExtISM/D%distance**2,1d23*fstar2/D%distance**2
		else
		   write(30,*) tel%lam2,1d23*flux*ExtISM/D%distance**2,
     &                         1d23*scatflux*ExtISM/D%distance**2,1d23*fstar2/D%distance**2,
     &                         1d23*fluxQ*ExtISM/D%distance**2
		endif
		close(unit=30)
	else if(tel%kind(1:5).eq.'IMAGE') then
		readmcscat=.false.
		call TracePath(image,angle,tel%nphi,tel%nr,tel%lam1)
		if(tel%opening.ne.0d0) call opendisk(image,tel%opening)
		call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		ExtISM=Reddening(tel%lam1,compute_dlam(tel%lam1),Av)
		do i=1,tel%nfov
			if(tel%scaletype.eq.1) then
				image%rscale=1d0
				image%zscale=1d0
			else if(tel%scaletype.eq.2) then
				tel%fov(i)=tel%fov(i)*D%distance/parsec
				image%rscale=parsec/D%distance
				image%zscale=ExtISM*1d3*((tel%npixel/tel%fov(i))*D%distance/parsec)**2
			endif
			if(i.eq.1) call RImage(image,tel)
			call MakeImage(image,tel,tel%fov(i)/2d0)
		enddo
	else if(tel%kind(1:5).eq.'IFU') then
		readmcscat=tel%readmcscat
		call TracePath(image,angle,tel%nphi,tel%nr,0.55d0)
		call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		ExtISM=Reddening(tel%lam1,compute_dlam(tel%lam1),Av)
		do i=1,tel%nfov
			if(tel%scaletype.eq.1) then
				image%rscale=1d0
				image%zscale=1d0
			else if(tel%scaletype.eq.2) then
				tel%fov(i)=tel%fov(i)*D%distance/parsec
				image%rscale=parsec/D%distance
				image%zscale=ExtISM*1d3*((tel%npixel/tel%fov(i))*D%distance/parsec)**2
			endif
			call MakeImage(image,tel,tel%fov(i)/2d0)
		enddo
		do j=1,nlam
			if(lam(j).gt.tel%lam1.and.lam(j).lt.tel%lam2) then
				call TraceFlux(image,lam(j),flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
				do i=1,tel%nfov
					if(tel%scaletype.eq.1) then
						image%rscale=1d0
						image%zscale=1d0
					else if(tel%scaletype.eq.2) then
						image%rscale=parsec/D%distance
						image%zscale=ExtISM*1d3*((tel%npixel/tel%fov(i))*D%distance/parsec)**2
					endif
					call MakeImage(image,tel,tel%fov(i)/2d0)
				enddo
			endif
		enddo
		call TraceFlux(image,tel%lam2,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		ExtISM=Reddening(tel%lam2,compute_dlam(tel%lam2),Av)
		do i=1,tel%nfov
			if(tel%scaletype.eq.1) then
				image%rscale=1d0
				image%zscale=1d0
			else if(tel%scaletype.eq.2) then
				image%rscale=parsec/D%distance
				image%zscale=ExtISM*1d3*((tel%npixel/tel%fov(i))*D%distance/parsec)**2
			endif
			call MakeImage(image,tel,tel%fov(i)/2d0)
		enddo
	else if(tel%kind(1:10).eq.'VISIBILITY') then
		readmcscat=.false.
		call TracePath(image,angle,tel%nphi,tel%nr,0.55d0)
		write(specfile,'(a,"visibility",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		open(unit=30,file=specfile,RECL=6000)
		write(30,'("# column   1: wavelength")')
		write(30,'("# column   2: full disk")')
		do k=1,tel%nbaseline
			write(30,'("# column ",i3," and ",i3,": vis and phase for: baseline ",f10.3,", angle ",f10.3)') 
     1                  k+2,k+2+tel%nbaseline,tel%b(k),180d0*tel%theta(k)/pi
		enddo

		call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		ExtISM=Reddening(tel%lam1,compute_dlam(tel%lam1),Av)
		do k=1,tel%nbaseline
			call Visibility(image,tel%b(k),tel%theta(k),tel%lam1,V(k),phase(k))
		enddo
		write(30,*) tel%lam1,1d23*flux*ExtISM/D%distance**2,V(1:tel%nbaseline),phase(1:tel%nbaseline)
		do j=1,nlam_obs
			if(lam_obs(j).gt.tel%lam1.and.lam_obs(j).lt.tel%lam2) then
				call TraceFlux(image,lam_obs(j),spec(j),scatspec(j),specQ(j),tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
				do k=1,tel%nbaseline
					call Visibility(image,tel%b(k),tel%theta(k),lam_obs(j),V(k),phase(k))
				enddo
				write(30,*) lam_obs(j),1d23*spec(j)*ExtISM/D%distance**2,V(1:tel%nbaseline),phase(1:tel%nbaseline)
				call flush(30)
			endif
		enddo
		call TraceFlux(image,tel%lam2,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		ExtISM=Reddening(tel%lam2,compute_dlam(tel%lam2),Av)
		do k=1,tel%nbaseline
			call Visibility(image,tel%b(k),tel%theta(k),tel%lam2,V(k),phase(k))
		enddo
		write(30,*) tel%lam2,1d23*flux*ExtISM/D%distance**2,V(1:tel%nbaseline),phase(1:tel%nbaseline)
		close(unit=30)
	else if(tel%kind(1:7).eq.'BASEVIS') then
c is still without interstellar extinction
c       Gijsexp: 
		readmcscat=.false.
		call TracePath(image,angle,tel%nphi,tel%nr,tel%lam1)
		write(specfile,'(a,"basevis",i1,f3.1,a,".dat")') 
     &                  outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		open(unit=30,file=specfile,RECL=6000)
		write(30,'("# column   1: baseline")')
		write(30,'("# column   2: full disk (only last wavelength)")')
		do j=1,tel%nlam
		   write(30,'("# column ",i3," and ",i3,": vis and phase for wavelength ",f10.3,", angle ",f10.3)') 
     &		j+2,j+2+tel%nlam,tel%lam(j),180d0*tel%theta(j)/pi
		enddo
c$$$		write(30,'("# column ",i3," and further: phase")') tel%nlam+3
		Call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		do j=1,tel%nlam
		   if (j.ne.1.and.tel%lam(j).ne.tel%lam(j-1)) then
		      call TraceFlux(image,tel%lam(j),flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		   endif

		   nbase=100
		   do k=1,nbase
c		   basegrid(k)=tel%b(1)*(tel%b(2)/tel%b(1))**((k-1d0)/(nbase-1d0))  ! log grid
		      basegrid(k)=tel%b(1)+(tel%b(2)-tel%b(1))*((k-1d0)/(nbase-1d0)) ! linear grid
		      call Visibility(image,basegrid(k),tel%theta(j),tel%lam(j),V2(j,k),phase2(j,k))
		   enddo
		enddo
		do k=1,nbase
		   write(30,*) basegrid(k),1d23*flux/D%distance**2,V2(1:tel%nlam,k),phase2(1:tel%nbaseline,k)
		enddo
		close(unit=30)
c       End add
	else if(tel%kind(1:4).eq.'FWHM') then
c is still without interstellar extinction
		readmcscat=.false.
		if(tel%scaletype.eq.1) then
			image%rscale=1d0
			image%zscale=1d0
		else if(tel%scaletype.eq.2) then
			image%rscale=parsec/D%distance
			image%zscale=1d3*(D%distance/parsec)**2
		endif
		call TracePath(image,angle,tel%nphi,tel%nr,tel%lam1)
		call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		call FWHM(image,tel%npixel,tel%nint,FWHM1,FWHM2,tel%width*(2d0*sqrt(-2d0*log(0.5d0)))/parsec*D%Distance)
	else if(tel%kind(1:8).eq.'SPECFWHM') then
c is still without interstellar extinction
		readmcscat=.false.
		if(tel%scaletype.eq.1) then
			image%rscale=1d0
			image%zscale=1d0
		else if(tel%scaletype.eq.2) then
			image%rscale=parsec/D%distance
			image%zscale=1d3*(D%distance/parsec)**2
		endif
		call TracePath(image,angle,tel%nphi,tel%nr,0.55d0)
		write(specfile,'(a,"FWHM",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		open(unit=30,file=specfile,RECL=6000)
		write(30,'("# column   1: wavelength")')
		write(30,'("# column   2: long axis")')
		write(30,'("# column   3: short axis")')

		call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		call FWHM(image,tel%npixel,tel%nint,FWHM1,FWHM2,tel%width*(2d0*sqrt(-2d0*log(0.5d0)))/parsec*D%Distance)
		write(30,*) tel%lam1,FWHM1*image%rscale,FWHM2*image%rscale,flux
		do j=1,nlam
			if(lam(j).gt.tel%lam1.and.lam(j).lt.tel%lam2) then
				call TraceFlux(image,lam(j),spec(j),scatspec(j),specQ(j),tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
				call FWHM(image,tel%npixel,tel%nint,FWHM1,FWHM2,tel%width*(2d0*sqrt(-2d0*log(0.5d0)))/parsec*D%Distance)
				write(30,*) lam(j),FWHM1*image%rscale,FWHM2*image%rscale,spec(j)
				call flush(30)
			endif
		enddo
		call TraceFlux(image,tel%lam2,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		call FWHM(image,tel%npixel,tel%nint,FWHM1,FWHM2,tel%width*(2d0*sqrt(-2d0*log(0.5d0)))/parsec*D%Distance)
		write(30,*) tel%lam2,FWHM1*image%rscale,FWHM2*image%rscale,flux
		close(unit=30)
	else if(tel%kind(1:4).eq.'TELFWHM') then
c is still without interstellar extinction
		readmcscat=.false.
c		tel%fov(1)=D%R(D%nR)*2d0
		if(tel%scaletype.eq.1) then
			image%rscale=1d0
			image%zscale=1d0
		else if(tel%scaletype.eq.2) then
			tel%fov(1)=tel%fov(1)*D%distance/parsec
			image%rscale=parsec/D%distance
			image%zscale=1d3*(D%distance/parsec)**2
		endif
		call TracePath(image,angle,tel%nphi,tel%nr,tel%lam1)
		call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		call TELFWHM(image,tel,tel%fov(1)/2d0,FWHM1,FWHM2)
	else if(tel%kind(1:11).eq.'TELSPECFWHM') then
c is still without interstellar extinction
		readmcscat=tel%readmcscat
		call TracePath(image,angle,tel%nphi,tel%nr,0.55d0)
		call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		write(specfile,'(a,"FWHM",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		open(unit=30,file=specfile,RECL=6000)
		write(30,'("# column   1: wavelength")')
		write(30,'("# column   2: East-West axis")')
		write(30,'("# column   3: North-South axis")')
c		tel%fov(1)=D%R(D%nR)*2d0
		if(tel%scaletype.eq.1) then
			image%rscale=1d0
			image%zscale=1d0
		else if(tel%scaletype.eq.2) then
			tel%fov(1)=tel%fov(1)*D%distance/parsec
			image%rscale=parsec/D%distance
			image%zscale=1d3*((tel%npixel/tel%fov(1))*D%distance/parsec)**2
		endif
		call TELFWHM(image,tel,tel%fov(1)/2d0,FWHM1,FWHM2)
		write(30,*) tel%lam1,FWHM1*image%rscale,FWHM2*image%rscale,flux
		call flush(30)
		do j=1,nlam
			if(lam(j).gt.tel%lam1.and.lam(j).lt.tel%lam2) then
				call TraceFlux(image,lam(j),flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
				call TELFWHM(image,tel,tel%fov(1)/2d0,FWHM1,FWHM2)
				write(30,*) lam(j),FWHM1*image%rscale,FWHM2*image%rscale,spec(j)
				call flush(30)
			endif
		enddo
		call TraceFlux(image,tel%lam2,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		call TELFWHM(image,tel,tel%fov(1)/2d0,FWHM1,FWHM2)
		write(30,*) tel%lam2,FWHM1*image%rscale,FWHM2*image%rscale,flux
		close(unit=30)
	else if(tel%kind(1:8).eq.'TAU1TEMP') then
c is still without interstellar extinction
		write(specfile,'(a,"tau1temp",i1,i1,i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%lam1)/1000d0)
     &			,int((tel%lam1-1000d0*int(tel%lam1)/1000d0)/100d0)
     &			,int((tel%lam1-100d0*int(tel%lam1/100d0))/10d0)
     &			,tel%lam1-10d0*int((tel%lam1/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		call tau1temp(specfile,tel%lam1)		
	else if(tel%kind(1:4).eq.'LINE') then
		if(tel%trans_nr2.lt.tel%trans_nr1) tel%trans_nr2=tel%trans_nr1
		allocate(velo(tel%nvelo))
		allocate(velo_flux(tel%nvelo))
		allocate(velo_flux_R(tel%nvelo))
		readmcscat=.false.
		write(specfile,'(a,"line",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		open(unit=30,file=specfile,RECL=6000)
		write(specfile,'(a,"specline",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		open(unit=31,file=specfile,RECL=6000)
		call TracePathLine(image,angle,tel%nphi,tel%nr,tel%nt,tel%nstar,tel%lam1)
		do j=tel%trans_nr1,tel%trans_nr2
			call TraceLine(image,j,velo_flux,tel%Nphot,tel%NphotAngle,tel%linefile,tel%popfile,velo,tel%dvelo,tel%nvelo,
     &								tel%abun,poptype,i_up,i_low)
			write(*,'("Writing transition number: ",i5," (",i3,"-",i3,")")') j,i_up,i_low
			write(9,'("Writing transition number: ",i5," (",i3,"-",i3,")")') j,i_up,i_low
			write(30,'("# transition nr ",i)') j
			write(30,'("# up, low       ",i,i)') i_up,i_low
			write(30,'("# population    ",a)') trim(poptype)
			ExtISM=Reddening(image%lam,compute_dlam(image%lam),Av)
c			Resolution=3000
c			do i=1,tel%nvelo
c				velo_flux_R(i)=0d0
c				tot=0d0
c				do k=1,tel%nvelo
c					temp=exp(-((velo(i)-velo(k))*Resolution/299792.458)**2)
c					tot=tot+temp
c					velo_flux_R(i)=velo_flux_R(i)+velo_flux(k)*temp
c				enddo
c				velo_flux_R(i)=velo_flux_R(i)/tot
c			enddo
			do i=1,tel%nvelo
				write(30,*) velo(i),1d23*velo_flux(i)*ExtISM/D%distance**2!,1d23*velo_flux_R(i)*ExtISM/D%distance**2
			enddo
			if(tel%trans_nr1.ne.tel%trans_nr2) write(30,*)
			do i=tel%nvelo,1,-1
				write(31,*) image%lam*sqrt((1d0+velo(i)/clight)/(1d0-velo(i)/clight)),1d23*velo_flux(i)*ExtISM/D%distance**2,
c     &											1d23*velo_flux_R(i)*ExtISM/D%distance**2,
     &											trim(poptype),i_up,i_low
			enddo
		enddo
		close(unit=30)
		close(unit=31)
		deallocate(velo)
		deallocate(velo_flux)
		deallocate(velo_flux_R)
	else
		stop
	endif
	
	if(allocated(image%p)) then
		do i=1,image%nr
			do j=1,image%nphi
				deallocate(image%p(i,j)%i)
				deallocate(image%p(i,j)%j)
				deallocate(image%p(i,j)%k)
				deallocate(image%p(i,j)%irg)
				deallocate(image%p(i,j)%v)
				deallocate(image%p(i,j)%phi1)
				deallocate(image%p(i,j)%phi2)
				deallocate(image%p(i,j)%jphi1)
				deallocate(image%p(i,j)%jphi2)
				deallocate(image%p(i,j)%rad)
				if(allocated(image%p(i,j)%velo1)) deallocate(image%p(i,j)%velo1)
				if(allocated(image%p(i,j)%velo2)) deallocate(image%p(i,j)%velo2)
			enddo
		enddo
		deallocate(image%p)
	endif
	if(allocated(image%R)) then
		deallocate(image%R)
		deallocate(image%Phi)
		deallocate(image%image)
		if(scat_how.eq.2) then
			deallocate(image%imageQ)
			deallocate(image%imageU)
			deallocate(image%imageV)
		endif
	endif

	if(scattering.and.(tel%Nphot+tel%NphotAngle).ne.0.and..not.storescatt) then
		do i=0,D%nR
		do j=1,D%nTheta-1
			deallocate(C(i,j)%scattfield)
			if(scat_how.eq.2) then
				deallocate(C(i,j)%scattQ)
				deallocate(C(i,j)%scattU)
				deallocate(C(i,j)%scattV)
			endif
		enddo
		enddo
	endif
	do i=0,D%nR
	do j=1,D%nTheta-1
		if(allocated(C(i,j)%thetrg)) deallocate(C(i,j)%thetrg)
	enddo
	enddo

	call cpu_time(stoptime)
	tottime=stoptime-starttime
	write(*,'("Observation time:",f11.4,"  s")') tottime
	write(9,'("Observation time:",f11.4,"  s")') tottime


	return
	end
	
	
	
	subroutine determineoutput(output,tel)
	use Parameters
	IMPLICIT NONE
	character*500 output,key,value
	character*500 line
	type(Telescope) tel(MAXOBS),def
	logical truefalse,setdef
	real*8 deflam
	real*8,allocatable :: temp(:)
	integer i,j

	inquire(file=output,exist=truefalse)
	if(.not.truefalse) then
		write(*,'("Observation file not found")')
		write(9,'("Observation file not found")')
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		stop
	endif

	write(*,'("Output observations:")')
	write(9,'("Output observations:")')

	open(unit=20,file=output,RECL=6000)

	setdef=.true.

	def%angle=D%IA
	deflam=(lam(1)+lam(nlam))/2d0
	def%lam1=lam(1)
	def%lam2=lam(nlam)
	def%angle1=1
	def%angle2=89
	def%nangle=45
	def%nbaseline=1
	allocate(def%b(1))
	allocate(def%theta(1))
	def%b(1)=100d0
	def%theta(1)=0d0
	def%nphi=45
	def%nr=2
	def%nt=1
	def%nstar=30
	def%Nphot=1000
	def%NphotAngle=1000
	def%npixel=500
	def%nint=1000
	def%nfov=1
	allocate(def%fov(1))
	def%fov(1)=D%R(D%nR)*2d0
	def%width=0d0
	def%snoise=0d0
	def%opening=0d0
	def%usepol=.true.
	def%readmcscat=.false.
	def%traceinverse=.false.
	def%nexits=5
	def%flag=''
	def%scaletype=1	! 1=Jy/pixel & AU, 2=mJy/arcsec & arcsec
	def%D=-1d0
	def%D2=-1d0
	def%spider=-1d0
	def%mask=1d0
	def%wmask=0.1d0
	def%texp=-1d0
	def%dlam=0.1d0
	def%RON=0d0
	def%psffile=' '
	def%Ptelescope=0d0
	def%APtelescope=0d0
	def%iwa=0d0
	def%iprad=0d0
	def%owa=0d0
	def%strehl=0d0
	def%fluxcontr=outfluxcontr
	allocate(def%trace(ngrains))
	do i=1,ngrains
		def%trace(i)=Grain(i)%trace
	enddo
	def%tracestar=tracestar
	def%traceemis=traceemis
	def%tracescat=tracescat
	def%tracegas=tracegas
	def%fastobs=fastobs
	def%nvelo=101
	def%dvelo=1d0
	def%abun=1d-4
	def%trans_nr1=1
	def%trans_nr2=1
	def%popfile=' '
	def%nlam_obs=-1
	
1	call ignorestar(20)
	read(20,'(a500)',end=2) line

	key=line(1:index(line,'=')-1)
	value=line(index(line,'=')+1:len_trim(line))

	if(value(1:1).eq.'"'.or.value(1:1).eq."'") then
		value=value(2:len_trim(value)-1)
	endif

	do i=1,len_trim(key)
		if(iachar(key(i:i)).ge.65.and.iachar(key(i:i)).le.90) then
			key(i:i)=achar(iachar(key(i:i))+32)
		endif
	enddo

	if(key.eq.'type') then
		setdef=.false.
		nobs=nobs+1
		if(nobs.gt.MAXOBS) then
			write(*,'("Maximum nr of observations reached (",i5,")")') MAXOBS
			write(9,'("Maximum nr of observations reached (",i5,")")') MAXOBS
			write(*,'("--------------------------------------------------------")')
			write(9,'("--------------------------------------------------------")')
			stop
		endif		
		read(value,*) tel(nobs)%kind
		write(*,'("Type: ",a)') tel(nobs)%kind(1:len_trim(tel(nobs)%kind))
		write(9,'("Type: ",a)') tel(nobs)%kind(1:len_trim(tel(nobs)%kind))
		allocate(tel(nobs)%fov(def%nfov))
		if(tel(nobs)%kind.eq.'SPECTRUM') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=def%lam1
			tel(nobs)%lam2=def%lam2
			tel(nobs)%readmcscat=def%readmcscat
		else if(tel(nobs)%kind.eq.'SPECANGLES') then
			tel(nobs)%lam1=def%lam1
			tel(nobs)%lam2=def%lam2
			tel(nobs)%angle1=def%angle1
			tel(nobs)%angle2=def%angle2
			tel(nobs)%nangle=def%nangle
			tel(nobs)%readmcscat=def%readmcscat
		else if(tel(nobs)%kind.eq.'IMAGE') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=deflam
			tel(nobs)%nfov=def%nfov
			tel(nobs)%fov=def%fov
			tel(nobs)%npixel=def%npixel
			tel(nobs)%nint=def%nint
			tel(nobs)%width=def%width
			tel(nobs)%D=def%D
			tel(nobs)%D2=def%D2
			tel(nobs)%spider=def%spider
			tel(nobs)%mask=def%mask
			tel(nobs)%wmask=def%wmask
			tel(nobs)%psffile=def%psffile
			tel(nobs)%opening=def%opening
			tel(nobs)%scaletype=def%scaletype			
			tel(nobs)%Ptelescope=def%Ptelescope
			tel(nobs)%APtelescope=def%APtelescope
			tel(nobs)%iwa=def%iwa
			tel(nobs)%owa=def%owa
			tel(nobs)%strehl=def%strehl
		else if(tel(nobs)%kind.eq.'IFU') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=def%lam1
			tel(nobs)%lam2=def%lam2
			tel(nobs)%nfov=def%nfov
			tel(nobs)%fov=def%fov
			tel(nobs)%npixel=def%npixel
			tel(nobs)%nint=def%nint
			tel(nobs)%width=def%width
			tel(nobs)%D=def%D
			tel(nobs)%D2=def%D2
			tel(nobs)%spider=def%spider
			tel(nobs)%mask=def%mask
			tel(nobs)%wmask=def%wmask
			tel(nobs)%psffile=def%psffile
			tel(nobs)%opening=def%opening
			tel(nobs)%scaletype=def%scaletype
			tel(nobs)%Ptelescope=def%Ptelescope
			tel(nobs)%APtelescope=def%APtelescope
			tel(nobs)%iwa=def%iwa
			tel(nobs)%owa=def%owa
			tel(nobs)%strehl=def%strehl
		else if(tel(nobs)%kind.eq.'VISIBILITY') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=def%lam1
			tel(nobs)%lam2=def%lam2
			tel(nobs)%nbaseline=def%nbaseline
			allocate(tel(nobs)%b(100))
			allocate(tel(nobs)%theta(100))
			tel(nobs)%b(1)=def%b(1)
			tel(nobs)%theta(1)=def%theta(1)
			tel(nobs)%theta(1)=pi*tel(nobs)%theta(1)/180d0
		else if(tel(nobs)%kind(1:7).eq.'BASEVIS') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=deflam
			tel(nobs)%nbaseline=def%nbaseline
			tel(nobs)%nlam=def%nlam
			allocate(tel(nobs)%b(100))
			allocate(tel(nobs)%theta(100))
			allocate(tel(nobs)%lam(100))
			tel(nobs)%b(1)=100d0
			tel(nobs)%theta(1)=0d0
			tel(nobs)%theta(1)=pi*tel(nobs)%theta(1)/180d0
		else if(tel(nobs)%kind.eq.'FWHM') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=def%lam1
			tel(nobs)%lam2=def%lam2
			tel(nobs)%npixel=def%npixel
			tel(nobs)%nint=def%nint
			tel(nobs)%width=def%width
			tel(nobs)%D=def%D
			tel(nobs)%D2=def%D2
			tel(nobs)%spider=def%spider
			tel(nobs)%mask=def%mask
			tel(nobs)%wmask=def%wmask
			tel(nobs)%psffile=def%psffile
			tel(nobs)%scaletype=def%scaletype
		else if(tel(nobs)%kind.eq.'SPECFWHM') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=def%lam1
			tel(nobs)%lam2=def%lam2
			tel(nobs)%npixel=def%npixel
			tel(nobs)%nint=def%nint
			tel(nobs)%width=def%width
			tel(nobs)%D=def%D
			tel(nobs)%D2=def%D2
			tel(nobs)%spider=def%spider
			tel(nobs)%mask=def%mask
			tel(nobs)%wmask=def%wmask
			tel(nobs)%psffile=def%psffile
			tel(nobs)%scaletype=def%scaletype
		else if(tel(nobs)%kind.eq.'TELFWHM') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=deflam
			tel(nobs)%nfov=def%nfov
			tel(nobs)%fov=def%fov
			tel(nobs)%npixel=def%npixel
			tel(nobs)%nint=def%nint
			tel(nobs)%width=def%width
			tel(nobs)%D=def%D
			tel(nobs)%D2=def%D2
			tel(nobs)%spider=def%spider
			tel(nobs)%mask=def%mask
			tel(nobs)%wmask=def%wmask
			tel(nobs)%psffile=def%psffile
			tel(nobs)%opening=def%opening
			tel(nobs)%scaletype=def%scaletype			
			tel(nobs)%Ptelescope=def%Ptelescope
			tel(nobs)%APtelescope=def%APtelescope
			tel(nobs)%iwa=def%iwa
			tel(nobs)%owa=def%owa
			tel(nobs)%strehl=def%strehl
		else if(tel(nobs)%kind.eq.'TELSPECFWHM') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=def%lam1
			tel(nobs)%lam2=def%lam2
			tel(nobs)%nfov=def%nfov
			tel(nobs)%fov=def%fov
			tel(nobs)%npixel=def%npixel
			tel(nobs)%nint=def%nint
			tel(nobs)%width=def%width
			tel(nobs)%D=def%D
			tel(nobs)%D2=def%D2
			tel(nobs)%spider=def%spider
			tel(nobs)%mask=def%mask
			tel(nobs)%wmask=def%wmask
			tel(nobs)%psffile=def%psffile
			tel(nobs)%opening=def%opening
			tel(nobs)%scaletype=def%scaletype
			tel(nobs)%Ptelescope=def%Ptelescope
			tel(nobs)%APtelescope=def%APtelescope
			tel(nobs)%iwa=def%iwa
			tel(nobs)%owa=def%owa
			tel(nobs)%strehl=def%strehl
		else if(tel(nobs)%kind.eq.'TAU1TEMP') then
			tel(nobs)%lam1=deflam
		else if(tel(nobs)%kind.eq.'LINE') then
			tel(nobs)%angle=def%angle
			tel(nobs)%lam1=deflam
			tel(nobs)%nfov=def%nfov
			tel(nobs)%fov=def%fov
			tel(nobs)%npixel=def%npixel
			tel(nobs)%nint=def%nint
			tel(nobs)%width=def%width
			tel(nobs)%D=def%D
			tel(nobs)%D2=def%D2
			tel(nobs)%spider=def%spider
			tel(nobs)%mask=def%mask
			tel(nobs)%wmask=def%wmask
			tel(nobs)%psffile=def%psffile
			tel(nobs)%opening=def%opening
			tel(nobs)%scaletype=def%scaletype			
			tel(nobs)%Ptelescope=def%Ptelescope
			tel(nobs)%APtelescope=def%APtelescope
			tel(nobs)%iwa=def%iwa
			tel(nobs)%owa=def%owa
			tel(nobs)%strehl=def%strehl
			tel(nobs)%linefile=def%linefile
			tel(nobs)%popfile=def%popfile
			tel(nobs)%trans_nr1=def%trans_nr1
			tel(nobs)%trans_nr2=def%trans_nr2
			tel(nobs)%nvelo=def%nvelo
			tel(nobs)%dvelo=def%dvelo
			tel(nobs)%abun=def%abun
		endif
		tel(nobs)%nphi=def%nphi
		tel(nobs)%nr=def%nr
		tel(nobs)%nt=def%nt
		tel(nobs)%nstar=def%nstar
		tel(nobs)%Nphot=def%Nphot
		tel(nobs)%NphotAngle=def%NphotAngle
		tel(nobs)%usepol=def%usepol
		tel(nobs)%traceinverse=def%traceinverse
		tel(nobs)%nexits=def%nexits
		tel(nobs)%flag=def%flag
		tel(nobs)%dlam=def%dlam
		tel(nobs)%texp=def%texp
		tel(nobs)%RON=def%RON
		tel(nobs)%fluxcontr=def%fluxcontr
		tel(nobs)%snoise=def%snoise
		tel(nobs)%iprad=def%iprad
		allocate(tel(nobs)%trace(ngrains))
		do i=1,ngrains
			tel(nobs)%trace(i)=def%trace(i)
		enddo
		tel(nobs)%tracestar=def%tracestar
		tel(nobs)%traceemis=def%traceemis
		tel(nobs)%tracescat=def%tracescat
		tel(nobs)%tracegas=def%tracegas
		tel(nobs)%nlam_obs=def%nlam_obs
		tel(nobs)%fastobs=def%fastobs
	endif
	if(setdef) then
		if(key.eq.'nphi') read(value,*) def%nphi
		if(key.eq.'nrad') read(value,*) def%nr
		if(key.eq.'ntheta') read(value,*) def%nt
		if(key.eq.'nstar') read(value,*) def%nstar
		if(key.eq.'nphot') read(value,*) def%Nphot
		if(key.eq.'nphotstar') read(value,*) def%NphotAngle
		if(key.eq.'ninc') read(value,*) def%nangle
		if(key.eq.'lam') read(value,*) deflam
		if(key.eq.'lam1') read(value,*) def%lam1
		if(key.eq.'lam2') read(value,*) def%lam2
		if(key.eq.'nlam') read(value,*) def%nlam_obs
		if(key.eq.'inc') read(value,*) def%angle
		if(key.eq.'inc1') read(value,*) def%angle1
		if(key.eq.'inc2') read(value,*) def%angle2
		if(key.eq.'npix') read(value,*) def%npixel
		if(key.eq.'nint') read(value,*) def%nint
		if(key.eq.'d') read(value,*) def%D
		if(key.eq.'d2') read(value,*) def%D2
		if(key.eq.'spider') read(value,*) def%spider
		if(key.eq.'mask') read(value,*) def%mask
		if(key.eq.'wmask') read(value,*) def%wmask
		if(key.eq.'psffile') def%psffile=value
		if(key.eq.'ptel') read(value,*) def%Ptelescope
		if(key.eq.'aptel') read(value,*) def%APtelescope
		if(key.eq.'dlam') read(value,*) def%dlam
		if(key.eq.'texp') read(value,*) def%texp
		if(key.eq.'ron') read(value,*) def%RON
		if(key.eq.'width') read(value,*) def%width
		if(key.eq.'snoise') read(value,*) def%snoise
		if(key.eq.'scaletype') read(value,*) def%scaletype
		if(key.eq.'opening') read(value,*) def%opening
		if(key.eq.'usepol') read(value,*) def%usepol
		if(key.eq.'readmcscat') read(value,*) def%readmcscat
		if(key.eq.'inversetrace') read(value,*) def%traceinverse
		if(key.eq.'nexits') read(value,*) def%nexits
		if(key.eq.'flag') read(value,*) def%flag
		if(key.eq.'iwa') read(value,*) def%iwa
		if(key.eq.'owa') read(value,*) def%owa
		if(key.eq.'strehl') read(value,*) def%strehl
		if(key.eq.'iprad') read(value,*) def%iprad
		if(key.eq.'fluxcontr') read(value,*) def%fluxcontr
		if(key(1:3).eq.'fov') then
			if(len_trim(key).ne.3) then
				read(key(4:len_trim(key)),*) i
				if(i.gt.def%nfov) then
					allocate(temp(i))
					do j=1,def%nfov
						temp(j)=def%fov(j)
					enddo
					deallocate(def%fov)
					allocate(def%fov(i))
					do j=1,def%nfov
						def%fov(j)=temp(j)
					enddo
					read(value,*) def%fov(i)
					def%nfov=i
					deallocate(temp)
				else
					read(value,*) def%fov(i)
				endif
			else
				read(value,*) def%fov(1)
			endif			
		endif
		if(key.eq.'tracestar') read(value,*) def%tracestar
		if(key.eq.'traceemis') read(value,*) def%traceemis
		if(key.eq.'tracescat') read(value,*) def%tracescat
		if(key.eq.'tracegas') read(value,*) def%tracegas
		if(key(1:5).eq.'trace'.and.key.ne.'tracestar'.and.key.ne.'traceemis'.and.key.ne.'tracescat'.and.key.ne.'tracegas') then
			read(key(6:len_trim(key)),*) i
			if(i.le.100) read(value,*) def%trace(i)
		endif
		if(key.eq.'linefile') def%linefile=value
		if(key.eq.'popfile') def%popfile=value
		if(key.eq.'trans_nr') read(value,*) def%trans_nr1
		if(key.eq.'trans_nr1') read(value,*) def%trans_nr1
		if(key.eq.'trans_nr2') read(value,*) def%trans_nr2
		if(key.eq.'nvelo') read(value,*) def%nvelo
		if(key.eq.'dvelo') read(value,*) def%dvelo
		if(key.eq.'abun') read(value,*) def%abun
		if(key.eq.'fastobs') read(value,*) def%fastobs
	else
		if(key.eq.'nphi') read(value,*) tel(nobs)%nphi
		if(key.eq.'nrad') read(value,*) tel(nobs)%nr
		if(key.eq.'ntheta') read(value,*) tel(nobs)%nt
		if(key.eq.'nstar') read(value,*) tel(nobs)%nstar
		if(key.eq.'nphot') read(value,*) tel(nobs)%Nphot
		if(key.eq.'nphotstar') read(value,*) tel(nobs)%NphotAngle
		if(key.eq.'ninc') read(value,*) tel(nobs)%nangle
		if(key.eq.'lam') read(value,*) tel(nobs)%lam1
		if(key.eq.'lam1') read(value,*) tel(nobs)%lam1
		if(key.eq.'lam2') read(value,*) tel(nobs)%lam2
		if(key.eq.'nlam') read(value,*) tel(nobs)%nlam_obs
		if(key.eq.'inc') read(value,*) tel(nobs)%angle
		if(key.eq.'inc1') read(value,*) tel(nobs)%angle1
		if(key.eq.'inc2') read(value,*) tel(nobs)%angle2
		if(key.eq.'npix') read(value,*) tel(nobs)%npixel
		if(key.eq.'nint') read(value,*) tel(nobs)%nint
		if(key.eq.'d') read(value,*) tel(nobs)%D
		if(key.eq.'d2') read(value,*) tel(nobs)%D2
		if(key.eq.'spider') read(value,*) tel(nobs)%spider
		if(key.eq.'mask') read(value,*) tel(nobs)%mask
		if(key.eq.'wmask') read(value,*) tel(nobs)%wmask
		if(key.eq.'psffile') tel(nobs)%psffile=value
		if(key.eq.'ptel') read(value,*) tel(nobs)%Ptelescope
		if(key.eq.'aptel') read(value,*) tel(nobs)%APtelescope
		if(key.eq.'dlam') read(value,*) tel(nobs)%dlam
		if(key.eq.'texp') read(value,*) tel(nobs)%texp
		if(key.eq.'ron') read(value,*) tel(nobs)%RON
		if(key.eq.'width') read(value,*) tel(nobs)%width
		if(key.eq.'snoise') read(value,*) tel(nobs)%snoise
		if(key.eq.'scaletype') read(value,*) tel(nobs)%scaletype
		if(key.eq.'opening') read(value,*) tel(nobs)%opening
		if(key.eq.'usepol') read(value,*) tel(nobs)%usepol
		if(key.eq.'readmcscat') read(value,*) tel(nobs)%readmcscat
		if(key.eq.'inversetrace') read(value,*) tel(nobs)%traceinverse
		if(key.eq.'nexits') read(value,*) tel(nobs)%nexits
		if(key.eq.'flag') read(value,*) tel(nobs)%flag
		if(key.eq.'iwa') read(value,*) tel(nobs)%iwa
		if(key.eq.'owa') read(value,*) tel(nobs)%owa
		if(key.eq.'strehl') read(value,*) tel(nobs)%strehl
		if(key.eq.'iprad') read(value,*) tel(nobs)%iprad
		if(key.eq.'fluxcontr') read(value,*) tel(nobs)%fluxcontr
		if(key(1:4).eq.'base'.and.(tel(nobs)%kind.eq.'VISIBILITY'
     &		   .or.tel(nobs)%kind.eq.'BASEVIS')) then
			read(key(5:len_trim(key)),*) i
			if(i.gt.tel(nobs)%nbaseline) tel(nobs)%nbaseline=i 
			read(value,*) tel(nobs)%b(i)
		endif
		if(key(1:5).eq.'angle'.and.(tel(nobs)%kind.eq.'VISIBILITY'
     &		   .or.tel(nobs)%kind.eq.'BASEVIS')) then
			read(key(6:len_trim(key)),*) i
			if(i.gt.tel(nobs)%nbaseline) tel(nobs)%nbaseline=i 
			read(value,*) tel(nobs)%theta(i)
			tel(nobs)%theta(i)=pi*tel(nobs)%theta(i)/180d0
		endif
c       Gijsexp: allow more than two wavelength/angle combo's for basevis 
		if(key(1:3).eq.'lam'.and.tel(nobs)%kind.eq.'BASEVIS') then
			read(key(4:len_trim(key)),*) i
			if(i.gt.tel(nobs)%nlam) tel(nobs)%nlam=i 
			read(value,*) tel(nobs)%lam(i)
		endif
		if(key(1:3).eq.'fov') then
			if(len_trim(key).ne.3) then
				read(key(4:len_trim(key)),*) i
				if(i.gt.tel(nobs)%nfov) then
					allocate(temp(i))
					do j=1,tel(nobs)%nfov
						temp(j)=tel(nobs)%fov(j)
					enddo
					deallocate(tel(nobs)%fov)
					allocate(tel(nobs)%fov(i))
					do j=1,tel(nobs)%nfov
						tel(nobs)%fov(j)=temp(j)
					enddo
					read(value,*) tel(nobs)%fov(i)
					tel(nobs)%nfov=i
					deallocate(temp)
				else
					read(value,*) tel(nobs)%fov(i)
				endif
			else
				read(value,*) tel(nobs)%fov(1)
			endif			
		endif
		if(key.eq.'tracestar') read(value,*) tel(nobs)%tracestar
		if(key.eq.'traceemis') read(value,*) tel(nobs)%traceemis
		if(key.eq.'tracescat') read(value,*) tel(nobs)%tracescat
		if(key.eq.'tracegas') read(value,*) tel(nobs)%tracegas
		if(key(1:5).eq.'trace'.and.key.ne.'tracestar'.and.key.ne.'traceemis'.and.key.ne.'tracescat'.and.key.ne.'tracegas') then
			read(key(6:len_trim(key)),*) i
			if(i.le.100) read(value,*) tel(nobs)%trace(i)
		endif
		if(key.eq.'linefile') tel(nobs)%linefile=value
		if(key.eq.'popfile') tel(nobs)%popfile=value
		if(key.eq.'trans_nr') read(value,*) tel(nobs)%trans_nr1
		if(key.eq.'trans_nr1') read(value,*) tel(nobs)%trans_nr1
		if(key.eq.'trans_nr2') read(value,*) tel(nobs)%trans_nr2
		if(key.eq.'nvelo') read(value,*) tel(nobs)%nvelo
		if(key.eq.'dvelo') read(value,*) tel(nobs)%dvelo
		if(key.eq.'abun') read(value,*) tel(nobs)%abun
		if(key.eq.'fastobs') read(value,*) tel(nobs)%fastobs
	endif
	
	goto 1

2	close(unit=20)

	deallocate(def%b)
	deallocate(def%theta)
	
	if(2*(tel(nobs)%nvelo/2).eq.tel(nobs)%nvelo) tel(nobs)%nvelo=tel(nobs)%nvelo+1

	if(tel(nobs)%fastobs.and..not.fastobs) tel(nobs)%fastobs=.false.

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	
	return
	end
	



	


	subroutine opendisk(image,opening)
	use Parameters
	IMPLICIT NONE
	type(RPhiImage) image
	real*8 opening
	integer i,j
	
	do i=1,image%nr
		call tellertje(i,image%nr)
		do j=1,image%nphi
			call openpath(image%p(i,j),opening)
		enddo
	enddo
	
	return
	end
	
	subroutine openpath(p,opening)
	use Parameters
	type(path) p
	real*8 opening,dphi
	integer i
	
	do i=1,p%n
		dphi=abs(p%phi2(i)-p%phi1(i))
		if(abs(p%phi1(i)).le.opening.and.abs(p%phi2(i)).le.opening) then
			p%v(i)=0d0
			p%phi1(i)=opening*p%phi1(i)/abs(p%phi1(i))
			p%jphi1(i)=real(NPHISCATT)*p%phi1(i)/360d0
			p%jphi1(i)=p%jphi1(i)+1
			if(p%jphi1(i).lt.1) p%jphi1(i)=1-p%jphi1(i)
			if(p%jphi1(i).gt.NPHISCATT/2) p%jphi1(i)=NPHISCATT+1-p%jphi1(i)
			p%phi2(i)=opening*p%phi2(i)/abs(p%phi2(i))
			p%jphi2(i)=real(NPHISCATT)*p%phi2(i)/360d0
			p%jphi2(i)=p%jphi2(i)+1
			if(p%jphi2(i).lt.1) p%jphi2(i)=1-p%jphi2(i)
			if(p%jphi2(i).gt.NPHISCATT/2) p%jphi2(i)=NPHISCATT+1-p%jphi2(i)
		else if(dphi.ne.0d0) then
		if(abs(p%phi1(i)).le.opening) then
			if(p%phi1(i).ge.0d0) then
				p%phi1(i)=opening
			else
				p%phi1(i)=-opening
			endif
			p%jphi1(i)=real(NPHISCATT)*p%phi1(i)/360d0
			p%jphi1(i)=p%jphi1(i)+1
			if(p%jphi1(i).lt.1) p%jphi1(i)=1-p%jphi1(i)
			if(p%jphi1(i).gt.NPHISCATT/2) p%jphi1(i)=NPHISCATT+1-p%jphi1(i)
		endif
		if(abs(p%phi2(i)).le.opening) then
			if(p%phi2(i).ge.0d0) then
				p%phi2(i)=opening
			else
				p%phi2(i)=-opening
			endif
			p%jphi2(i)=real(NPHISCATT)*p%phi2(i)/360d0
			p%jphi2(i)=p%jphi2(i)+1
			if(p%jphi2(i).lt.1) p%jphi2(i)=1-p%jphi2(i)
			if(p%jphi2(i).gt.NPHISCATT/2) p%jphi2(i)=NPHISCATT+1-p%jphi2(i)
		endif
		p%v(i)=p%v(i)*abs(p%phi2(i)-p%phi1(i))/dphi
		endif
	enddo
	
	return
	end
	

		
c-----------------------------------------------------------------------
c This subroutine outputs the tau=1 temperature of the disk for
c radiation at wavelength lam0 to the file filename
c-----------------------------------------------------------------------
	subroutine tau1temp(filename,lam0)
	use Parameters
	IMPLICIT NONE
	character*500 filename
	real*8 lam0,tau1(D%nR)
	integer i,jtau1(D%nR)
	
	call calctau1R(jtau1,lam0)
	do i=1,D%nR-1
	   tau1(i)=C(i,jtau1(i))%T
	enddo

	open(unit=80,file=filename,RECL=6000)
	write(80,'("# temperature at the tau=1 surface, midplane, opt thin")')
	write(80,'("# wavelength: ",f7.1," micron")') lam0
	do i=1,D%nR-1
		write(80,*) D%R_av(i)/AU,(tau1(i)),C(i,D%nTheta-1)%T,C(i,1)%T
	enddo
	close(unit=80)
		
	return
	end

c-----------------------------------------------------------------------
c This subroutine returns the location of the tau=1 surface, at each
c radius, at wavelength lam0 
c-----------------------------------------------------------------------
	subroutine calctau1R(jtau1,lam0)
	use Parameters
	IMPLICIT NONE
	integer jtau1(D%nR)
	real*8 r(D%nR),ct,tau
	real*8 p,p0,p1,Kext,lam0
	logical escape,hitstar
	type(photon) phot
	integer i,j,l,n,ii,nl,iopac
	
	phot%nr=1
	tautot=0d0
	do i=1,D%nR-1
		r(i)=D%R_av(i)/AU
		phot%x=r(i)
		phot%y=0d0
		phot%z=sqrt(D%R(D%nR)**2-phot%x**2)
		ct=phot%z/D%R(D%nR)
		do j=1,D%nTheta-1
			if(ct.lt.D%Theta(j).and.ct.gt.D%Theta(j+1)) then
				phot%i=D%nR-1
				phot%j=j
			endif
		enddo
		phot%edgenr=2
		phot%onEdge=.true.
		phot%vx=0d0
		phot%vy=0d0
		phot%vz=-1d0
		tau=1d0
		phot%lam=lam0
		do j=1,nlam-1
			if(phot%lam.ge.lam(j).and.phot%lam.le.lam(j+1)) then
				phot%ilam1=j
				phot%ilam2=j+1
				phot%wl1=(lam(j+1)-phot%lam)/(lam(j+1)-lam(j))
				phot%wl2=(phot%lam-lam(j))/(lam(j+1)-lam(j))
			endif
		enddo
		phot%nu=nu(phot%ilam1)*phot%wl1+nu(phot%ilam2)*phot%wl2
		phot%E=1d0
		do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				Grain(ii)%KabsL(iopac)=Grain(ii)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kabs(iopac,phot%ilam2)*phot%wl2
					Grain(ii)%KextL(iopac)=Grain(ii)%Kext(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kext(iopac,phot%ilam2)*phot%wl2
			enddo
		enddo
		call trace2d(phot,tau,escape,hitstar,.false.)
		if (phot%z.ge.0d0) then
		   jtau1(i)=phot%j
		else
		   jtau1(i)=D%nTheta-1
		endif
	enddo
		
	return
	end
	

		
