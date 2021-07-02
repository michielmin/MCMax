	subroutine Observe(tel)
	use NAG
	use Parameters
	IMPLICIT NONE
	type(Telescope) tel
	type(RPhiImage) image
	real*8 spec(nlam),scatspec(nlam),flux,scatflux,R0,angle,FWHM1,FWHM2
	real*8 arrinterpol,fstar1,fstar2
	real*8 starttime,stoptime,tottime,fluxQ,specQ(nlam),clight,Resolution
	real*8,allocatable :: V(:,:),basegrid(:),phase(:,:)
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
	vis_norm_fov=tel%vis_norm_fov

	tracestar=tel%tracestar
	traceemis=tel%traceemis
	tracescat=tel%tracescat
	tracegas=tel%tracegas

	if(.not.tracescat) then
		tel%Nphot=0
		tel%NphotAngle=0
	endif

	if(tel%nlam_obs.lt.0.and.tel%lamfile.eq.' ') then
		nlam_obs=nlam+2
		allocate(lam_obs(nlam_obs))
		nlam_obs=1
		lam_obs(1)=tel%lam1
		do i=1,nlam
			if(lam(i).gt.tel%lam1.and.lam(i).lt.tel%lam2) then
				nlam_obs=nlam_obs+1
				lam_obs(nlam_obs)=lam(i)
			endif
		enddo
		nlam_obs=nlam_obs+1
		lam_obs(nlam_obs)=tel%lam2
	else if(tel%lamfile.ne.' ') then
		open(unit=25,file=tel%lamfile)
		nlam_obs=0
1		read(25,*,end=2,err=1) tot
		nlam_obs=nlam_obs+1
		goto 1
2		continue
		allocate(lam_obs(nlam_obs))
		rewind(25)
		do i=1,nlam_obs
3			read(25,*,err=3) lam_obs(i)
		enddo
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
		if(tel%noangle) then
			write(specfile,'(a,"spectrum",a,".dat")') outdir(1:len_trim(outdir))
     &			,tel%flag(1:len_trim(tel%flag))
		else
			write(specfile,'(a,"spectrum",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		endif
		open(unit=30,file=specfile,RECL=6000)
		do j=1,nlam_obs
			call TraceFlux(image,lam_obs(j),spec(j),scatspec(j),specQ(j),tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
			ExtISM=Reddening(lam_obs(j),compute_dlam(lam_obs(j)),Av)
			fstar1=arrinterpol(lam_obs(j),lam,D%Fstar,nlam,1)
			if(scat_how.ne.2) then
				write(30,*) lam_obs(j),1d23*spec(j)*ExtISM/D%distance**2,
     &                               1d23*scatspec(j)*ExtISM/D%distance**2,1d23*fstar1/D%distance**2
			else
				write(30,*) lam_obs(j),1d23*spec(j)*ExtISM/D%distance**2,
     &                               1d23*scatspec(j)*ExtISM/D%distance**2,1d23*fstar1/D%distance**2,
     &                               1d23*specQ(j)*ExtISM/D%distance**2
			endif
			call flush(30)
		enddo
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
	else if(tel%kind(1:3).eq.'IFU') then
		readmcscat=tel%readmcscat
		call TracePath(image,angle,tel%nphi,tel%nr,0.55d0)
		do j=1,nlam_obs
			call TraceFlux(image,lam_obs(j),flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
			ExtISM=Reddening(lam_obs(j),compute_dlam(lam_obs(j)),Av)
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
		enddo
	else if(tel%kind(1:10).eq.'VISIBILITY') then
		readmcscat=.false.
		allocate(V(nlam_obs,tel%nbaseline))
		allocate(phase(nlam_obs,tel%nbaseline))
		call TracePath(image,angle,tel%nphi,tel%nr,0.55d0)
		if(tel%fits) then
			if(tel%noangle) then
				write(specfile,'(a,"visibility",a,".fits.gz")') outdir(1:len_trim(outdir))
     &			,tel%flag(1:len_trim(tel%flag))
			else
				write(specfile,'(a,"visibility",i1,f3.1,a,".fits.gz")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
			endif
		else
			if(tel%noangle) then
				write(specfile,'(a,"visibility",a,".dat")') outdir(1:len_trim(outdir))
     &			,tel%flag(1:len_trim(tel%flag))
			else
				write(specfile,'(a,"visibility",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
			endif
			open(unit=30,file=specfile,RECL=6000)
			write(30,'("# column   1: wavelength")')
			write(30,'("# column   2: full disk")')
			do k=1,tel%nbaseline
				write(30,'("# column ",i3," and ",i3,": vis and phase for: baseline ",f10.3,", angle ",f10.3)') 
     1                  k+2,k+2+tel%nbaseline,tel%b(k),180d0*tel%theta(k)/pi
			enddo
		endif
		do j=1,nlam_obs
			call TraceFlux(image,lam_obs(j),spec(j),scatspec(j),specQ(j),tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
			ExtISM=Reddening(lam_obs(j),compute_dlam(lam_obs(j)),Av)
			spec(j)=1d23*spec(j)*ExtISM/D%distance**2
			do k=1,tel%nbaseline
				call Visibility(image,tel%b(k),tel%theta(k),lam_obs(j),V(j,k),phase(j,k),tel%D)
			enddo
			if(.not.tel%fits) then
				write(30,*) lam_obs(j),spec(j),V(j,1:tel%nbaseline),phase(j,1:tel%nbaseline)
				call flush(30)
			endif
		enddo
		if(tel%fits) then
			call visibilityfits(specfile,lam_obs,tel%b,180d0*tel%theta/pi,spec,V,phase,nlam_obs,tel%nbaseline)
		else
			close(unit=30)
		endif
	else if(tel%kind(1:7).eq.'BASEVIS') then
c is still without interstellar extinction
c       Gijsexp: 
		nbase=100
		allocate(basegrid(nbase))
		allocate(V(tel%nlam,nbase))
		allocate(phase(tel%nlam,nbase))

		readmcscat=.false.
		call TracePath(image,angle,tel%nphi,tel%nr,tel%lam1)
		if(tel%noangle) then
			write(specfile,'(a,"basevis",a,".dat")') 
     &                  outdir(1:len_trim(outdir))
     &			,tel%flag(1:len_trim(tel%flag))
		else
			write(specfile,'(a,"basevis",i1,f3.1,a,".dat")') 
     &                  outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		endif
		open(unit=30,file=specfile,RECL=6000)
		write(30,'("# column   1: baseline")')
		write(30,'("# column   2: full disk (only last wavelength)")')
		do j=1,tel%nlam
		   write(30,'("# column ",i3," and ",i3,": vis and phase for wavelength ",f10.3,", angle ",f10.3)') 
     &		j+2,j+2+tel%nlam,tel%lam(j),180d0*tel%theta(j)/pi
		enddo
c$$$		write(30,'("# column ",i3," and further: phase")') tel%nlam+3
c		Call TraceFlux(image,tel%lam1,flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		do j=1,tel%nlam
		   if (j.eq.1.or.(tel%lam(j).ne.tel%lam(j-1))) then
		      call TraceFlux(image,tel%lam(j),flux,scatflux,fluxQ,tel%Nphot,tel%NphotAngle,tel%opening,angle,tel%mask,tel%wmask)
		   endif

		   do k=1,nbase
c		   basegrid(k)=tel%b(1)*(tel%b(2)/tel%b(1))**((k-1d0)/(nbase-1d0))  ! log grid
		      basegrid(k)=tel%b(1)+(tel%b(2)-tel%b(1))*((k-1d0)/(nbase-1d0)) ! linear grid
		      call Visibility(image,basegrid(k),tel%theta(j),tel%lam(j),V(j,k),phase(j,k),tel%D)
		   enddo
		enddo
		do k=1,nbase
		   write(30,*) basegrid(k),1d23*flux/D%distance**2,V(1:tel%nlam,k),phase(1:tel%nlam,k)
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
		if(tel%noangle) then
			write(specfile,'(a,"FWHM",a,".dat")') outdir(1:len_trim(outdir))
     &			,tel%flag(1:len_trim(tel%flag))
		else
			write(specfile,'(a,"FWHM",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		endif
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
		if(tel%noangle) then
			write(specfile,'(a,"FWHM",a,".dat")') outdir(1:len_trim(outdir))
     &			,tel%flag(1:len_trim(tel%flag))
		else
			write(specfile,'(a,"FWHM",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		endif
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
c     GFORTRAN i to i0
			write(30,'("# transition nr ",i0)') j
			write(30,'("# up, low       ",i0,i0)') i_up,i_low
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
	else if(tel%kind(1:11).eq.'VISIBLEMASS') then
		readmcscat=tel%readmcscat
		call TracePath(image,angle,tel%nphi,tel%nr,0.55d0)
		if(tel%noangle) then
			write(specfile,'(a,"visiblemass",a,".dat")') outdir(1:len_trim(outdir))
     &			,tel%flag(1:len_trim(tel%flag))
		else
			write(specfile,'(a,"visiblemass",i1,f3.1,a,".dat")') outdir(1:len_trim(outdir))
     &			,int((tel%angle)/10d0),tel%angle-10d0*int((tel%angle/10d0))
     &			,tel%flag(1:len_trim(tel%flag))
		endif
		open(unit=30,file=specfile,RECL=6000)
		call TraceVisibleMass(image,tel%lam1,flux)
		write(30,*) tel%lam1,flux
		do j=1,nlam_obs
		   if(lam_obs(j).gt.tel%lam1.and.lam_obs(j).lt.tel%lam2) then
				call TraceVisibleMass(image,lam_obs(j),flux)
				write(30,*) lam_obs(j),flux
				call flush(30)
			endif
		enddo
		call TraceVisibleMass(image,tel%lam2,flux)
		write(30,*) tel%lam2,flux
		close(unit=30)
	else if(tel%kind(1:12).eq.'OPTICALDEPTH') then
		readmcscat=.false.
		call TracePath(image,angle,tel%nphi,tel%nr,tel%lam1)
		call TraceOpticalDepth(image,tel%lam1)
		do i=1,tel%nfov
			if(tel%scaletype.eq.1) then
				image%rscale=1d0
				image%zscale=1d0
			else if(tel%scaletype.eq.2) then
				tel%fov(i)=tel%fov(i)*D%distance/parsec
				image%rscale=parsec/D%distance
				image%zscale=1d0
			endif
			call MakeSimpleImage(image,tel,tel%fov(i)/2d0,'OD')
		enddo
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
	def%lamfile=' '
	def%fits=.false.
	def%noangle=.false.
	def%vis_norm_fov=vis_norm_fov
	
1	call ignorestar(20)
c	GFORTRAN removed the format string, otherwise ther are problems with gfortran and EOF
c   I dont't think the format string is needed here
c	Michiel: Format string is needed to avoid complications when there is a space between
c            the '=' sign and the value
	read(20,'(a500)',end=2) line
c	read(20,*,end=2) line

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
		else if(tel(nobs)%kind.eq.'VISIBLEMASS') then
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
		else if(tel(nobs)%kind.eq.'OPTICALDEPTH') then
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
			tel(nobs)%width=def%width
			tel(nobs)%D=def%D
			tel(nobs)%D2=def%D2
			tel(nobs)%spider=def%spider
			tel(nobs)%vis_norm_fov=def%vis_norm_fov
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
			tel(nobs)%width=def%width
			tel(nobs)%D=def%D
			tel(nobs)%D2=def%D2
			tel(nobs)%spider=def%spider
			tel(nobs)%vis_norm_fov=def%vis_norm_fov
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
		tel(nobs)%lamfile=def%lamfile
		tel(nobs)%fits=def%fits
		tel(nobs)%noangle=def%noangle
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
		if(key.eq.'fits') read(value,*) def%fits
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
		if(key.eq.'lamfile') def%lamfile=value
		if(key.eq.'popfile') def%popfile=value
		if(key.eq.'trans_nr') read(value,*) def%trans_nr1
		if(key.eq.'trans_nr1') read(value,*) def%trans_nr1
		if(key.eq.'trans_nr2') read(value,*) def%trans_nr2
		if(key.eq.'nvelo') read(value,*) def%nvelo
		if(key.eq.'dvelo') read(value,*) def%dvelo
		if(key.eq.'abun') read(value,*) def%abun
		if(key.eq.'fastobs') read(value,*) def%fastobs
		if(key.eq.'noangle') read(value,*) def%noangle
		if(key.eq.'vis_norm_fov') read(value,*) def%vis_norm_fov
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
		if(key.eq.'fits') read(value,*) tel(nobs)%fits
		if(key(1:4).eq.'base'.and.(tel(nobs)%kind.eq.'VISIBILITY'
     &		   .or.tel(nobs)%kind.eq.'BASEVIS')) then
			read(key(5:len_trim(key)),*) i
			if(i.gt.tel(nobs)%nbaseline) then
				allocate(temp(i))
				do j=1,tel(nobs)%nbaseline
					temp(j)=tel(nobs)%b(j)
				enddo
				deallocate(tel(nobs)%b)
				allocate(tel(nobs)%b(i))
				do j=1,tel(nobs)%nbaseline
					tel(nobs)%b(j)=temp(j)
				enddo
				read(value,*) tel(nobs)%b(i)

				do j=1,tel(nobs)%nbaseline
					temp(j)=tel(nobs)%theta(j)
				enddo
				deallocate(tel(nobs)%theta)
				allocate(tel(nobs)%theta(i))
				do j=1,tel(nobs)%nbaseline
					tel(nobs)%theta(j)=temp(j)
				enddo

				tel(nobs)%nbaseline=i
				deallocate(temp)
			else
				read(value,*) tel(nobs)%b(i)
			endif
		endif
		if(key(1:5).eq.'angle'.and.(tel(nobs)%kind.eq.'VISIBILITY'
     &		   .or.tel(nobs)%kind.eq.'BASEVIS')) then
			read(key(6:len_trim(key)),*) i
			if(i.gt.tel(nobs)%nbaseline) then
				allocate(temp(i))
				do j=1,tel(nobs)%nbaseline
					temp(j)=tel(nobs)%theta(j)
				enddo
				deallocate(tel(nobs)%theta)
				allocate(tel(nobs)%theta(i))
				do j=1,tel(nobs)%nbaseline
					tel(nobs)%theta(j)=temp(j)
				enddo
				read(value,*) tel(nobs)%theta(i)

				do j=1,tel(nobs)%nbaseline
					temp(j)=tel(nobs)%b(j)
				enddo
				deallocate(tel(nobs)%b)
				allocate(tel(nobs)%b(i))
				do j=1,tel(nobs)%nbaseline
					tel(nobs)%b(j)=temp(j)
				enddo

				tel(nobs)%nbaseline=i
				deallocate(temp)
			else
				read(value,*) tel(nobs)%theta(i)
			endif
			tel(nobs)%theta(i)=pi*tel(nobs)%theta(i)/180d0
		endif
c       Gijsexp: allow more than two wavelength/angle combo's for basevis 
		if(key(1:3).eq.'lam'.and.tel(nobs)%kind.eq.'BASEVIS'.and.key.ne.'lamfile') then
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
		if(key.eq.'lamfile') tel(nobs)%lamfile=value
		if(key.eq.'popfile') tel(nobs)%popfile=value
		if(key.eq.'trans_nr') read(value,*) tel(nobs)%trans_nr1
		if(key.eq.'trans_nr1') read(value,*) tel(nobs)%trans_nr1
		if(key.eq.'trans_nr2') read(value,*) tel(nobs)%trans_nr2
		if(key.eq.'nvelo') read(value,*) tel(nobs)%nvelo
		if(key.eq.'dvelo') read(value,*) tel(nobs)%dvelo
		if(key.eq.'abun') read(value,*) tel(nobs)%abun
		if(key.eq.'fastobs') read(value,*) tel(nobs)%fastobs
		if(key.eq.'noangle') read(value,*) tel(nobs)%noangle
		if(key.eq.'vis_norm_fov') read(value,*) tel(nobs)%vis_norm_fov
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
			p%phi1(i)=opening*sign(1d0,p%phi1(i))
			p%jphi1(i)=real(NPHISCATT)*p%phi1(i)/360d0
			p%jphi1(i)=p%jphi1(i)+1
			if(p%jphi1(i).lt.1) p%jphi1(i)=1-p%jphi1(i)
			if(p%jphi1(i).gt.NPHISCATT/2) p%jphi1(i)=NPHISCATT+1-p%jphi1(i)
			p%phi2(i)=opening*sign(1d0,p%phi2(i))
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
	

	subroutine TraceVisibleMass(image,lam0,vmass)	
	use Parameters
	IMPLICIT NONE
	type(RPhiImage) image
	real*8 lam0,vmass,tau_e,wl1,wl2,w1,w2,tau0,Kext,dens(0:D%nR,0:D%nTheta)
	real*8 ww(ngrains,image%nr,image%nphi),www(ngrains),ww0(ngrains)
	integer i,j,k,ip,jp,ilam1,ilam2,iopac,ii

	write(*,*) lam0

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

!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,j,ii,iopac)
!$OMP& SHARED(dens,C,Grain,wl1,wl2,ilam1,ilam2,ngrains,D)
!$OMP DO
	do i=0,D%nR-1
	do j=1,D%nTheta-1
		C(i,j)%Kext=0d0
		dens(i,j)=0d0
		do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				C(i,j)%Kext=C(i,j)%Kext+(wl1*Grain(ii)%Kext(iopac,ilam1)
     &			+wl2*Grain(ii)%Kext(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
			enddo
			if(Grain(ii)%rv.lt.20d-4.and.C(i,j)%T.gt.30d0) then
				dens(i,j)=dens(i,j)+C(i,j)%dens*C(i,j)%w(ii)
			endif
		enddo
	enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL


	ww=0d0
!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,j,tau0,ip,jp,Kext,tau_e,ii)
!$OMP& SHARED(image,dens,C,ww,ngrains)
!$OMP DO
	do i=1,image%nr
	do j=1,image%nphi
		tau0=0d0
		image%image(i,j)=0d0
		do k=1,image%p(i,j)%n

			ip=image%p(i,j)%i(k)
			jp=image%p(i,j)%j(k)

			Kext=C(ip,jp)%Kext

			tau_e=image%p(i,j)%v(k)*C(ip,jp)%dens*Kext*AU
			if((tau0+tau_e).gt.1d0) then
				image%image(i,j)=image%image(i,j)+image%p(i,j)%v(k)*dens(ip,jp)*AU*(1d0-tau0)/tau_e
				do ii=1,ngrains
					ww(ii,i,j)=ww(ii,i,j)+image%p(i,j)%v(k)*C(ip,jp)%dens*C(ip,jp)%w(ii)*AU*(1d0-tau0)/tau_e
				enddo
				goto 1
			else
				image%image(i,j)=image%image(i,j)+image%p(i,j)%v(k)*dens(ip,jp)*AU
				do ii=1,ngrains
					ww(ii,i,j)=ww(ii,i,j)+image%p(i,j)%v(k)*C(ip,jp)%dens*C(ip,jp)%w(ii)*AU
				enddo
			endif
			tau0=tau0+tau_e
		enddo
1		continue
	enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	vmass=0d0
	www=0d0
	do i=1,image%nr-1
	do k=1,image%nPhi
		w1=2d0*pi*abs(image%R(i))*AU**2/real(image%nPhi)
		w2=2d0*pi*abs(image%R(i+1))*AU**2/real(image%nPhi)
		vmass=vmass+(image%R(i+1)-image%R(i))*
     &		(w1*image%image(i,k)+w2*image%image(i+1,k))/2d0
		do ii=1,ngrains
			www(ii)=www(ii)+(image%R(i+1)-image%R(i))*
     &		(w1*ww(ii,i,k)+w2*ww(ii,i+1,k))/2d0
		enddo
	enddo
	enddo

	ww0=0d0
	do i=1,D%nR-1
		do j=1,D%nTheta-1
			ww0=ww0+C(i,j)%dens*C(i,j)%V*C(i,j)%w
		enddo
	enddo

	write(85,*) lam0,www(1:ngrains)/ww0(1:ngrains)

	vmass=vmass/Msun
	
	return
	end
	


	subroutine visibilityfits(filename,lam,base,theta,spec,V,phase,nlam,nbase)
	IMPLICIT NONE
	character*500 filename
	integer nlam,nbase
	real*8 lam(nlam),spec(nlam),base(nbase),theta(nbase),V(nlam,nbase),phase(nlam,nbase)
	real*8 tb(nbase,2)

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

C     initialize parameters about the FITS image (64-bit reals)
	simple=.true.
	bitpix=-64
	extend=.true.
	group=1
	fpixel=1
	naxis=2
	naxes(1)=nbase
	naxes(2)=2
	nelements=naxes(1)*naxes(2)

C		write the required header keywords
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	call ftpkyj(unit,'nlam',nlam,' ',status)
	call ftpkyj(unit,'nbase',nbase,' ',status)

	!------------------------------------------------------------------------------
	! HDU 0: baseline grid
	!------------------------------------------------------------------------------

	do i=1,nbase
		tb(i,1)=base(i)
		tb(i,2)=theta(i)
	enddo

C		Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,tb(1:nbase,1:2),status)


	!------------------------------------------------------------------------------
	! HDU 1: wavelength grid
	!------------------------------------------------------------------------------

C		create new hdu
	call ftcrhd(unit, status)
C     initialize parameters about the FITS image (64-bit reals)
	simple=.true.
	bitpix=-64
	extend=.true.
	group=1
	fpixel=1
	naxis=1
	naxes(1)=nlam
	naxes(2)=1
	nelements=naxes(1)*naxes(2)

C		write the required header keywords
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
C		Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,lam(1:nlam),status)


	!------------------------------------------------------------------------------
	! HDU 2: full disk spectrum
	!------------------------------------------------------------------------------

C		create new hdu
	call ftcrhd(unit, status)
C     initialize parameters about the FITS image (64-bit reals)
	simple=.true.
	bitpix=-64
	extend=.true.
	group=1
	fpixel=1
	naxis=1
	naxes(1)=nlam
	naxes(2)=1
	nelements=naxes(1)*naxes(2)

C		write the required header keywords
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
C		Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,spec(1:nlam),status)


	!------------------------------------------------------------------------------
	! HDU 3: visibilities
	!------------------------------------------------------------------------------

C		create new hdu
	call ftcrhd(unit, status)
C     initialize parameters about the FITS image (64-bit reals)
	simple=.true.
	bitpix=-64
	extend=.true.
	group=1
	fpixel=1
	naxis=2
	naxes(1)=nlam
	naxes(2)=nbase
	nelements=naxes(1)*naxes(2)

C		write the required header keywords
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
C		Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,V(1:nlam,1:nbase),status)


	!------------------------------------------------------------------------------
	! HDU 4: phases
	!------------------------------------------------------------------------------

C		create new hdu
	call ftcrhd(unit, status)
C     initialize parameters about the FITS image (64-bit reals)
	simple=.true.
	bitpix=-64
	extend=.true.
	group=1
	fpixel=1
	naxis=2
	naxes(1)=nlam
	naxes(2)=nbase
	nelements=naxes(1)*naxes(2)

C		write the required header keywords
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
C		Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,phase(1:nlam,1:nbase),status)





C     close the file and free the unit number
	call ftclos(unit, status)
	call ftfiou(unit, status)
	
	
	return
	end
	
	


