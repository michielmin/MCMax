c	2008-01-21: MM changed the refinement of the grid to happen
c				also when only the vertical structure is iterated
c	2008-02-20: Added the possibility for mixed aggregates with a 
c				powerlaw gradient in composition
c	2008-05-23: - NaN error problem solved.
c				- I think the 500 gridpoint limit is solved.
c				- Backwarming included.
c				- Removed the absolute maximum taustep of 5 to decrease the number of refinement cells needed.
c				- Added a limit to the optical depth increase of a factor of 2 throughout the disk.
c	2008-07-03:	Fixed a bug which prevented the Random Walk procedure to be
c				used when thermal contact was off.
c	2008-07-11:	Fixed a bug introduced about two weeks ago causing the
c				absorption to be twice as high.
c	2008-08-28:	Added the possibility of tracing PAH emission.
c	2008-09-02:	Added the possibility of multi-photon PAH excitation
c				Multi-photon PAHs are now treated fully selfconsistent
c	2008-09-03:	Several bugfixes (among others in the regridtheta subroutine.)
c	2008-09-23:	Changed the fweight to be spatially varying. According to suggestions
c				by Mihkel. Now is a powerlaw between the optically thin evaporation region
c				and the backwarming evaporation region.
c	2009-04-22:	Changed the gasdens output to the denstemp file to include the gas2dust ratio
c	2009-04-22:	Ngrains is now output to the denstemp file

	program MCMax
	use Parameters
	use NAG
	IMPLICIT NONE
	character*500 input,denstempfile,tau1file,tmp,tau1fileR,version
	character*500 comm,gasdenstempfile
	character*500 photfile,output,pshfile,logfile,kappafile,errfile
	character*500 denstempfileND,radgridfile,qhpfile
	character*500,allocatable :: denstempfileP(:)
	character*10 date,time
	real*8 determineT,tottime,errT,errRho,radialtau,tau,radtau,Rmax,Mass,tau0
	real*8 starttime,stoptime,starttime0,determineTP,errE,Rdes,Rmin
	real*8,allocatable :: DustMass(:)
	integer i,j,ii,Nphot,ncor,niter,ntau1,l,iT,NphotFirst,NFirst,iopac
	integer number_invalid
	integer,allocatable :: icor(:),jcor(:)
	logical finaliter,prevalloc,coralloc,truefalse,lastiter,ComputeLRF_backup
	type(Photon) phot
	type(Cell),allocatable :: Cprev(:,:)
	type(Telescope) tel(MAXOBS)
	real*8 T,kappa,KappaGas,w(100)
	
c The version number (please update frequently!)
	call VersionDateTime(version)

	arraysallocated=.false.
	prevalloc=.false.
	coralloc=.false.
	
	write(outdir,'("./")')
	i=1
3	call get_command_argument(i,tmp)
	if(tmp.eq.' ') goto 4
	if(tmp.eq.'-o') then
		call get_command_argument(i+1,tmp)
		write(outdir,'(a,"/")') tmp(1:len_trim(tmp))
		write(tmp,'("mkdir -p ",a)') outdir(1:len_trim(outdir)-1)
		call system(tmp)
	endif
	i=i+1
	goto 3
4	continue

	call date_and_time(date,time)
	write(logfile,'(a,"log.dat")') outdir(1:len_trim(outdir))
	inquire(file=logfile,exist=truefalse)
	if(truefalse) then
		write(comm,'("cat ",a,"log.dat >> ",a,"logall.dat")') outdir(1:len_trim(outdir)),outdir(1:len_trim(outdir))
		write(*,*) comm(1:len_trim(comm))
		call system(comm)
	endif
		open(unit=9,file=logfile,RECL=6000)
		write(*,'("Creating log file:  ",a)') logfile(1:len_trim(logfile))
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		write(*,'("Date: ",a2,"/",a2,"/",a4)') date(7:8),date(5:6),date(1:4)
		write(*,'("Time: ",a2,":",a2,":",a6)') time(1:2),time(3:4),time(5:10)
		write(9,'("Date: ",a2,"/",a2,"/",a4)') date(7:8),date(5:6),date(1:4)
		write(9,'("Time: ",a2,":",a2,":",a6)') time(1:2),time(3:4),time(5:10)
c	endif
	write(*,'("Version: ",a)') version(1:len_trim(version))
	write(9,'("Version: ",a)') version(1:len_trim(version))
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	call cpu_time(starttime)
	starttime0=starttime
	nemit=0
	nmaxinteract=0

	idum=-42
	
	call get_command_argument(1,input)
	call get_command_argument(2,tmp)
	read(tmp,*) Nphot

		
	call initialize(input,Nphot,NphotFirst,NFirst,niter)	

	nobs=0
	i=0
2	call get_command_argument(3+i,output)
	if(output.ne.' ') then
		if(output(1:2).eq.'-o'.or.output(1:2).eq.'-s'.or.output(1:2).eq.'-p') then 
			i=i+2
			goto 2
		endif
		call determineoutput(output,tel)
		i=i+1
		goto 2
	endif
	if(nobs.eq.0) then
		write(*,'("No output observations")')
		write(9,'("No output observations")')
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
	endif

	if(.not.tracestar) then
		do i=1,nobs
			if(tel(i)%flag.eq.' ') then
				tel(i)%flag='NOSTAR'
			else
				tel(i)%flag(len_trim(tel(i)%flag)+1:len_trim(tel(i)%flag)+8)='_NOSTAR'
			endif
		enddo
	endif
	

	Mass=0d0
	do i=0,D%nR-1
	do j=1,D%nTheta-1
		Mass=Mass+C(i,j)%gasdens*C(i,j)%V*gas2dust
	enddo
	enddo
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Total mass in the disk:  ",e10.4," Msun")') Mass/Msun
	write(9,'("Total mass in the disk:  ",e10.4," Msun")') Mass/Msun
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	call cpu_time(stoptime)
	tottime=stoptime-starttime
	write(*,'("Initialize time:   ",f11.4,"  s")') tottime
	write(9,'("Initialize time:   ",f11.4,"  s")') tottime
	write(9,*)
	call flush(9)

	finaliter=.false.
	lastiter=.false.
	if(.not.tcontact) allocate(denstempfileP(ngrains))

	allocate(DustMass(0:ngrains))

1	continue

	DustMass=0d0
	do i=0,D%nR-1
	do j=1,D%nTheta-1
		DustMass(0)=DustMass(0)+C(i,j)%dens*C(i,j)%V
		do ii=1,ngrains
			DustMass(ii)=DustMass(ii)+C(i,j)%dens*C(i,j)%V*C(i,j)%w(ii)
		enddo
	enddo
	enddo
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Dust mass in the disk:      ",e10.4," Msun")') DustMass(0)/Msun
	write(9,'("Dust mass in the disk:      ",e10.4," Msun")') DustMass(0)/Msun
	do ii=1,ngrains
	write(*,'("Dust mass in component",i4,": ",e10.4," Msun")') ii,DustMass(ii)/Msun
	write(9,'("Dust mass in component",i4,": ",e10.4," Msun")') ii,DustMass(ii)/Msun
	enddo
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	radtau=radialtau(1.0d0,tau,ntau1,D%nTheta-1)
	write(*,'("Radial optical depth at 1 micron:   ",f15.3)') radtau
	write(*,'("First cell:                         ",f15.3)') tau
	write(*,'("Number of cells from 0.2 - 1.0:     ",i15)') ntau1
	write(9,'("Radial optical depth at 1 micron:   ",f15.3)') radtau
	write(9,'("First cell:                         ",f15.3)') tau
	write(9,'("Number of cells from 0.2 - 1.0:     ",i15)') ntau1
	radtau=radialtau(0.55d0,tau,ntau1,D%nTheta-1)
	write(*,'("Radial optical depth in the visual: ",f15.3)') radtau
	write(*,'("First cell:                         ",f15.3)') tau
	write(*,'("Number of cells from 0.2 - 1.0:     ",i15)') ntau1
	write(9,'("Radial optical depth in the visual: ",f15.3)') radtau
	write(9,'("First cell:                         ",f15.3)') tau
	write(9,'("Number of cells from 0.2 - 1.0:     ",i15)') ntau1
	radtau=radialtau(0.55d0,tau,ntau1,1)
	write(*,'("Upwards optical depth in the visual:",f15.3)') radtau
	write(9,'("Upwards optical depth in the visual:",f15.3)') radtau
	tau0=radtau
	do i=1,D%nTheta-1
		radtau=radialtau(0.55d0,tau,ntau1,i)
		if(radtau.lt.tau0) tau0=radtau
	enddo
	write(*,'("Minimum optical depth in the visual:",f15.3)') tau0
	write(9,'("Minimum optical depth in the visual:",f15.3)') tau0

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	radtau=radialtau(9.8d0,tau,ntau1,D%nTheta-1)
	write(*,'("Radial optical depth at 9.8 micron: ",f15.3)') radtau
	write(9,'("Radial optical depth at 9.8 micron: ",f15.3)') radtau
	radtau=radialtau(9.8d0,tau,ntau1,1)
	write(*,'("Upwards optical depth at 9.8 micron:",f15.3)') radtau
	write(9,'("Upwards optical depth at 9.8 micron:",f15.3)') radtau

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	if((struct_iter.or.tdes_iter.or.use_qhp.or.use_topac).and.Nphot.gt.0.and..not.lastiter.and.maxiter.gt.0) then
		if(niter.lt.100) then
		if(outputfits) then
			write(denstempfile,'(a,"denstemp",i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
			write(denstempfileND,'(a,"denstempND",i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
			write(gasdenstempfile,'(a,"denstempGas",i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
			write(photfile,'(a,"Nphotons",i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
			write(errfile,'(a,"errors",i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
		else
			write(denstempfile,'(a,"denstemp",i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
			write(denstempfileND,'(a,"denstempND",i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
			write(gasdenstempfile,'(a,"denstempGas",i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
			write(photfile,'(a,"Nphotons",i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
			write(errfile,'(a,"errors",i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
		endif
		write(radgridfile,'(a,"radgrid",i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
		write(tau1file,'(a,"height",i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
		write(tau1fileR,'(a,"heightR",i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
		write(pshfile,'(a,"scaleheight",i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/10,niter-10*(niter/10)
		else
		if(outputfits) then
			write(denstempfile,'(a,"denstemp",i1,i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
			write(denstempfileND,'(a,"denstempND",i1,i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
			write(gasdenstempfile,'(a,"denstempGas",i1,i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
			write(photfile,'(a,"Nphotons",i1,i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
			write(errfile,'(a,"errors",i1,i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
		else
			write(denstempfile,'(a,"denstemp",i1,i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
			write(denstempfileND,'(a,"denstempND",i1,i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
			write(gasdenstempfile,'(a,"denstempGas",i1,i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
			write(photfile,'(a,"Nphotons",i1,i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
			write(errfile,'(a,"errors",i1,i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
		endif
		write(radgridfile,'(a,"radgrid",i1,i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
		write(tau1file,'(a,"height",i1,i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
		write(tau1fileR,'(a,"heightR",i1,i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
		write(pshfile,'(a,"scaleheight",i1,i1,i1,".dat")') outdir(1:len_trim(outdir)),niter/100,(niter-100*(niter/100))/10
     &																,niter-10*(niter/10)
		endif
		write(*,'("Writing tau=1 height to:      ",a)') tau1file(1:len_trim(tau1file))
		write(9,'("Writing tau=1 height to:      ",a)') tau1file(1:len_trim(tau1file))
		write(*,'("Writing tau=1 height (radial):",a)') tau1fileR(1:len_trim(tau1file))
		write(9,'("Writing tau=1 height (radial):",a)') tau1fileR(1:len_trim(tau1file))
		write(*,'("Writing radial grid to:       ",a)') radgridfile(1:len_trim(radgridfile))
		write(9,'("Writing radial grid to:       ",a)') radgridfile(1:len_trim(radgridfile))
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		if(niter.eq.0) then
			write(*,'("Obtaining initial guess  (",i4,"/",i4,")")') niter,maxiter
			write(9,'("Obtaining initial guess  (",i4,"/",i4,")")') niter,maxiter
		else
			write(*,'("Iterating disk structure (",i4,"/",i4,")")') niter,maxiter
			write(9,'("Iterating disk structure (",i4,"/",i4,")")') niter,maxiter
		endif
	else
		if(outputfits) then
			write(denstempfile,'(a,"denstemp.fits.gz")')outdir(1:len_trim(outdir))
			write(denstempfileND,'(a,"denstempND.fits.gz")')outdir(1:len_trim(outdir))
			write(gasdenstempfile,'(a,"denstempGas.fits.gz")')outdir(1:len_trim(outdir))
			write(photfile,'(a,"Nphotons.fits.gz")')outdir(1:len_trim(outdir))
			write(errfile,'(a,"errors.fits.gz")')outdir(1:len_trim(outdir))
		else
			write(denstempfile,'(a,"denstemp.dat")')outdir(1:len_trim(outdir))
			write(denstempfileND,'(a,"denstempND.dat")')outdir(1:len_trim(outdir))
			write(gasdenstempfile,'(a,"denstempGas.dat")')outdir(1:len_trim(outdir))
			write(photfile,'(a,"Nphotons.dat")')outdir(1:len_trim(outdir))
			write(errfile,'(a,"errors.dat")')outdir(1:len_trim(outdir))
		endif
		write(radgridfile,'(a,"radgrid.dat")') outdir(1:len_trim(outdir))
		write(tau1file,'(a,"height.dat")')outdir(1:len_trim(outdir))
		write(tau1fileR,'(a,"heightR.dat")')outdir(1:len_trim(outdir))
		write(pshfile,'(a,"scaleheight.dat")')outdir(1:len_trim(outdir))
		write(*,'("Writing tau=1 height to:      ",a)') tau1file(1:len_trim(tau1file))
		write(9,'("Writing tau=1 height to:      ",a)') tau1file(1:len_trim(tau1file))
		write(*,'("Writing tau=1 height (radial):",a)') tau1fileR(1:len_trim(tau1file))
		write(9,'("Writing tau=1 height (radial):",a)') tau1fileR(1:len_trim(tau1file))
		write(*,'("Writing radial grid to:       ",a)') radgridfile(1:len_trim(radgridfile))
		write(9,'("Writing radial grid to:       ",a)') radgridfile(1:len_trim(radgridfile))
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		if(niter.ne.0) then
			write(*,'("Final computation")')
			write(9,'("Final computation")')
			if(NphotFinal.ne.0) Nphot=NphotFinal
		endif
	endif

	if(.not.tcontact) then
		do ii=1,ngrains
			write(denstempfileP(ii),'(a,"P",i1,i1,".dat")') denstempfile(1:len_trim(denstempfile)-4),ii/10,ii-10*(ii/10)
		enddo
	endif
	
	if(Nphot.le.0) then
		call reset2d()
		call ReadTemperatures()
		finaliter=.true.
	endif

	if(Nphot.gt.0.and.use_qhp) then
c		call PAHMCMax(niter)
		call Stochastic(niter)
		call destroyQHP(TdesQHP)
	endif

	call tau1height(tau1file)
	call tau1heightR(tau1fileR)
	open(unit=20,file=radgridfile,RECL=100)
	do i=1,D%nR
		write(20,*) D%R(i)
	enddo
	close(unit=20)


	write(output,'(a,"composition.dat")') outdir(1:len_trim(outdir))
	open(unit=20,file=output,RECL=6000)
	do i=1,D%nR-1
		w(1:ngrains)=0d0
		do j=1,D%nTheta-1
			w(1:ngrains)=w(1:ngrains)+C(i,j)%mass*C(i,j)%w(1:ngrains)
		enddo
		write(20,*) D%R_av(i)/AU,(C(i,D%nTheta-1)%w(ii),ii=1,ngrains),(w(ii)/sum(w(1:ngrains)),ii=1,ngrains)
	enddo
	close(unit=20)

	call cpu_time(starttime)

	if(arraysallocated.and..not.prevalloc) then
		allocate(Cprev(0:D%nR,0:D%nTheta))
		do i=0,D%nR
		do j=0,D%nTheta
			allocate(Cprev(i,j)%w(ngrains))
			allocate(Cprev(i,j)%w0(ngrains))
			allocate(Cprev(i,j)%wopac(ngrains,ngrains2))
			if(.not.tcontact.or.tdes_iter) then
				allocate(Cprev(i,j)%TP(ngrains))
				allocate(Cprev(i,j)%EJvP(ngrains))
			endif
			if(use_qhp) then
				allocate(Cprev(i,j)%QHP(nqhp,nlam))
				allocate(Cprev(i,j)%tdistr(nqhp,NTQHP))
				allocate(Cprev(i,j)%Tqhp(nqhp))
				allocate(Cprev(i,j)%EJvQHP(nqhp))
			endif
			if(computeLRF) then
				allocate(Cprev(i,j)%LRF(nlam))
				allocate(Cprev(i,j)%nLRF(nlam))
			endif
		enddo
		enddo
		prevalloc=.true.
	endif
	if(arraysallocated) Cprev(1:D%nR,1:D%nTheta)=C(1:D%nR,1:D%nTheta)

	if(Nphot.gt.0) then

	ComputeLRF_backup=ComputeLRF
	if(((struct_iter.or.tdes_iter.or.use_topac).and.Nphot.gt.0.and..not.lastiter.and.maxiter.gt.0).and..not.use_qhp) then
		ComputeLRF=.false.
	endif

	if(niter.eq.int(NFirst/2+1)) overflow=.false.
	if(niter.le.NFirst) then
		call MCRadiation(NphotFirst,niter,.true.)
	else if(viscous.and.tdes_iter.and.niter.le.nBW.and.nBW.lt.maxiter.and..not.lastiter) then
		call MCRadiation(Nphot/10,niter,.true.)
	else
		call MCRadiation(Nphot,niter,fastviscous)
	endif
	ComputeLRF=ComputeLRF_backup

	if(.not.coralloc) then
		allocate(icor((D%nTheta+1)*(D%nR+1)))
		allocate(jcor((D%nTheta+1)*(D%nR+1)))
		coralloc=.true.
	endif

	call cpu_time(stoptime)
	tottime=stoptime-starttime

	write(*,'("Trace time:        ",f11.4,"  s")') tottime
	write(*,'("Time per package:  ",f11.4," ms")') 1000d0*tottime/real(Nphot)
	write(*,'("Time per emission: ",f11.4," ms")') 1000d0*tottime/real(nemit)

	write(9,'("Trace time:        ",f11.4,"  s")') tottime
	write(9,'("Time per package:  ",f11.4," ms")') 1000d0*tottime/real(Nphot)
	write(9,'("Time per emission: ",f11.4," ms")') 1000d0*tottime/real(nemit)

	starttime=stoptime

	endif

	if(Nphot.gt.0) then
	do j=1,D%nTheta-1
		do i=1,D%nR-1
			phot%i=i
			phot%j=j
			C(i,j)%EJv=C(i,j)%EJv/C(i,j)%V
c
c  changed by AJ van Marle
c   was not allocated if use_qhp == .false.		
			if( use_qhp ) C(i,j)%EJvQHP=C(i,j)%EJvQHP/C(i,j)%V
c
c
c			
			C(i,j)%EJv2=C(i,j)%EJv2/C(i,j)%V**2
			phot%E=C(i,j)%EJv
			C(i,j)%T=determineT(phot)
			if(C(i,j)%T.lt.1d0) C(i,j)%T=1d0
			if(C(i,j)%Ni.gt.1) then
				C(i,j)%dEJv=(C(i,j)%EJv
     &				+sqrt(abs(C(i,j)%EJv2-C(i,j)%EJv**2/real(C(i,j)%Ni))))
     &				/sqrt(real(C(i,j)%Ni))
				C(i,j)%dT=0.25d0*C(i,j)%T*C(i,j)%dEJv/C(i,j)%EJv
			else
				C(i,j)%dT=real(TMAX)*dT
				C(i,j)%dEJv=0d0
				do ii=1,ngrains
					do iopac=1,Grain(ii)%nopac
						C(i,j)%dEJv=C(i,j)%dEJv+C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)*Grain(ii)%Kp(iopac,TMAX)
					enddo
				enddo
			endif
			if(number_invalid(C(i,j)%T).ne.0) then
				print*,i,j,C(i,j)%Ni
				C(i,j)%T=1d0
			endif
			if(number_invalid(C(i,j)%dT).ne.0) C(i,j)%dT=dT
			if(.not.tcontact.or.tdes_iter) then
			do ii=1,ngrains
				phot%i=i
				phot%j=j
				C(i,j)%EJvP(ii)=C(i,j)%EJvP(ii)/C(i,j)%V
				phot%E=C(i,j)%EJvP(ii)
				C(i,j)%TP(ii)=determineTP(phot,ii)
				if(C(i,j)%TP(ii).lt.1d0) C(i,j)%TP(ii)=1d0
				if(number_invalid(C(i,j)%TP(ii)).ne.0.or.C(i,j)%randomwalk) C(i,j)%TP(ii)=C(i,j)%T
			enddo
			endif			
		enddo
		C(0,j)%T=D%Tstar*0.5d0*(1d0-sqrt(1d0-(D%Rstar/(D%R(1)*AU))**2))
		if(C(0,j)%T.lt.1d0) C(0,j)%T=1d0
		C(0,j)%dT=dT
	enddo
	if(convection.and.niter.gt.NFirst) call MakeAdiabatic(2d0/7d0)
	endif

	do j=1,D%nTheta-1
		do i=1,D%nR-1
			C(i,j)%TMC=C(i,j)%T
		enddo
	enddo
	if(Tsmooth) call do_smooth()

	if((NphotDiffuse.gt.0.or.dTDiffuse.lt.1d0).and.(Nphot.ne.0.or.forcediff)) then
		call Diffuse(NPhotDiffuse,dTDiffuse)
	endif
	if(Nphot.gt.0) then
		if(computeTgas.or.viscous) then
			do i=1,D%nR-1
				call tellertje(i,D%nR-1)
				do j=1,D%nTheta-1
					call determineTgas(i,j)
				enddo
			enddo
		endif
	endif


	if(tdes_iter) then
		do i=0,D%nR-1
		do j=1,D%nTheta-1
			if(C(i,j)%Ni.gt.10) then
				C(i,j)%KextLRF=C(i,j)%KextLRF/C(i,j)%ILRF
			else
				iT=C(i,j)%T/dT
				call DiffCoeffCell(i,j,iT)
				C(i,j)%KextLRF=C(i,j)%KDext
			endif
		enddo
		enddo
	endif

	
	write(kappafile,'(a,"kappas.dat")') outdir(1:len_trim(outdir))
	open(unit=20,file=kappafile,RECL=6000)
	do i=1,nlam
		write(20,*) lam(i),(Grain(j)%Kabs(1,i),j=1,ngrains),(Grain(j)%Ksca(1,i),j=1,ngrains)
	enddo
	close(unit=20)

	if(Nphot.gt.0.or.forcediff) then
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		write(*,'("Writing density and temperature structure to: ",a)') denstempfile(1:len_trim(denstempfile))
		write(9,'("Writing density and temperature structure to: ",a)') denstempfile(1:len_trim(denstempfile))
		call outputstruct(denstempfile,(/'DENS   ','TEMP   ','COMP   ','GASDENS','DENS0  '/),5,0)
	
		if(useTgas.or.viscous) then
			write(*,'("--------------------------------------------------------")')
			write(9,'("--------------------------------------------------------")')
			write(*,'("Writing gas density and temperature structure to: ",a)') gasdenstempfile(1:len_trim(gasdenstempfile))
			write(9,'("Writing gas density and temperature structure to: ",a)') gasdenstempfile(1:len_trim(gasdenstempfile))
			call outputstruct(gasdenstempfile,(/'GASDENS','GASTEMP'/),2,0)
		endif

		if(NphotDiffuse.gt.0.or.dTDiffuse.lt.1d0) then
			write(*,'("--------------------------------------------------------")')
			write(9,'("--------------------------------------------------------")')
			write(*,'("Writing density and MC temp. structure to:    ",a)') denstempfileND(1:len_trim(denstempfileND))
			write(9,'("Writing density and MC temp. structure to: ",a)') denstempfileND(1:len_trim(denstempfileND))
			call outputstruct(denstempfileND,(/'DENS   ','TEMPMC ','COMP   ','GASDENS'/),4,0)
		endif

		if(.not.tcontact) then
			do ii=1,ngrains
				write(*,'("Writing density and temperature structure to: ",a)') denstempfileP(ii)(1:len_trim(denstempfileP(ii)))
				write(9,'("Writing density and temperature structure to: ",a)') denstempfileP(ii)(1:len_trim(denstempfileP(ii)))
				call outputstruct(denstempfileP(ii),(/'DENSP  ','TEMPP  '/),2,ii)
			enddo
		endif

		write(*,'("Writing photon statistics to:     ",a)') photfile(1:len_trim(photfile))
		write(9,'("Writing photon statistics to:     ",a)') photfile(1:len_trim(photfile))
		call outputstruct(photfile,(/'NPHOT  ','VOLUME '/),2,0)

		write(*,'("Writing errors to:                ",a)') errfile(1:len_trim(errfile))
		write(9,'("Writing errors to:                ",a)') errfile(1:len_trim(errfile))
		call outputstruct(errfile,(/'dTEMP  ','dEJv   '/),2,0)

		if(use_qhp) then
			do ii=1,ngrains
				if(Grain(ii)%qhp) then
					if(outputfits) then
						write(qhpfile,'(a,"QHPemis",i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),ii/10,ii-10*(ii/10)
					else
						write(qhpfile,'(a,"QHPemis",i1,i1,".dat")') outdir(1:len_trim(outdir)),ii/10,ii-10*(ii/10)
					endif
					write(*,'("Writing QHP emission to:          ",a)') qhpfile(1:len_trim(qhpfile))
					write(9,'("Writing QHP emission to:          ",a)') qhpfile(1:len_trim(qhpfile))
					call outputstruct(qhpfile,(/'LAM    ','QHPEJv ','QHPRF  '/),3,ii)
				endif
			enddo
		endif
	endif
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Writing pressure scaleheight to:  ",a)') pshfile(1:len_trim(pshfile))
	write(9,'("Writing pressure scaleheight to:  ",a)') pshfile(1:len_trim(pshfile))
	call scaleheight(pshfile)
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')



c -------------------------------------------------------------
c Compute the temperature in the midplane according to Shakura-Sunyaev disk
c -------------------------------------------------------------
	call ShakuraSunyaev()
c -------------------------------------------------------------
c -------------------------------------------------------------

	if((struct_iter.or.tdes_iter.or.use_qhp.or.use_topac).and.Nphot.ne.0.and..not.lastiter.and.maxiter.gt.0) then
		niter=niter+1
		if(niter.gt.niter0) then
			call TempAverage(dble(1d0/(niter-niter0)))
		else if(niter.gt.1) then
			call TempAverage(f_weight)
		else
			call TempAverage(1d0)
		endif
		if(struct_iter.and.gsd) call GrainsizeDistribution()  ! GijsExp
		if(struct_iter) call DiskStructure()
		do i=1,ngrains
			if(Grain(i)%parttype.eq.3) call MakeMixAggregates(i)
		enddo
		call DestroyDustR()
		if(tdes_iter) then
			if(niter.le.nBW.or.nBW.lt.0) then
				call BackWarming(1d0)
				tauincrease=1d30
				do ii=1,6	!15
					write(*,'("Iteration ",i2," of 6")') ii
					call OpticallyThin(.true.)
					call DestroyDustT(C,Rdes,Rmin,1d10,.true.)
					write(*,'("Optically thin destruction R: ",f13.6," AU")') Rmin
					write(9,'("Optically thin destruction R: ",f13.6," AU")') Rmin
					if(Rdes.gt.D%R(1)) then
						write(*,'("No dust within:               ",f13.6," AU")') Rdes
						write(9,'("No dust within:               ",f13.6," AU")') Rdes
					endif
					if(gridrefine) then
						Rmax=D%Rout
						call RegridR(Rdes,Rmax)
						RmaxRefine=Rmax
					endif
				enddo
				tauincrease=1.25d0
				call OpticallyThin(.true.)
				call DestroyDustT(C,Rdes,Rmin,1d10,.true.)
			else
c				if(((niter/2)*2).eq.(niter)) then
					call BackWarming(1d0)
					call DestroyDustT(C,Rdes,Rmin,0.1d0,.false.)
c				endif
			endif
			write(*,'("Optically thin destruction R: ",f13.6," AU")') Rdes
			write(9,'("Optically thin destruction R: ",f13.6," AU")') Rdes
			if(Rdes.gt.D%R(1)) then
				write(*,'("No dust within:               ",f13.6," AU")') Rdes
				write(9,'("No dust within:               ",f13.6," AU")') Rdes
			endif
		endif
		call compareCells(Cprev,C,errT,errRho,errE)
c Regridding now only when error larger than 2 times epsiter
		if(gridrefine.and.(errT.gt.epsiter*2d0)) then
			Rmax=RmaxRefine
			do i=1,10
				call RegridR(Rdes,Rmax)
				radtau=radialtau(0.55d0,tau,ntau1,D%nTheta-1)
				if(ntau1.gt.3) goto 10
			enddo
10			continue
		endif
		if (use_topac) call Topac()
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
c		write(*,'("Error on the density structure:     ",f5.1," sigma")') errRho
c		write(9,'("Error on the density structure:     ",f5.1," sigma")') errRho
		write(*,'("Error on the energy structure:      ",f5.1," sigma")') errE
		write(9,'("Error on the energy structure:      ",f5.1," sigma")') errE
		write(*,'("Error on the temperature structure: ",f5.1," sigma")') errT
		write(9,'("Error on the temperature structure: ",f5.1," sigma")') errT
		lastiter=.false.
		if(((errT.lt.epsiter.and.finaliter).and..not.(tdes_iter.and.niter.le.(nBW+2))).or.niter.gt.maxiter) then
			lastiter=.true.
			goto 1
		endif
		if(errT.lt.epsiter) then
			finaliter=.true.
		else
			finaliter=.false.
		endif
		call InitRandomWalk()
		goto 1
	else
!c -------------------------------------------------------------
!c Make the observations
!c -------------------------------------------------------------
		if(exportProDiMo.and.Nphot.gt.0) call DoExportProdimo()
		if(use_obs_TMC) then
		do j=1,D%nTheta-1
			do i=1,D%nR-1
				C(i,j)%T=C(i,j)%TMC
			enddo
		enddo
		endif
		do i=1,nobs
			call Observe(tel(i))
		enddo
	endif


	call cpu_time(stoptime)
	tottime=stoptime-starttime0

	write(*,'("Iterations for disk structure:",i11)') niter-1
	write(9,'("Iterations for disk structure:",i11)') niter-1
	write(*,'("Total runtime:                ",f11.4," s")') tottime
	write(9,'("Total runtime:                ",f11.4," s")') tottime
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	call date_and_time(date,time)
	write(*,'("Date: ",a2,"/",a2,"/",a4)') date(7:8),date(5:6),date(1:4)
	write(*,'("Time: ",a2,":",a2,":",a6)') time(1:2),time(3:4),time(5:10)
	write(9,'("Date: ",a2,"/",a2,"/",a4)') date(7:8),date(5:6),date(1:4)
	write(9,'("Time: ",a2,":",a2,":",a6)') time(1:2),time(3:4),time(5:10)
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	close(unit=9)

	call FreeAllMem()
	if(prevalloc) then
		do i=1,D%nR
		do j=1,D%nTheta
			deallocate(Cprev(i,j)%w)
			deallocate(Cprev(i,j)%w0)
			deallocate(Cprev(i,j)%wopac)
			if(.not.tcontact.or.tdes_iter) then
				deallocate(Cprev(i,j)%TP)
				deallocate(Cprev(i,j)%EJvP)
			endif
			if(use_qhp) then
				deallocate(Cprev(i,j)%QHP)
				deallocate(Cprev(i,j)%tdistr)
				deallocate(Cprev(i,j)%Tqhp)
				deallocate(Cprev(i,j)%EJvQHP)
			endif
			if(computeLRF) then
				deallocate(Cprev(i,j)%LRF)
				deallocate(Cprev(i,j)%nLRF)
			endif
		enddo
		enddo
		deallocate(Cprev)
	endif
	if(coralloc) then
		deallocate(icor)
		deallocate(jcor)
	endif
	do i=1,nobs
		if(tel(i)%kind.eq.'VISIBILITY'.or.tel(i)%kind.eq.'BASEVIS') then
			deallocate(tel(i)%b)
			deallocate(tel(i)%theta)
		endif
	enddo
	if(.not.tcontact) deallocate(denstempfileP)
	deallocate(DustMass)

	end
	
