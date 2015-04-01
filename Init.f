!c History:
!c MK 20071003:	introducing 'tfit' keyword for exponential Tevp dependence.
!c 		Find all changes with search string 'Mihkelexp'
!c-----------------------------------------------------------------------
!c Initialize the disk using input file.
!c Set are:
!c - The wavelength grid.
!c - The spatial grid.
!c - The opcaities and phase functions.
!c - The composition in the disk.
!c - The blackbodies are tabulated.
!c - All arrays are allocated.
!c - The planck mean opacities are computed
!c-----------------------------------------------------------------------
	subroutine initialize(input,Nphot,NphotFirst,NFirst,iter0)
	use Parameters
	use InputOutput
	IMPLICIT NONE
	integer i,j,k,ii,jj,scale_R,ia,Nphot,iter,iter0,NphotFirst,NFirst,number_invalid,njj
	real*8,allocatable :: w(:),ww(:,:),spec(:),dBB(:),rtemp(:)
	real*8 T,Planck,Vtot,clight,MassTot,thet,tot,Tmax_R,RTmax,stretch,Kext,Kabs
	real*8 rd,zd,f1,f2,hr,r,z,lam1,lam2,warg(MAXPART),scale,Luminosity,MassTot0
	real*8 powslope(MAXPART),powrad0(MAXPART),tau550,Kext550,wl1,wl2,vexp,vexp1,vexp2,bexp
	real*8 gap1(100),gap2(100),gap(100),Rdes,Rmin,zlam1,zlam2,dummy
	real*8 TdesA(MAXPART),TdesB(MAXPART),int1,int2,f_weight_backup,theta
	real*8 gapshape(100),gaproundpow(100)
	real*8 asym(MAXPART),asym2(MAXPART),wasym2(MAXPART),Pmax(MAXPART),tdes_fast(MAXPART),powinner(MAXPART)
	integer npow,powclose(MAXPART),powfar(MAXPART),ilam1,ilam2,ngap,nzlam,nr,nt
	character*500 particlefile,gridfile,input,partarg(MAXPART),thetagridfile
	character*500 densfile,opacityfile,surfdensfile,output,compositionfile,gasdensfile
	character*500 scalesh,shscalefile,starfile,radfile,settlefile(MAXPART)
	character*20 denstype,abuntype,lamtype,scattype,startype,material(MAXPART),startiter
	character*20 shtype(MAXPART)
	character*20 gaproundtype(MAXPART),roundtype(MAXPART)
	character*1000 line
	character*500 key,value,file,keyzone
	logical truefalse,opacity_set,arg_abun,force_vert_gf(MAXPART)
	logical mdustscale,settle(MAXPART),trace(MAXPART),denscomposition
	integer coupledabun(MAXPART),coupledfrac(MAXPART),info,nrhs,nl,parttype(MAXPART),nRfix,iopac
	real*8 frac(MAXPART),f,phi,shscalevalue,part_shscale(MAXPART),shpow,minrad(MAXPART),maxrad(MAXPART)
	real*8 RoundOff,roundwidth(MAXPART),roundpow(MAXPART),roundpeak(MAXPART) !Gijsexp
	real*8 powmix(MAXPART),radmix(MAXPART),DiffCoeff,Tmix(MAXPART),Rfix(100),rimscale,rimwidth
	real*8 rgrain(MAXPART),rgrain_edges(101),rhograin(MAXPART) ! Gijsexp
	real*8,allocatable :: abunA(:,:),abunB(:),surfacedens(:),F11(:,:)
	real*8,allocatable :: wfunc(:),g(:),waver(:),DiffC(:),zonedens(:,:,:,:)
	integer,allocatable :: IPIV(:)
	real*8 psettR0,psettpow,psettphi0,KappaGas,shaperad(MAXPART),mu,mu0,rho,minR2,maxR2
	real*8 MeixA,MeixB,MeixC,MeixD,MeixE,MeixF,MeixG,MeixRsw,timeshift,MeixRin
	real*8 radtau,tau,reprocess,tot2,mrn_index0,dens1,dens2,T_IRF,F_IRF,wedgeopen
	integer ntau1,NUV,N1UV,iz,mrn_ngrains0
	type(DiskZone) ZoneTemp(10) ! maximum of 10 Zones
	real*8 computepart_amin(MAXPART),computepart_amax(MAXPART),computepart_apow(MAXPART),computepart_fmax(MAXPART)
	real*8 computepart_porosity(MAXPART),mrn_tmp_rmin,mrn_tmp_rmax,mrn_tmp_index,maxtauV
	real*8 computepart_abun(MAXPART,MAXPART),computepart_norm_abun(MAXPART),computepart_fcarbon(MAXPART)
	real*8,allocatable :: computepart_T(:,:)
	integer computepart_ngrains(MAXPART),computepart_nT(MAXPART),computepart_nsubgrains(MAXPART)
	logical computepart_blend(MAXPART)
	character*500,allocatable :: computepart_Tfile(:,:)
	character*20 computepart_standard(MAXPART)
	integer nused,npert
	real*8 fpert,mintheta(MAXPART),maxtheta(MAXPART)

	allocate(computepart_T(MAXPART,50))
	allocate(computepart_Tfile(MAXPART,50))

	startype='PLANCK'
	D%Tstar=10000d0
	D%Rstar=2.1d0
	D%Mstar=2.5d0
	D%logg=4d0
	D%distance=100d0
	D%Av=0d0
	adjustAv=.false.
	
	D%Tstar2=-1d0
	D%Rstar2=-1d0

	D%Rin=0.5d0
	D%Rout=1000d0
	D%nR=70
	D%nTheta=60
	nspan=7
	nlev=7
	ntspan=0
	ntlev=0
	nRfix=0

	denstype='POW'
	D%denspow=1d0
	D%denspow2=1d0
	D%Rpow2=1d10
	D%Rpow3=0d0		! Gijsexp
	D%Rexp=1d10		! Gijsexp
	D%gamma_exp=-1d0
	D%Mtot=1d-15
	mdustscale=.false.
	scale_R=0
	Tmax_R=1500d0
	D%shpow=1.3
	D%sh1AU=0.03
	tau550=100d0
	reprocess=-1d0
	D%Mdot=1d-8		!Msun/yr
	vexp1=10d0		!km/s
	vexp2=-10d0		!don't use vexp2
	bexp=0.5d0		!beta velocity law
	D%Minfall=-1d0		!defaults to the value of Mdot

	struct_iter=.false.
	raditer=.false.
	gravstable=.false.
	reducemdot=.false.
	epsiter=3d0
	maxiter=10
	niter0=1000

	lamtype(1:3)='LOG'
	nlam=100
	lam1=0.1
	lam2=10000d0
	zlam1=5d0
	zlam2=40d0
	nzlam=0

	nlamHR=0

	NphotFirst=100000
	NFirst=-1

	NphotUV=250000
	maxlamUV=0.3d0
	UV_PAH=.true.
	PAHion=.true.

	scattype(1:9)='ISOTROPIC'
	scat_how=1
	scattering=.true.
	storescatt=.true.
	maxinteract=1000000
	NphotDiffuse=25
	dTDiffuse=1d0
	use_obs_TMC=.true.
	opacity_set=.false.
	overflow=.false.
	tcontact=.true.
	tdes_iter=.false.

!c Default values are those of olivine
	TdesA(1:MAXPART)=4.4491649
	TdesB(1:MAXPART)=0.35676055
	force_vert_gf(1:MAXPART)=.true.
	tdes_fast(1:MAXPART)=0d0

	arg_abun=.false.
	warg=1d0
	ngrains=0
	npow=0

	powslope(1:MAXPART)=1d0
	powrad0(1:MAXPART)=1d0
	powinner(1:MAXPART)=1d0	! Gijsexp

	shell1D=.false.
	FLD=.false.
	scset=.false.		! Gijsexp
 	scsetsave=.true.	! Gijsexp
 	scseteq=.true.		! Gijsexp
	mpset=.false.		! Gijsexp
	mpstr=.false.		! Gijsexp
	
	fixmpset=.false.

	mrn=.false.		! Gijsexp: calculate grain size distribution 
	mrn_index=3.5d0		! Gijsexp
	mrn_rmin=1d-6		! Gijsexp: 0.01 micron
	mrn_rmax=1d-1       ! Gijsexp: 1 mm
	mrn_ngrains=0		! Gijsexp: include grains up to mrn_ngrains (if >0)

	gsd=.false.		! Gijsexp: calculate grain size distribution 
	gsd_full=.false.	! Gijsexp
	gsd_plot=.false.	! Gijsexp: write excessive output files
	gsd_rmin=2.5d-6		! Gijsexp
	gsd_rmax=1d1		! Gijsexp
	gsd_xi=1.8d0		! Gijsexp: fragmentation slope
	gsd_vfrag=100d0		! Gijsexp: fragmentation velocity in cm/s
	gsd_diag=-1		! Gijsexp

	deadzone=.false.	! Gijsexp
	deadcolumn=50d0		! Gijsexp: gas column in gr/cm^2
	deadtemp=1d3		! Gijsexp
	deadalpha=1d-5		! Gijsexp: alpha in deadzone (near midplane)

	ngap=0
	gap(1:100)=0.1d0
	gapshape(1:100)=0d0
	gaproundtype(1:100)=' '   ! 
	gaproundpow(1:100)=0d0    ! r^p, NOT r^-p like SDP
	gap1(1:100)=0d0
	gap2(1:100)=0d0
	
	coupledabun(1:MAXPART)=0
	coupledfrac(1:MAXPART)=0
	frac(1:MAXPART)=0d0
	RNDW=.false.
	factRW=10d0
	scalesh='NO'
	
	wedgeopen=5d0

	minrad(1:MAXPART)=0d0
	maxrad(1:MAXPART)=1d50
	mintheta(1:MAXPART)=-360d0
	maxtheta(1:MAXPART)=360d0
	shaperad(1:MAXPART)=0d0
	roundtype(1:MAXPART)=' '	! Gijsexp: round off rim at minrad 
	roundwidth(1:MAXPART)=0d0	! Gijsexp
	roundpow(1:MAXPART)=10d0	! Gijsexp
	roundpeak(1:MAXPART)=0d0	! Gijsexp
	settle(1:MAXPART)=.false.
	settlefile(1:MAXPART)=' '
	part_shscale(1:MAXPART)=1d0
	shtype(1:MAXPART)='DISK'
	gridrefine=.true.
	thgridrefine=.true.
	nruns=1

	psettR0=2d10
	psettpow=1d0
	psettphi0=1d0

	compositionfile=' '
	gasdensfile=' '
	radfile=' '

	forcediff=.false.

	etrace=.false.
	iTD(0:TMAX)=0

	powmix(1:MAXPART)=0.65
	radmix(1:MAXPART)=1d0
	Tmix(1:MAXPART)=real(TMAX)*dT*2d0

	f_weight=0.5d0
	NphotFinal=0
	viscous=.false.
	fastviscous=.false.
	g2d_heat=.false.	! don't compute the heating of the dust by the gas by default
	
	forcefirst=.false.	! do not force the first interaction in the radiative transport
	
	NsigDiskstructure=1
	
	nBW=-1
	use_qhp=.false.
	computeLRF=.false.
	TdesQHP=-1d0
	
	use_topac=.false.
	
	D%mu0max=2d0

	multiwav=.false.
			
	rimscale=0d0
	rimwidth=3d0
	
	alphavis=0.01
	alphavispow=0d0
	getalpha=.false.        ! Gijsexp retrieve alpha for current surface density
	prandtl=-1		! relates alphaturb to alphavis (if gt 0)
	alphaturb=1d-4		! Gijsexp
	qturb=0.5		! Gijsexp
	lifetime=1d8		! Gijsexp  (large, equilibrium dust settling)
	gas2dust=100d0		! Gijsexp
	thinparticle=-2	! Gijsexp (composition in optically thin layers)
	
	D%PA=0d0
	D%IA=45d0
	
	startiter=' '
	
	tracestar=.true.
	traceemis=.true.
	tracescat=.true.
	tracegas=.false.
	trace(1:MAXPART)=.true.
	
	asym(1:MAXPART)=2d0
	asym2(1:MAXPART)=2d0
	wasym2(1:MAXPART)=0d0
	Pmax(1:MAXPART)=1d0
	
	radpress=.false.

	convection=.false.
	
	nplanets=0
	
	dimstar=1d0
	
	nRWinter=0
	
	computeTgas=.false.
	useTgas=.false.
	inner_gas=.false.
	Rinner_gas=1d0
	
	nspike=0		!number of angles made isotropic
	
	outfluxcontr=.false.	! output contributions of the flux when observing
	topac_interpol=.true.	! interpolate temperature dependent opacities
	
	MeixA=8d0
	MeixB=2d0
	MeixC=2d0
	MeixD=1d0
	MeixE=3d0
	MeixF=1d0
	MeixG=0d0
	MeixRsw=1800d0
	MeixRin=-1d0
	timeshift=0d0
	
	Tsmooth=.false.
	
	outputfits=.false.
	
	exportProDiMo=.false.
	runProDiMo=.false.
	ProDiModir='filesProDiMo'
	maxthetaProDiMo=75d0
	prodimo1zone=.false.
	exportFLiTs=.false.
	runscript=.false.
	scriptname='postprocessing.sh'

	computepart_nT=0
	computepart_amin=-1d0
	computepart_amax=-1d0
	computepart_apow=1d5
	computepart_fmax=0d0
	computepart_blend=.true.
	computepart_porosity=0.25d0
	computepart_abun=-1d0
	computepart_ngrains=1
	computepart_nsubgrains=1
	computepart_norm_abun=-1d0
	computepart_standard(1:MAXPART)='FILE'
	computepart_fcarbon=0.2
	
	maxruntime=-1
	emptylower=.false.
	
	particledir=outdir
	
	multicore=.true.

	NMAX_CONVOLUTION=8000
	
	denscomposition=.false.
	
	f_ne=1d0
	qhp_solver=0
	
	UVdes=.false.
	gammaUVdes=0d0

	lnkloglog=.false.
	
	fastobs=.false.
	
	abun_in_name=0
	
	pfstop=0.25d0
	
c	Initialize the 10 temp zones with defaults
	do i=1,10
		ZoneTemp(i)%fix_struct=.false.
		ZoneTemp(i)%sizedis=.false.
		allocate(ZoneTemp(i)%abun(MAXPART))
		allocate(ZoneTemp(i)%inc_grain(MAXPART))
		ZoneTemp(i)%inc_grain(1:MAXPART)=.true.
		ZoneTemp(i)%abun(1:MAXPART)=1d0
		ZoneTemp(i)%Rsh=1d0
		ZoneTemp(i)%sh=0.01d0
		ZoneTemp(i)%shpow=1.1
		ZoneTemp(i)%Rexp=1d200
		ZoneTemp(i)%a_min=-1d0
		ZoneTemp(i)%a_max=-1d0
		ZoneTemp(i)%a_pow=1d5
		ZoneTemp(i)%pertA=0d0
		ZoneTemp(i)%pertR=0d0
		ZoneTemp(i)%pertR0=0d0
		ZoneTemp(i)%gamma_exp=-1d0
		ZoneTemp(i)%maxtauV=-1d0
		ZoneTemp(i)%roundtype='NONE'
		ZoneTemp(i)%roundwidth=0d0
		ZoneTemp(i)%roundindex=0.30
		ZoneTemp(i)%roundscalemin=1d-60
		ZoneTemp(i)%fPAH=0d0
		ZoneTemp(i)%pertA=0d0
		ZoneTemp(i)%pertR=10d0
		ZoneTemp(i)%pertR0=0d0
		ZoneTemp(i)%Mconnect=0
		ZoneTemp(i)%Sconnect=0
		ZoneTemp(i)%Rconnect=-1d0
		ZoneTemp(i)%alphaturb=-1d0
	enddo
	nzones=0
	
c	Interstellar Radiation Field (IRF)
	T_BG=2.7d0
	T_IRF=20000d0
	T_ISMdust=2.7d0
	F_IRF=1d0
	use_IRF=.false.
	bg_correct=.false.
	
	ntau1_lam=1		! allways one extra tau=1 surface
	do i=1,100
	   tau1_lam(i)=0.55
	enddo

	do i=1,100
		material(i)='UNKNOWN'
		rgrain(i)=1d-5 ! Gijsexp, default grain size in cm, 0.1 micron
		rhograin(i)=3.0	! Gijsexp
	enddo

c	call system("rm -f " // trim(outdir) // "/prodimo_extra.in")

	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,'("Input file not found !")')
		write(*,'("using defaults")')
		write(9,'("Input file not found !")')
		write(9,'("using defaults")')
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		goto 20
	endif
	open(unit=20,file=input,RECL=1000)

	write(line,'("cp ",a," ",a,"input.dat")') trim(input),trim(outdir)
	call system(line)
	write(line,'(a,"input.dat")') trim(outdir)
	if(Nphot.gt.0) then
		open(unit=21,file=line,RECL=1000,access='APPEND')
		write(21,'("*** command line keywords ***")')
	endif

	ia=1
10	if(ia.eq.1) then
		call ignorestar(20)
		read(20,'(a1000)',end=20,err=20) line
		goto 30
	endif

20	call get_command_argument(2+ia,line)
	if(line(1:2).eq.'-s') then
		ia=ia+1
		call get_command_argument(2+ia,line)
		print*,line(1:len_trim(line))
		if(Nphot.gt.0) write(21,'(a)') trim(line)
		ia=ia+1
		goto 30
	endif
	if(line(1:2).eq.'-p') nplanets=nplanets+1
	if(line.eq.' ') goto 40
	ia=ia+1
	goto 20

30	continue

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
	
	call Keywordlist(key)
	
	if(key.eq.'startype') startype=value
	if(key.eq.'starfile') starfile=value
	if(key.eq.'tstar') read(value,*) D%Tstar
	if(key.eq.'logg') read(value,*) D%logg
	if(key.eq.'rstar') read(value,*) D%Rstar
	if(key.eq.'mstar') read(value,*) D%Mstar
	if(key.eq.'lstar') then
		read(value,*) D%Lstar
		D%Lstar=D%Lstar*Luminosity(5778d0,Rsun)
	endif
	if(key.eq.'distance') read(value,*) D%distance
	if(key.eq.'av') read(value,*) D%Av
	if(key.eq.'adjustav') read(value,*) adjustAv
	if(key.eq.'posangle') read(value,*) D%PA
	if(key.eq.'incangle') read(value,*) D%IA

	if(key.eq.'tstar2') read(value,*) D%Tstar2
	if(key.eq.'rstar2') read(value,*) D%Rstar2

	if(key.eq.'rin') read(value,*) D%Rin
	if(key.eq.'rout') read(value,*) D%Rout
	if(key.eq.'nrad') read(value,*) D%nR
	if(key.eq.'ntheta') read(value,*) D%nTheta
	if(key.eq.'nspan') read(value,*) nspan
	if(key.eq.'nlev') read(value,*) nlev
	if(key.eq.'ntspan') read(value,*) ntspan
	if(key.eq.'ntlev') read(value,*) ntlev
	if(key.eq.'gridrefine') read(value,*) gridrefine
	if(key.eq.'denstype') read(value,*) denstype
	if(key.eq.'radfile') then
		radfile=value
	endif
	if(key.eq.'densfile') then
		densfile=value
	endif
	if(key.eq.'gasdensfile') then
		gasdensfile=value
	endif
	if(key.eq.'denspow') read(value,*) D%denspow
	if(key.eq.'denspow2') read(value,*) D%denspow2
	if(key.eq.'rpow2') read(value,*) D%Rpow2
	if(key.eq.'rpow3') read(value,*) D%Rpow3
	if(key.eq.'rexp') read(value,*) D%Rexp
	if(key.eq.'gamma_exp') read(value,*) D%gamma_exp
	if(key.eq.'shpow') read(value,*) D%shpow
	if(key.eq.'sh1au') read(value,*) D%sh1AU
	if(key.eq.'mdust') then
		read(value,*) D%Mtot
		mdustscale=.true.
	endif
	if(key.eq.'tau550') read(value,*) tau550
	if(key.eq.'reprocess') read(value,*) reprocess
	if(key.eq.'mdot') read(value,*) D%Mdot
	if(key.eq.'minfall') read(value,*) D%Minfall
	if(key.eq.'mu0max') read(value,*) D%mu0max
	if(key.eq.'vexp'.or.key.eq.'vexp1') read(value,*) vexp1
	if(key.eq.'vexp2') read(value,*) vexp2
	if(key.eq.'bexp') read(value,*) bexp
	if(key.eq.'scaledisk') read(value,*) scale_R
	if(key.eq.'tmax') read(value,*) Tmax_R
	if(key.eq.'forcefirst') read(value,*) forcefirst

	if(key.eq.'iter') read(value,*) struct_iter
	if(key.eq.'raditer') read(value,*) raditer
	if(key.eq.'gravstable') read(value,*) gravstable
	if(key.eq.'reducemdot') read(value,*) reducemdot
	if(key.eq.'epsiter') read(value,*) epsiter
	if(key.eq.'maxiter') read(value,*) maxiter
	if(key.eq.'niter0') read(value,*) niter0
	if(key.eq.'startiter') read(value,*) startiter
	if(key.eq.'abun_in_name') read(value,*) abun_in_name
	if(key.eq.'lamgrid') then
		gridfile=value
		lamtype='FILE'
	endif
	if(key.eq.'lam1') then
		read(value,*) lam1
		lamtype='LOG'
	endif
	if(key.eq.'lam2') then
		read(value,*) lam2
		lamtype='LOG'
	endif
	if(key.eq.'nlam') then
		read(value,*) nlam
		lamtype='LOG'
	endif
	if(key.eq.'zlam1') read(value,*) zlam1
	if(key.eq.'zlam2') read(value,*) zlam2
	if(key.eq.'nzlam') read(value,*) nzlam
c	if(key.eq.'opacity') then
c		read(value,*) opacityfile
c		opacity_set=.true.
c		opacity_file_set=.true.
c	endif
	if(key.eq.'scattype') then
		read(value,*) scattype
		if(scattype.eq.'NONE') then
			scat_how=0
		else if(scattype.eq.'ISOTROPIC') then
			scat_how=1
		else if(scattype.eq.'FULL') then
			scat_how=2
		else if(scattype.eq.'RAYLEIGH') then
			scat_how=2
		else
			write(*,'("Unrecognized scattype")')
			write(9,'("Unrecognized scattype")')
			stop
		endif
	endif
	if(key.eq.'tcontact') read(value,*) tcontact
	if(key(1:5).eq.'tdesa') then
		read(key(6:len_trim(key)),*) i
		if(i.le.100) read(value,*) TdesA(i)
	endif
	if(key(1:5).eq.'tdesb') then
		read(key(6:len_trim(key)),*) i
		if(i.le.100) read(value,*) TdesB(i)
	endif
	if(key(1:7).eq.'forcegf') then
		read(key(8:len_trim(key)),*) i
		if(i.le.100) read(value,*) force_vert_gf(i)
	endif
	if(key(1:8).eq.'tdesfast') then
		read(key(9:len_trim(key)),*) i
		if(i.le.100) read(value,*) tdes_fast(i)
	endif

	if(key.eq.'tdesiter') read(value,*) tdes_iter
	if(key.eq.'1d') read(value,*) shell1D
	if(key.eq.'halo') read(value,*) haloswitch

	if(key.eq.'fweight') read(value,*) f_weight

	if(key.eq.'storescatt') read(value,*) storescatt
	if(key.eq.'maxinteract') then
		read(value,*) maxinteract
		if(maxinteract.ne.0) overflow=.true.
	endif
	if(key.eq.'randomwalk') read(value,*) RNDW
	if(key.eq.'etrace') read(value,*) etrace
	if(key.eq.'factrw') read(value,*) factRW
	if(key.eq.'nphotdiffuse') read(value,*) NphotDiffuse
	if(key.eq.'dtdiffuse') read(value,*) dTDiffuse
	if(key.eq.'fld') read(value,*) FLD

	if(key.eq.'nphotfinal') read(value,*) NphotFinal
	if(key.eq.'nphotfirst') read(value,*) NphotFirst
	if(key.eq.'nfirst') read(value,*) NFirst
	if(key.eq.'nphotuv') read(value,*) NphotUV
	if(key.eq.'maxlamuv') read(value,*) maxlamUV
	if(key.eq.'uv_pah') read(value,*) UV_PAH
	if(key.eq.'pahion') read(value,*) PAHion

	if(key.eq.'idum') read(value,*) idum

	if(key.eq.'innergas') read(value,*) inner_gas
	if(key.eq.'rinnergas') read(value,*) Rinner_gas		!inner radius of the gas disk in Stellar radii

	if(key.eq.'obstmc') read(value,*) use_obs_TMC
	if(key.eq.'fastobs') read(value,*) fastobs
	if(key.eq.'tracestar') read(value,*) tracestar
	if(key.eq.'dimstar') read(value,*) dimstar
	if(key.eq.'traceemis') read(value,*) traceemis
	if(key.eq.'tracescat') read(value,*) tracescat
	if(key.eq.'tracegas') read(value,*) tracegas
	if(key(1:5).eq.'trace'.and.key.ne.'tracestar'.and.key.ne.'traceemis'.and.key.ne.'tracescat'.and.key.ne.'tracegas') then
		read(key(6:len_trim(key)),*) i
		if(i.le.100) read(value,*) trace(i)
	endif

	if(key.eq.'radpress') read(value,*) radpress

	if(key.eq.'nmax_conv') read(value,*) NMAX_CONVOLUTION
	
	if(key.eq.'fluxcontr') read(value,*) outfluxcontr
	if(key.eq.'topac_interpol') read(value,*) topac_interpol	

c keyword abundances set (maximum number is 99)
	if(key(1:4).eq.'abun'.and.key.ne.'abun_in_name') then
		arg_abun=.true.
		read(key(5:len_trim(key)),*) i
		if(i.le.100) then
c			if(iregion.eq.0) then
				read(value,*) warg(i)
c			else
c				read(value,*) wregion(iregion,i)
c			endif
		endif
	endif
	if(key(1:3).eq.'mat') then
		arg_abun=.true.
		read(key(4:len_trim(key)),*) i
		if(i.le.100) material(i)=value
	endif
	if(key(1:4).eq.'part') then
		opacity_set=.true.
		read(key(5:len_trim(key)),*) i
		if(i.le.100) then
			partarg(i)=value
			if(i.gt.ngrains) ngrains=i
			parttype(i)=1
		endif
	endif
	if(key(1:8).eq.'fitspart') then
		opacity_set=.true.
		read(key(9:len_trim(key)),*) i
		if(i.le.100) then
			partarg(i)=value
			if(i.gt.ngrains) ngrains=i
			parttype(i)=6
		endif
	endif
	if(key(1:4).eq.'opac'.and.key.ne.'opacity') then
		opacity_set=.true.
		read(key(5:len_trim(key)),*) i
		if(i.le.100) then
			partarg(i)=value
			if(i.gt.ngrains) ngrains=i
			parttype(i)=2
		endif
	endif
	if(key(1:6).eq.'aggmix') then
		opacity_set=.true.
		read(key(7:len_trim(key)),*) i
		if(i.le.100) then
			partarg(i)=value
			if(i.gt.ngrains) ngrains=i
			parttype(i)=3
		endif
	endif
	if(key(1:5).eq.'topac'.and.key.ne.'topac_interpol') then
		opacity_set=.true.
		read(key(6:len_trim(key)),*) i
		if(i.le.100) then
			partarg(i)=value
			if(i.gt.ngrains) ngrains=i
			parttype(i)=5
		endif
		use_topac=.true.
	endif

	if(key.eq.'dirparticle') write(particledir,'(a,"/")') trim(value)

	if(key(1:11).eq.'computepart') then
		read(key(12:index(key,":")-1),*) i
		opacity_set=.true.
		if(i.gt.ngrains) ngrains=i
		parttype(i)=7
		write(keyzone,'(a)') key(index(key,":")+1:len_trim(key))
		if(keyzone.eq.'file') then
			partarg(i)=value
		else if(keyzone(1:5).eq.'tfile') then
			computepart_nT(i)=computepart_nT(i)+1
			computepart_Tfile(i,computepart_nT(i))=value
			read(keyzone(6:len_trim(keyzone)),*) computepart_T(i,computepart_nT(i))
			use_topac=.true.
		else if(keyzone.eq.'amin') then
			read(value,*) computepart_amin(i)
		else if(keyzone.eq.'amax') then
			read(value,*) computepart_amax(i)
		else if(keyzone.eq.'apow') then
			read(value,*) computepart_apow(i)
		else if(keyzone.eq.'fmax') then
			read(value,*) computepart_fmax(i)
		else if(keyzone.eq.'ngrains') then
			read(value,*) computepart_ngrains(i)
		else if(keyzone.eq.'nsubgrains') then
			read(value,*) computepart_nsubgrains(i)
		else if(keyzone.eq.'blend') then
			read(value,*) computepart_blend(i)
		else if(keyzone.eq.'porosity') then
			read(value,*) computepart_porosity(i)
		else if(keyzone(1:4).eq.'abun') then
			read(keyzone(5:len_trim(keyzone)),*) j
			read(value,*) computepart_abun(i,j)
		else if(keyzone.eq.'norm_abun') then
			read(value,*) computepart_norm_abun(i)
		else if(keyzone.eq.'standard') then
			computepart_standard(i)=value
		else if(keyzone.eq.'fcarbon') then
			read(value,*) computepart_fcarbon(i)
		else
			write(*,'("ComputePart keyword not understood:",a)') trim(keyzone)
			write(9,'("ComputePart keyword not understood:",a)') trim(keyzone)
			stop
		endif
	endif

	if(key(1:10).eq.'computepah') then
		read(key(11:index(key,":")-1),*) i
		opacity_set=.true.
		if(i.gt.ngrains) ngrains=i
		parttype(i)=8
		write(keyzone,'(a)') key(index(key,":")+1:len_trim(key))
		if(keyzone.eq.'amin') then
			read(value,*) computepart_amin(i)
		else if(keyzone.eq.'amax') then
			read(value,*) computepart_amax(i)
		else if(keyzone.eq.'apow') then
			read(value,*) computepart_apow(i)
		else if(keyzone.eq.'ngrains') then
			read(value,*) computepart_ngrains(i)
		else
			write(*,'("ComputePAH keyword not understood:",a)') trim(keyzone)
			write(9,'("ComputePAH keyword not understood:",a)') trim(keyzone)
			stop
		endif
		use_qhp=.true.
	endif


	if(key(1:6).eq.'powmix') then
		read(key(7:len_trim(key)),*) i
		if(i.le.100) read(value,*) powmix(i)
	endif
	if(key(1:6).eq.'radmix') then
		read(key(7:len_trim(key)),*) i
		if(i.le.100) read(value,*) radmix(i)
		Tmix(i)=real(TMAX)*dT*2d0
	endif
	if(key(1:4).eq.'tmix') then
		read(key(5:len_trim(key)),*) i
		if(i.le.100) read(value,*) Tmix(i)
	endif
	if(key(1:3).eq.'qhp') then
		opacity_set=.true.
		read(key(4:len_trim(key)),*) i
		if(i.le.100) then
			partarg(i)=value
			if(i.gt.ngrains) ngrains=i
			parttype(i)=4
		endif
		use_qhp=.true.
	endif

	! gijsexp: override ngrains on command line
	if(key.eq.'ngrains') read(value,*) ngrains 

	if(key(1:6).eq.'minrad') then
		arg_abun=.true.
		read(key(7:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) minrad(i)
		endif
	endif
	if(key(1:6).eq.'maxrad') then
		arg_abun=.true.
		read(key(7:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) maxrad(i)
		endif
	endif
	if(key(1:8).eq.'mintheta') then
		arg_abun=.true.
		read(key(9:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) mintheta(i)
		endif
	endif
	if(key(1:8).eq.'maxtheta'.and.key.ne.'maxthetaprodimo') then
		arg_abun=.true.
		read(key(9:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) maxtheta(i)
		endif
	endif
	if(key(1:8).eq.'shaperad') then
		arg_abun=.true.
		read(key(9:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) shaperad(i)
		endif
	endif
	if(key(1:9).eq.'roundtype'.and.key(1:12).ne.'roundtypegap') then
		arg_abun=.true.
		read(key(10:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) roundtype(i)
		endif
	endif
	if(key(1:10).eq.'roundwidth') then
		arg_abun=.true.
		read(key(11:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) roundwidth(i)
		endif
	endif
	if(key(1:8).eq.'roundpow'.and.key(1:11).ne.'roundpowgap') then
		arg_abun=.true.
		read(key(9:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) roundpow(i)
		endif
	endif
	if(key(1:9).eq.'roundpeak') then
		arg_abun=.true.
		read(key(10:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) roundpeak(i)
		endif
	endif
        if(key(1:6).eq.'settle'.and.key(1:10).ne.'settlefile') then
		opacity_set=.true.
		read(key(7:len_trim(key)),*) i
		if(i.le.100) then
			settle(i)=.true.
			read(value,*) part_shscale(i)
		endif
	endif
	if(key(1:10).eq.'settlefile') then
		opacity_set=.true.
		read(key(11:len_trim(key)),*) i
		if(i.le.100) then
			settle(i)=.true.
			settlefile(i)=value
		endif
	endif
	if(key(1:6).eq.'shtype') then
		opacity_set=.true.
		read(key(7:len_trim(key)),*) i
		if(i.le.100) then
			settle(i)=.true.
			shtype(i)=value
		endif
	endif
	if(key(1:4).eq.'asym'.and.key(5:5).ne.'b') then
		arg_abun=.true.
		read(key(5:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) asym(i)
		endif
	endif
	if(key(1:5).eq.'asymb') then
		arg_abun=.true.
		read(key(6:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) asym2(i)
		endif
	endif
	if(key(1:6).eq.'wasymb') then
		arg_abun=.true.
		read(key(7:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) wasym2(i)
		endif
	endif
	if(key(1:4).eq.'pmax') then
		arg_abun=.true.
		read(key(5:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) Pmax(i)
		endif
	endif

	if(key(1:15).eq.'compositionfile') compositionfile=value

c powerlaw abundance gradients
	if(key(1:9).eq.'gradslope') then
		read(key(10:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) powslope(i)
			if(i.gt.npow) npow=i
		endif
	endif
	if(key(1:9).eq.'gradinner') then
		read(key(10:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) powinner(i)
			if(i.gt.npow) npow=i

			if (powinner(i).lt.0d0) powinner(i)=0d0
			if (powinner(i).gt.1d0) powinner(i)=1d0
		endif
	endif
	if(key(1:9).eq.'gradclose') then
		read(key(10:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) powclose(i)
			if(i.gt.npow) npow=i
		endif
	endif
	if(key(1:7).eq.'gradfar') then
		read(key(8:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) powfar(i)
			if(i.gt.npow) npow=i
		endif
	endif
	if(key(1:7).eq.'gradrad') then
		read(key(8:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) powrad0(i)
			if(i.gt.npow) npow=i
		endif
	endif

c gaps in the density structure
	if(key(1:8).eq.'startgap') then
		read(key(9:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) gap1(i)
			if(i.gt.ngap) ngap=i
		endif
	endif
	if(key(1:6).eq.'endgap') then
		read(key(7:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) gap2(i)
			if(i.gt.ngap) ngap=i
		endif
	endif
	if(key(1:3).eq.'gap') then
		read(key(4:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) gap(i)
			if(i.gt.ngap) ngap=i
		endif
	endif
	if(key(1:8).eq.'shapegap') then
		read(key(9:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) gapshape(i)
			if(i.gt.ngap) ngap=i
		endif
	endif
	if(key(1:12).eq.'roundtypegap') then
		read(key(13:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) gaproundtype(i)
			if(i.gt.ngap) ngap=i
		endif
	endif
	if(key(1:11).eq.'roundpowgap') then
		read(key(12:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) gaproundpow(i)
			if(i.gt.ngap) ngap=i
		endif
	endif

	if(key(1:4).eq.'rfix') then
		read(key(5:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) Rfix(i)
			if(i.gt.nRfix) nRfix=i
		endif
	endif

	if(key(1:11).eq.'coupledabun') then
		read(key(12:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) coupledabun(i)
		endif
	endif
	if(key(1:11).eq.'coupledfrac') then
		read(key(12:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) coupledfrac(i)
		endif
	endif
	if(key(1:4).eq.'frac') then
		read(key(5:len_trim(key)),*) i
		if(i.le.100) then
			read(value,*) frac(i)
		endif
	endif

	if(key.eq.'tdesqhp') read(value,*) TdesQHP

	if(key.eq.'scalesh') read(value,*) scalesh
	if(key.eq.'shscale') read(value,*) shscalevalue
	if(key.eq.'shpow') read(value,*) shpow
	if(key.eq.'shscalefile') shscalefile=value

	if(key.eq.'wedgeopen') read(value,*) wedgeopen

c keys for parameterized settling of the grains
	if(key.eq.'psettr0') read(value,*) psettr0
	if(key.eq.'psettpow') read(value,*) psettpow
	if(key.eq.'psettphi0') then
		read(value,*) psettphi0
		settle=.true.
	endif
c----------------------------------------------

	if(key.eq.'rimscale') read(value,*) rimscale
	if(key.eq.'rimwidth') read(value,*) rimwidth
	
	if(key.eq.'nruns') read(value,*) nruns

c		set the maximum runtime for the radiative transfer in seconds
	if(key.eq.'maxruntime') read(value,*) maxruntime
	
c		possibility to empty the lower half of the density space
	if(key.eq.'emptylower') read(value,*) emptylower

	if(key.eq.'viscous') read(value,*) viscous
	if(key.eq.'fastviscous') read(value,*) fastviscous
	if(key.eq.'alphaviscous') read(value,*) alphavis
	if(key.eq.'getalpha') read(value,*) getalpha
	if(key.eq.'tgas') read(value,*) computeTgas
	if(key.eq.'g2d_heat') read(value,*) g2d_heat
	
	if(key.eq.'convection') read(value,*) convection
	
	if(key.eq.'forcediff') read(value,*) forcediff

	if(key.eq.'multiwav') read(value,*) multiwav

	if(key.eq.'nbw') read(value,*) nBW
	if(key.eq.'nspike') read(value,*) nspike

	! Wavelength(s) of tau=1 surface (vertical)
	if(key(1:8).eq.'tau1_lam') then
		read(key(9:len_trim(key)),*) i
		if(i.le.100) then
		   read(value,*) tau1_lam(i)
		endif
		if (i .gt. ntau1_lam) ntau1_lam=i
	endif


C       Gijsexp, read in parameters for s.c. settling
	if(key.eq.'scset') read(value,*) scset
	if(key.eq.'mpset') read(value,*) mpset
	if(key.eq.'mpstr') read(value,*) mpstr ! isothermal vertical structure
	if(key.eq.'fixmpset') read(value,*) fixmpset
	if(key.eq.'scseteq') read(value,*) scseteq
	if(key.eq.'scsetsave') read(value,*) scsetsave
	if(key.eq.'lifetime') read(value,*) lifetime
	if(key.eq.'alphaturb') read(value,*) alphaturb
	if(key.eq.'gas2dust') read(value,*) gas2dust
	if(key.eq.'qturb') read(value,*) qturb
	if(key.eq.'thinparticle') read(value,*) thinparticle
	if(key(1:6).eq.'rgrain') then
	   read(key(7:len_trim(key)),*) i
	   if(i.le.100) then
	      read(value,*) rgrain(i)
	      rgrain(i)=rgrain(i)*1d-4
	   endif
	endif
	if(key(1:8).eq.'rhograin') then
	   read(key(9:len_trim(key)),*) i
	   if(i.le.100) then
	      read(value,*) rhograin(i)
	   endif
	endif
	
	! parameters for grain size distribution
	if(key.eq.'mrn') read(value,*) mrn	
	if(key.eq.'mrn_index'.or.key.eq.'apow') read(value,*) mrn_index
	if(key.eq.'mrn_rmin'.or.key.eq.'amin') then
	   read(value,*) mrn_rmin
	   mrn_rmin = mrn_rmin * 1d-4
	endif
	if(key.eq.'mrn_rmax'.or.key.eq.'amax') then
	   read(value,*) mrn_rmax
	   mrn_rmax = mrn_rmax * 1d-4
	endif
	if(key.eq.'mrn_ngrains') read(value,*) mrn_ngrains
	
	! parameters for grain size distribution
	if(key.eq.'gsd') read(value,*) gsd	
	if(key.eq.'gsd_full') read(value,*) gsd_full
	if(key.eq.'gsd_plot') read(value,*) gsd_plot	
	if(key.eq.'gsd_rmin') then
	   read(value,*) gsd_rmin
	   gsd_rmin = gsd_rmin * 1d-4
	endif
	if(key.eq.'gsd_rmax') then
	   read(value,*) gsd_rmax 
	   gsd_rmax = gsd_rmax * 1d-4
	endif
	if(key.eq.'gsd_xi') read(value,*) gsd_xi	
	if(key.eq.'gsd_vfrag') read(value,*) gsd_vfrag	
	if(key.eq.'gsd_diag') read(value,*) gsd_diag

	!dead zone near the midplane
	if(key.eq.'deadzone') read(value,*) deadzone
	if(key.eq.'deadcolumn') read(value,*) deadcolumn
	if(key.eq.'deadtemp') read(value,*) deadtemp
	if(key.eq.'deadalpha') read(value,*) deadalpha

	if(key.eq.'prandtl') read(value,*) prandtl

	if(key.eq.'lnkloglog') read(value,*) lnkloglog

	!have a hot/shocked radius region (hot gas)
	if(key.eq.'hotgasminrad') read(value,*) HotGasMinRad
	if(key.eq.'hotgasmaxrad') read(value,*) HotGasMaxRad
	if(key.eq.'hotgast') read(value,*) HotGasT

	!Meixner density distribution (see Meixner et al. 2002)
	if(key.eq.'meixa') read(value,*) MeixA
	if(key.eq.'meixb') read(value,*) MeixB
	if(key.eq.'meixc') read(value,*) MeixC
	if(key.eq.'meixd') read(value,*) MeixD
	if(key.eq.'meixe') read(value,*) MeixE
	if(key.eq.'meixf') read(value,*) MeixF
	if(key.eq.'meixg') read(value,*) MeixG
	if(key.eq.'meixrsw') read(value,*) MeixRsw
	if(key.eq.'meixrin') read(value,*) MeixRin
	if(key.eq.'timeshift') read(value,*) timeshift
	
	!Interstellar Radiation Field
	if(key.eq.'irf') read(value,*) use_IRF
	if(key.eq.'t_bg') read(value,*) T_BG
	if(key.eq.'t_irf') read(value,*) T_IRF
	if(key.eq.'t_ismdust') read(value,*) T_ISMdust
	if(key.eq.'f_irf') read(value,*) F_IRF
	if(key.eq.'bg_correct') read(value,*) bg_correct

	!scale the electron density with respect to the ISM value
	if(key.eq.'f_ne') read(value,*) f_ne
	if(key.eq.'solver_qhp') read(value,*) qhp_solver ! 0=Kees, 1=MC

	if(key.eq.'exportprodimo') read(value,*) exportProDiMo
	if(key.eq.'runprodimo') read(value,*) runProDiMo
	if(key.eq.'dirprodimo') write(ProDiModir,'(a,"/")') trim(value)
	if(key.eq.'prodimo1zone') read(value,*) prodimo1zone
	if(key.eq.'runscript') read(value,*) runscript
	if(key.eq.'scriptname') scriptname=value
	if(key.eq.'maxthetaprodimo') read(value,*) maxthetaProDiMo
	if(key.eq.'exportflits') read(value,*) exportFLiTs
	if(key.eq.'tsmooth') read(value,*) Tsmooth

	if(key.eq.'gammauvdes') read(value,*) gammaUVdes
	if(key.eq.'uvdes') read(value,*) UVdes
	
	if(key(1:7).eq.'prodimo'.and.key(8:8).eq.':') call writeExtraProDiMo(line(9:len_trim(key)),trim(value))

	if(key.eq.'outputfits') read(value,*) outputfits
	if(key.eq.'multicore') read(value,*) multicore

	if(key(1:4).eq.'zone') then
		read(key(5:index(key,":")-1),*) i
		if(i.gt.nzones) nzones=i
		write(keyzone,'(a)') key(index(key,":")+1:len_trim(key))
		if(keyzone.eq.'rin') then
			read(value,*) ZoneTemp(i)%Rin
		else if(keyzone.eq.'rout') then
			read(value,*) ZoneTemp(i)%Rout
		else if(keyzone.eq.'denspow') then
			read(value,*) ZoneTemp(i)%denspow
		else if(keyzone.eq.'rexp') then
			read(value,*) ZoneTemp(i)%Rexp
		else if(keyzone.eq.'gamma_exp') then
			read(value,*) ZoneTemp(i)%gamma_exp
		else if(keyzone.eq.'sh') then
			read(value,*) ZoneTemp(i)%sh
		else if(keyzone.eq.'shpow') then
			read(value,*) ZoneTemp(i)%shpow
		else if(keyzone.eq.'rsh') then
			read(value,*) ZoneTemp(i)%Rsh
		else if(keyzone.eq.'mdust') then
			read(value,*) ZoneTemp(i)%Mdust
		else if(keyzone.eq.'maxtau') then
			read(value,*) ZoneTemp(i)%maxtauV
		else if(keyzone.eq.'amin') then
			read(value,*) ZoneTemp(i)%a_min
		else if(keyzone.eq.'amax') then
			read(value,*) ZoneTemp(i)%a_max
		else if(keyzone.eq.'apow') then
			read(value,*) ZoneTemp(i)%a_pow
		else if(keyzone.eq.'fix') then
			read(value,*) ZoneTemp(i)%fix_struct
		else if(keyzone.eq.'sizedis') then
			read(value,*) ZoneTemp(i)%sizedis
		else if(keyzone.eq.'fpah') then
			read(value,*) ZoneTemp(i)%fPAH
		else if(keyzone.eq.'perta') then
			read(value,*) ZoneTemp(i)%pertA
		else if(keyzone.eq.'pertr') then
			read(value,*) ZoneTemp(i)%pertR
		else if(keyzone.eq.'pertr0') then
			read(value,*) ZoneTemp(i)%pertR0
		else if(keyzone(1:4).eq.'abun') then
			read(keyzone(5:len_trim(keyzone)),*) j
			read(value,*) ZoneTemp(i)%abun(j)
		else if(keyzone(1:4).eq.'incl') then
			read(keyzone(5:len_trim(keyzone)),*) j
			read(value,*) ZoneTemp(i)%inc_grain(j)
		else if(keyzone.eq.'roundtype') then
			read(value,*) ZoneTemp(i)%roundtype
		else if(keyzone.eq.'roundwidth') then
			read(value,*) ZoneTemp(i)%roundwidth
		else if(keyzone.eq.'roundindex') then
			read(value,*) ZoneTemp(i)%roundindex
		else if(keyzone.eq.'roundscalemin') then
			read(value,*) ZoneTemp(i)%roundscalemin
		else if(keyzone.eq.'mconnect') then
			read(value,*) ZoneTemp(i)%Mconnect
		else if(keyzone.eq.'sconnect') then
			read(value,*) ZoneTemp(i)%Sconnect
		else if(keyzone.eq.'rconnect') then
			read(value,*) ZoneTemp(i)%Rconnect
		else if(keyzone.eq.'alphaturb') then
			read(value,*) ZoneTemp(i)%alphaturb
		else
			write(*,'("Zone keyword not understood:",a)') trim(keyzone)
			write(9,'("Zone keyword not understood:",a)') trim(keyzone)
			stop
		endif
	endif

C       End 

	goto 10
40	close(unit=20)
	close(unit=21)

57	j=ngrains
	do ii=1,j
		if(computepart_amin(ii).le.0d0) computepart_amin(ii)=mrn_rmin*1d4
		if(computepart_amax(ii).le.0d0) computepart_amax(ii)=mrn_rmax*1d4
		if(computepart_apow(ii).gt.100d0) computepart_apow(ii)=mrn_index
		if((parttype(ii).eq.7.or.parttype(ii).eq.8).and.computepart_ngrains(ii).gt.1) then
			do i=ngrains,ii+1,-1
				computepart_Tfile(i+computepart_ngrains(ii)-1,:)=computepart_Tfile(i,:)
				computepart_T(i+computepart_ngrains(ii)-1,:)=computepart_T(i,:)
				computepart_nT(i+computepart_ngrains(ii)-1)=computepart_nT(i)
				computepart_norm_abun(i+computepart_ngrains(ii)-1)=computepart_norm_abun(i)
				computepart_amin(i+computepart_ngrains(ii)-1)=computepart_amin(i)
 				computepart_amax(i+computepart_ngrains(ii)-1)=computepart_amax(i)
 				computepart_apow(i+computepart_ngrains(ii)-1)=computepart_apow(i)
				computepart_fmax(i+computepart_ngrains(ii)-1)=computepart_fmax(i)
				computepart_blend(i+computepart_ngrains(ii)-1)=computepart_blend(i)
				computepart_porosity(i+computepart_ngrains(ii)-1)=computepart_porosity(i)
				computepart_ngrains(i+computepart_ngrains(ii)-1)=computepart_ngrains(i)
				computepart_nsubgrains(i+computepart_ngrains(ii)-1)=computepart_nsubgrains(i)
				computepart_abun(i+computepart_ngrains(ii)-1,:)=computepart_abun(i,:)
				computepart_standard(i+computepart_ngrains(ii)-1)=computepart_standard(i)
				computepart_fcarbon(i+computepart_ngrains(ii)-1)=computepart_fcarbon(i)
				warg(i+computepart_ngrains(ii)-1)=warg(i)
				powslope(i+computepart_ngrains(ii)-1)=powslope(i)
				powrad0(i+computepart_ngrains(ii)-1)=powrad0(i)
				TdesA(i+computepart_ngrains(ii)-1)=TdesA(i)
				TdesB(i+computepart_ngrains(ii)-1)=TdesB(i)
				material(i+computepart_ngrains(ii)-1)=material(i)
				asym(i+computepart_ngrains(ii)-1)=asym(i)
				asym2(i+computepart_ngrains(ii)-1)=asym2(i)
				wasym2(i+computepart_ngrains(ii)-1)=wasym2(i)
				Pmax(i+computepart_ngrains(ii)-1)=Pmax(i)
				tdes_fast(i+computepart_ngrains(ii)-1)=tdes_fast(i)
				powinner(i+computepart_ngrains(ii)-1)=powinner(i)
				powclose(i+computepart_ngrains(ii)-1)=powclose(i)
				powfar(i+computepart_ngrains(ii)-1)=powfar(i)
				partarg(i+computepart_ngrains(ii)-1)=partarg(i)
				settlefile(i+computepart_ngrains(ii)-1)=settlefile(i)
				material(i+computepart_ngrains(ii)-1)=material(i)
				shtype(i+computepart_ngrains(ii)-1)=shtype(i)
				settle(i+computepart_ngrains(ii)-1)=settle(i)
				trace(i+computepart_ngrains(ii)-1)=trace(i)
				coupledabun(i+computepart_ngrains(ii)-1)=coupledabun(i)
				coupledfrac(i+computepart_ngrains(ii)-1)=coupledfrac(i)
				frac(i+computepart_ngrains(ii)-1)=frac(i)
				parttype(i+computepart_ngrains(ii)-1)=parttype(i)
				part_shscale(i+computepart_ngrains(ii)-1)=part_shscale(i)
				minrad(i+computepart_ngrains(ii)-1)=minrad(i)
				maxrad(i+computepart_ngrains(ii)-1)=maxrad(i)
				mintheta(i+computepart_ngrains(ii)-1)=mintheta(i)
				maxtheta(i+computepart_ngrains(ii)-1)=maxtheta(i)
				roundwidth(i+computepart_ngrains(ii)-1)=roundwidth(i)
				roundpow(i+computepart_ngrains(ii)-1)=roundpow(i)
				roundpeak(i+computepart_ngrains(ii)-1)=roundpeak(i)
				roundtype(i+computepart_ngrains(ii)-1)=roundtype(i)
				powmix(i+computepart_ngrains(ii)-1)=powmix(i)
				radmix(i+computepart_ngrains(ii)-1)=radmix(i)
				Tmix(i+computepart_ngrains(ii)-1)=Tmix(i)
				rgrain(i+computepart_ngrains(ii)-1)=rgrain(i)
				rhograin(i+computepart_ngrains(ii)-1)=rhograin(i)
			enddo
			do i=2,computepart_ngrains(ii)
				computepart_Tfile(i+ii-1,:)=computepart_Tfile(ii,:)
				computepart_T(i+ii-1,:)=computepart_T(ii,:)
				computepart_nT(i+ii-1)=computepart_nT(ii)
				computepart_norm_abun(i+ii-1)=computepart_norm_abun(ii)
				computepart_amin(i+ii-1)=10d0**(log10(computepart_amin(ii))
     &					+log10(computepart_amax(ii)/computepart_amin(ii))
     &					*real(i-1)/real(computepart_ngrains(ii)))
				computepart_amax(i+ii-1)=10d0**(log10(computepart_amin(ii))
     &					+log10(computepart_amax(ii)/computepart_amin(ii))
     &					*real(i)/real(computepart_ngrains(ii)))
				computepart_apow(i+ii-1)=computepart_apow(ii)
				computepart_fmax(i+ii-1)=computepart_fmax(ii)
				computepart_blend(i+ii-1)=computepart_blend(ii)
				computepart_porosity(i+ii-1)=computepart_porosity(ii)
				computepart_abun(i+ii-1,:)=computepart_abun(ii,:)
				computepart_ngrains(i+ii-1)=1
				computepart_nsubgrains(i+ii-1)=computepart_nsubgrains(ii)
				computepart_standard(i+ii-1)=computepart_standard(ii)
				computepart_fcarbon(i+ii-1)=computepart_fcarbon(ii)
				rgrain(i+ii-1)=sqrt(computepart_amin(i+ii-1)*computepart_amax(i+ii-1))*1d-4
				warg(i+ii-1)=warg(ii)
				powslope(i+ii-1)=powslope(ii)
				powrad0(i+ii-1)=powrad0(ii)
				TdesA(i+ii-1)=TdesA(ii)
				TdesB(i+ii-1)=TdesB(ii)
				material(i+ii-1)=material(ii)
				asym(i+ii-1)=asym(ii)
				asym2(i+ii-1)=asym2(ii)
				wasym2(i+ii-1)=wasym2(ii)
				Pmax(i+ii-1)=Pmax(ii)
				tdes_fast(i+ii-1)=tdes_fast(ii)
				powinner(i+ii-1)=powinner(ii)
				powclose(i+ii-1)=powclose(ii)
				powfar(i+ii-1)=powfar(ii)
				partarg(i+ii-1)=partarg(ii)
				settlefile(i+ii-1)=settlefile(ii)
				material(i+ii-1)=material(ii)
				shtype(i+ii-1)=shtype(ii)
				settle(i+ii-1)=settle(ii)
				trace(i+ii-1)=trace(ii)
				coupledabun(i+ii-1)=coupledabun(ii)
				coupledfrac(i+ii-1)=coupledfrac(ii)
				frac(i+ii-1)=frac(ii)
				parttype(i+ii-1)=parttype(ii)
				part_shscale(i+ii-1)=part_shscale(ii)
				minrad(i+ii-1)=minrad(ii)
				maxrad(i+ii-1)=maxrad(ii)
				mintheta(i+ii-1)=mintheta(ii)
				maxtheta(i+ii-1)=maxtheta(ii)
				roundwidth(i+ii-1)=roundwidth(ii)
				roundpow(i+ii-1)=roundpow(ii)
				roundpeak(i+ii-1)=roundpeak(ii)
				roundtype(i+ii-1)=roundtype(ii)
				powmix(i+ii-1)=powmix(ii)
				radmix(i+ii-1)=radmix(ii)
				Tmix(i+ii-1)=Tmix(ii)
				rhograin(i+ii-1)=rhograin(ii)
			enddo
			ngrains=ngrains+computepart_ngrains(ii)-1
			if(computepart_ngrains(ii).ne.1) then
				mrn_tmp_rmin=mrn_rmin
				mrn_tmp_rmax=mrn_rmax
				mrn_tmp_index=mrn_index
				mrn_rmin=computepart_amin(ii)
				mrn_rmax=computepart_amax(ii)
				mrn_index=computepart_apow(ii)

				computepart_amax(ii)=10d0**(log10(computepart_amin(ii))
     &				+log10(computepart_amax(ii)/computepart_amin(ii))
     &				*real(1)/real(computepart_ngrains(ii)))
				rgrain(ii)=sqrt(computepart_amin(ii)*computepart_amax(ii))*1d-4

				tot=warg(ii)
				allocate(w(ngrains))
				allocate(rtemp(ngrains))
				j=0
				do i=ii,ii+computepart_ngrains(ii)-1
					j=j+1
					rtemp(j)=rgrain(i)*1d4
				enddo
				mrn_ngrains0=mrn_ngrains
				mrn_ngrains=j
				w(1:ngrains)=0d0
				do i=j+1,ngrains
					rtemp(i)=(10d0+real(i))*mrn_rmax
				enddo
				call gsd_MRN(rtemp(1:ngrains),w(1:ngrains))
				mrn_ngrains=mrn_ngrains0
				j=0
				do i=ii,ii+computepart_ngrains(ii)-1
					j=j+1
					warg(i)=tot*w(j)
				enddo
				deallocate(w)
				deallocate(rtemp)
   	  			computepart_ngrains(ii)=1
				mrn_rmin=mrn_tmp_rmin
				mrn_rmax=mrn_tmp_rmax
				mrn_index=mrn_tmp_index
			else
				computepart_amax(ii)=10d0**(log10(computepart_amin(ii))
     &				+log10(computepart_amax(ii)/computepart_amin(ii))
     &				*real(1)/real(computepart_ngrains(ii)))
				rgrain(ii)=sqrt(computepart_amin(ii)*computepart_amax(ii))*1d-4
			endif
			goto 57
		endif
	enddo

	if(nzones.ne.0) then
		D%Mtot=0d0
		write(denstype,'("ZONES")')
		allocate(Zone(nzones))
		D%Rin=ZoneTemp(1)%Rin
		D%Rout=ZoneTemp(nzones)%Rout
		do i=1,nzones
			if(ZoneTemp(i)%Rin.lt.D%Rin) D%Rin=ZoneTemp(i)%Rin
			if(ZoneTemp(i)%Rout.gt.D%Rout) D%Rout=ZoneTemp(i)%Rout
		enddo

		do i=1,nzones
			allocate(Zone(i)%abun(ngrains))
			allocate(Zone(i)%inc_grain(ngrains))
			Zone(i)=ZoneTemp(i)
			if(Zone(i)%a_min.le.0d0) Zone(i)%a_min=mrn_rmin*1d4
			if(Zone(i)%a_max.le.0d0) Zone(i)%a_max=mrn_rmax*1d4
			if(Zone(i)%a_pow.gt.100d0) Zone(i)%a_pow=mrn_index
			if(Zone(i)%gamma_exp.lt.0d0) Zone(i)%gamma_exp=2d0-Zone(i)%denspow
			if(Zone(i)%alphaturb.lt.0d0) Zone(i)%alphaturb=alphaturb
			D%Mtot=D%Mtot+Zone(i)%Mdust
			if(Zone(i)%Rin.gt.D%Rin.and.Zone(i)%Rin.lt.D%Rout.and.minval(abs(Rfix(1:nRfix)-Zone(i)%Rin)).ne.0d0) then
				nRfix=nRfix+1
				Rfix(nRfix)=Zone(i)%Rin
			endif
			if(Zone(i)%Rout.gt.D%Rin.and.Zone(i)%Rout.lt.D%Rout.and.minval(abs(Rfix(1:nRfix)-Zone(i)%Rout)).ne.0d0) then
				nRfix=nRfix+1
				Rfix(nRfix)=Zone(i)%Rout
			endif
			if(mpset.or.scset.or..not.Zone(i)%fix_struct) struct_iter=.true.
			if(Zone(i)%Rin.lt.D%Rin) D%Rin=Zone(i)%Rin
			if(Zone(i)%Rout.gt.D%Rout) D%Rout=Zone(i)%Rout
		enddo
		call ConnectZones()
	endif

	if(D%Minfall.lt.0d0) D%Minfall=D%Mdot

	if(MeixRin.le.0d0) MeixRin=D%Rin
	if(MeixD.le.0d0) MeixD=MeixE
	if(MeixE.le.0d0) MeixE=MeixD
	if(viscous) computeTgas=.true.
	if(computeTgas.or.viscous.or.denstype.eq.'PRODIMO') useTgas=.true.

	if(runProDiMo) then
		exportProDiMo=.true.
		call writeExtraProDiMo('use_MCFOST_rgrid','.false.')
		call writeExtraProDiMo('readMCFOST','.true.')
		call writeExtraProDiMo('','forProDiMo.fits.gz')
		call writeExtraProDiMo('FLiTs','.true.')
	endif
	if(exportProDiMo.or.UVdes) computeLRF=.true.
	if(fastobs) then
		if(Nphot.le.0) then
			fastobs=.false.
		else
			computeLRF=.true.
			storescatt=.false.
		endif
	endif

	if(denstype.eq.'FILE'.and.densfile.eq.compositionfile) denscomposition=.true.

	if(nphotdiffuse.eq.0) use_obs_TMC=.false.
	if(Nphot.eq.0) then
		denstype='FILE'
		if(outputfits) then
			write(densfile,'(a,"denstemp.fits.gz")') outdir(1:len_trim(outdir))
		else
			write(densfile,'(a,"denstemp.dat")') outdir(1:len_trim(outdir))
		endif
		mdustscale=.false.
		struct_iter=.false.
	endif
	if(denstype.eq.'MASSLOSS'.or.denstype.eq.'INFALL') then
		mdustscale=.false.
	endif

	if(D%gamma_exp.lt.0d0) D%gamma_exp=2d0-D%denspow


	if(viscous) forcefirst=.false. !cannot use first interaction forcing with viscous heating
	
C       Gijsexp, disable scaleheight scaling when doing the s.c. settling
	if (scset) scalesh='NONE'
	if (mpset) scalesh='NONE'	
	if (fixmpset) scalesh='NONE'	
C	End

C       Gijsexp, can now also use roundpeak
	do i=1,100
	   if (roundpeak(i).ne.0d0.and.roundwidth(i).eq.0d0) then
	      roundwidth(i)=roundpeak(i)-minrad(i)
	   endif
	enddo
C	End

	!  if endgap(i) not supplied, use startgap(i+1) or outer disk
	do i=1,100
	   if (gap2(i).eq.0d0.and.i.lt.ngap) then
	      gap2(i)= gap1(i+1)
	   elseif (gap2(i).eq.0d0.and.i.eq.ngap) then
	      gap2(i)= D%Rout
	   endif
	enddo

	!  if startgap(i) not supplied, use endgap(i-1) or inner disk
	do i=1,100
	   if (gap1(i).eq.0d0.and.i.gt.1) then
	      gap1(i)= gap2(i-1)
	   elseif (gap1(i).eq.0d0.and.i.eq.1) then
	      gap1(i)= D%Rin
	   endif
	enddo

	if(arg_abun.and.sum(coupledabun(1:ngrains)).ne.0) then
	allocate(abunA(ngrains,ngrains))
	allocate(abunB(ngrains))
	allocate(IPIV(ngrains))
	do i=1,ngrains
		do j=1,ngrains
			abunA(i,j)=0d0
		enddo
		if(coupledabun(i).eq.0) then
			f=1d0
			call MakeMatrixCoupled(abunA,coupledabun,i,i,ngrains)
			abunA(i,i)=1d0
			abunB(i)=warg(i)
		else
			abunA(i,coupledabun(i))=frac(coupledfrac(i))
			abunA(i,i)=-(1d0-frac(coupledfrac(i)))
			abunB(i)=0d0
		endif		
	enddo
	NRHS=1
	call DGESV(ngrains,NRHS,abunA,ngrains,IPIV,abunB,ngrains,INFO)
	do i=1,ngrains
		warg(i)=abunB(i)
	enddo
	deallocate(abunA)
	deallocate(abunB)
	deallocate(IPIV)
	endif

	if(timeshift.ne.0d0.and.denstype.eq.'MEIXNER') then
		D%Rin=D%Rin+timeshift*vexp1*0.210944839d0
		D%Rout=D%Rout+timeshift*vexp1*0.210944839d0
	endif

	if((D%Rstar*Rsun/AU).ge.D%Rin) then
		write(*,'("Inner radius inside the star!")')
		write(*,'("Changing inner radius")')
		write(9,'("Inner radius inside the star!")')
		write(9,'("Changing inner radius")')
		D%Rin=1.01*D%Rstar*Rsun/AU
	endif

	if(startype.eq.'PLANCK') then
	write(*,'("Stellar temperature:  ",f14.3," K")') D%Tstar
	write(*,'("Stellar radius:       ",f14.3," Rsun")') D%Rstar
	write(9,'("Stellar temperature:  ",f14.3," K")') D%Tstar
	write(9,'("Stellar radius:       ",f14.3," Rsun")') D%Rstar
	D%Lstar=Luminosity(D%Tstar,D%Rstar*Rsun)
	else if(startype.eq.'FILE') then
	write(*,'("Star file:            ",a)') starfile(1:len_trim(starfile))
	write(9,'("Star file:            ",a)') starfile(1:len_trim(starfile))
	else if(startype.eq.'KURUCZ') then
	write(*,'("Stellar temperature:  ",f14.3," K")') D%Tstar
	write(*,'("Stellar log(g):       ",f14.3)') D%logg
	write(9,'("Stellar temperature:  ",f14.3," K")') D%Tstar
	write(9,'("Stellar log(g):       ",f14.3)') D%logg
	else
	write(*,'("Error in star type")')
	write(9,'("Error in star type")')
	stop
	endif
	write(*,'("Stellar luminosity:   ",f14.3," Lsun")') D%Lstar/Luminosity(5778d0,Rsun)
	write(9,'("Stellar luminosity:   ",f14.3," Lsun")') D%Lstar/Luminosity(5778d0,Rsun)
	write(*,'("Stellar mass:         ",f14.3," Msun")') D%Mstar
	write(9,'("Stellar mass:         ",f14.3," Msun")') D%Mstar
	
	write(*,'("Distance:             ",f14.3," parsec")') D%distance
	write(9,'("Distance:             ",f14.3," parsec")') D%distance
	write(*,'("Interstellar Av:      ",f14.3)') D%Av
	write(9,'("Interstellar Av:      ",f14.3)') D%Av

	if(raditer) then
	   write(*,'("Solving radial structure")')
	   write(9,'("Solving radial structure")')
	   if(gravstable) then
	      write(*,'("Keeping disk gravitationally stable")')
	      write(9,'("Keeping disk gravitationally stable")')
	      if(reducemdot) then
		 write(*,'("Reducing accretion rate if necessary")')
		 write(9,'("Reducing accretion rate if necessary")')
	      endif
	   endif
	endif
	if(viscous) then
	   write(*,'("Using viscous heating")')
	   write(9,'("Using viscous heating")')
	endif
	if(inner_gas) then
	write(*,'("Using inner accreting gas disk")')
	write(9,'("Using inner accreting gas disk")')
	write(*,'("Disk inwards to ",f5.2," stellar radii")') Rinner_gas
	write(9,'("Disk inwards to ",f5.2," stellar radii")') Rinner_gas
	write(*,'("Mass accretion:       ",e14.3," Msun/yr")') D%Mdot
	write(9,'("Mass accretion:       ",e14.3," Msun/yr")') D%Mdot
	endif
	if(viscous.or.raditer) then
	   write(*,'("Mass accretion:       ",e14.3," Msun/yr")') D%Mdot
	   write(9,'("Mass accretion:       ",e14.3," Msun/yr")') D%Mdot
	   if(D%Rexp.lt.1d10) then 
	      write(*,'("Exponential taper at:    ",f14.3," AU")') D%Rexp
	      write(9,'("Exponential taper at:    ",f14.3," AU")') D%Rexp
	   endif
	   if(prandtl.le.0) then
	      if(alphavispow.eq.0d0) then
		 write(*,'("Alpha:                ",e14.3)') alphavis
		 write(9,'("Alpha:                ",e14.3)') alphavis
	      else
		 write(*,'("Alpha:                ",e14.3,"*(R/AU)^",f5.2)') alphavis,alphavispow
		 write(9,'("Alpha:                ",e14.3,"*(R/AU)^",f5.2)') alphavis,alphavispow
	      endif
	   else
	      write(*,'("Using alphavis=alphaturb")')
	      write(9,'("Using alphavis=alphaturb")')
	   endif		! prandtl
	endif ! viscous

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	if(denstype.ne.'FILE'.and.radfile.eq.' ') then
	write(*,'("Inner radius:         ",f14.3," AU")') D%Rin
	write(*,'("Outer radius:         ",f14.3," AU")') D%Rout
	write(*,'("Nrad:                 ",i10)') D%nR
	write(*,'("Ntheta:               ",i10)') D%nTheta
	write(*,'("Nspan (radial):       ",i10)') nspan
	if(nlev*nspan.ge.D%nR-1) then
	nlev=(D%nR-1)/nspan
	write(*,'("Lowering Nlevels to:  ",i10)') nlev
	else
	write(*,'("Nlevels (radial):     ",i10)') nlev
	endif
	write(*,'("Nspan (theta):        ",i10)') ntspan
	if(ntlev*ntspan.ge.D%nTheta-1) then
	ntlev=(D%nTheta-1)/ntspan
	write(*,'("Lowering Nlevels to:  ",i10)') ntlev
	else
	write(*,'("Nlevels (theta):      ",i10)') ntlev
	endif
	if(D%Rin*AU.lt.D%Rstar*Rsun) then
		write(*,'("Inner radius is inside the star")')
		write(*,'("Rescaling disc inner radius")')
		write(9,'("Inner radius is inside the star")')
		write(9,'("Rescaling disc inner radius")')
		D%Rin=D%Rstar*Rsun*1.0001/AU
	endif
	write(9,'("Inner radius:         ",f14.3," AU")') D%Rin
	write(9,'("Outer radius:         ",f14.3," AU")') D%Rout
	write(9,'("Nr:                   ",i10)') D%nR
	write(9,'("Ntheta:               ",i10)') D%nTheta
	write(9,'("Nspan (radial):       ",i10)') nspan
	write(9,'("Nlevels (radial):     ",i10)') nlev
	write(9,'("Nspan (theta):        ",i10)') ntspan
	write(9,'("Nlevels (theta):      ",i10)') ntlev
	else if(radfile.ne.' ') then
	write(*,'("Radial grid file:     ",a)') radfile(1:len_trim(radfile))
	write(9,'("Radial grid file:     ",a)') radfile(1:len_trim(radfile))
	endif
	do i=1,nRfix
	write(*,'("Fixed radius:         ",f14.3," AU")') Rfix(i)
	write(9,'("Fixed radius:         ",f14.3," AU")') Rfix(i)
	enddo
	if(gridrefine) then
	write(*,'("Using automatic grid refinement")')
	write(9,'("Using automatic grid refinement")')
	else
	write(*,'("No automatic grid refinement")')
	write(9,'("No automatic grid refinement")')
	endif
	if(forcefirst) then
	write(*,'("Forcing first interaction")')
	write(9,'("Forcing first interaction")')
	endif
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	if(nruns.gt.1) then
	write(*,'("Number of runs in each iteration",i4)') nruns
	write(9,'("Number of runs in each iteration",i4)') nruns
	endif
	
	if(mdustscale) then
	write(*,'("Dust mass:            ",e18.3," Msun")') D%Mtot
	write(9,'("Dust mass:            ",e18.3," Msun")') D%Mtot
	endif
	write(*,'("Density type:         ",a10)') denstype(1:len_trim(denstype))
	write(9,'("Density type:         ",a10)') denstype(1:len_trim(denstype))
	if(shell1D) then
	write(*,'("1D dust shell")')
	write(9,'("1D dust shell")')
	endif
	if(denstype.eq.'FILE') then
	write(*,'("Density file:         ",a)') densfile(1:len_trim(densfile))
	write(9,'("Density file:         ",a)') densfile(1:len_trim(densfile))
	else if(denstype.eq.'SURFFILE') then
	write(*,'("Surface density file: ",a)') densfile(1:len_trim(densfile))
	write(9,'("Surface density file: ",a)') densfile(1:len_trim(densfile))
		if(gasdensfile.ne.' ') then
			write(*,'("Gas surface density:  ",a)') gasdensfile(1:len_trim(gasdensfile))
			write(9,'("Gas surface density:  ",a)') gasdensfile(1:len_trim(gasdensfile))
		endif
	else if(denstype.eq.'SHELLFILE') then
	write(*,'("Shell density file:   ",a)') densfile(1:len_trim(densfile))
	write(9,'("Shell density file:   ",a)') densfile(1:len_trim(densfile))
	else if(denstype.eq.'SHELL') then
	write(*,'("Powerlaw:             ",f14.3)') D%denspow
	write(9,'("Powerlaw:             ",f14.3)') D%denspow
	else if(denstype.eq.'WEDGE') then
	write(*,'("Density powerlaw:     ",f14.3)') D%denspow
	write(9,'("Density powerlaw:     ",f14.3)') D%denspow
	write(*,'("Wedge opening angle:  ",f14.3)') wedgeopen
	write(9,'("Wedge opening angle:  ",f14.3)') wedgeopen
	else if(denstype.eq.'POW') then
	write(*,'("Powerlaw:             ",f14.3)') D%denspow
	write(9,'("Powerlaw:             ",f14.3)') D%denspow
	else if(denstype.eq.'DOUBLEPOW') then
	write(*,'("Powerlaw for small R: ",f14.3)') D%denspow
	write(9,'("Powerlaw for small R: ",f14.3)') D%denspow
	write(*,'("Powerlaw for large R: ",f14.3)') D%denspow2
	write(9,'("Powerlaw for large R: ",f14.3)') D%denspow2
	write(*,'("Turnover point:       ",f14.3," AU")') D%Rpow2
	write(9,'("Turnover point:       ",f14.3," AU")') D%Rpow2
	if (D%Rpow3 .ne. 0d0) then ! Gijsexp
	   write(*,'("Turnover point 2:       ",f14.3," AU")') D%Rpow3
	   write(9,'("Turnover point 2:       ",f14.3," AU")') D%Rpow3
	endif
	else if(denstype.eq.'SIMILARITY') then
	write(*,'("Powerlaw for small R: ",f14.3)') D%denspow
	write(9,'("Powerlaw for small R: ",f14.3)') D%denspow
	write(*,'("Turnover point:       ",f14.3," AU")') D%Rexp
	write(9,'("Turnover point:       ",f14.3," AU")') D%Rexp
	write(*,'("Exponent:             ",f14.3," AU")') D%gamma_exp
	write(9,'("Exponent:             ",f14.3," AU")') D%gamma_exp
	else if(denstype.eq.'DOUBLEPOWSIM') then
	write(*,'("Powerlaw for small R: ",f14.3)') D%denspow
	write(9,'("Powerlaw for small R: ",f14.3)') D%denspow
	write(*,'("Powerlaw for large R: ",f14.3)') D%denspow2
	write(9,'("Powerlaw for large R: ",f14.3)') D%denspow2
	write(*,'("Turnover powerlaw:    ",f14.3," AU")') D%Rpow2
	write(9,'("Turnover powerlaw:    ",f14.3," AU")') D%Rpow2
	write(*,'("Turnover similarity:  ",f14.3," AU")') D%Rexp
	write(9,'("Turnover similarity:  ",f14.3," AU")') D%Rexp
	else if(denstype.eq.'MIN'.or.denstype.eq.'PINTE') then
	write(*,'("Powerlaw density:     ",f14.3)') D%denspow
	write(*,'("Powerlaw scaleheight: ",f14.3)') D%shpow
	write(*,'("Scaleheight at 1AU:   ",f14.3)') D%sh1AU
	write(9,'("Powerlaw density:     ",f14.3)') D%denspow
	write(9,'("Powerlaw scaleheight: ",f14.3)') D%shpow
	write(9,'("Scaleheight at 1AU:   ",f14.3)') D%sh1AU
	else if(denstype.eq.'PARAMETERIZED') then
	write(*,'("Powerlaw density:     ",f14.3)') D%denspow
	write(*,'("Powerlaw scaleheight: ",f14.3)') D%shpow
	write(*,'("Scaleheight at 1AU:   ",f14.3)') D%sh1AU
	write(*,'("Scaling inner rim:    1+",f6.3,"exp(-((Rin-R)/",f6.3,")^2)")') rimscale,rimwidth
	write(9,'("Powerlaw density:     ",f14.3)') D%denspow
	write(9,'("Powerlaw scaleheight: ",f14.3)') D%shpow
	write(9,'("Scaleheight at 1AU:   ",f14.3)') D%sh1AU
	write(9,'("Scaling inner rim:    1+",f6.3,"exp(-((Rin-R)/",f6.3,")^2)")') rimscale,rimwidth
	else if(denstype.eq.'PASCUCCI') then
	write(*,'("Powerlaw:             ",f14.3)') D%denspow
	write(*,'("Beta:                 ",f14.3)') D%denspow+1.125
	write(9,'("Powerlaw:             ",f14.3)') D%denspow
	write(9,'("Beta:                 ",f14.3)') D%denspow+1.125
	else if(denstype.eq.'MASSLOSS') then
	write(*,'("Mass loss:            ",e14.3," Msun/yr")') D%Mdot
	write(*,'("Escape velocity:      ",f14.3," km/s")') vexp1
	write(9,'("Mass loss:            ",e14.3," Msun/yr")') D%Mdot
	write(9,'("Escape velocity:      ",f14.3," km/s")') vexp1
	if(vexp2.gt.0d0) then
	write(*,'("    at infinity:      ",f14.3," km/s")') vexp2
	write(9,'("    at infinity:      ",f14.3," km/s")') vexp2
	write(*,'("    beta value:       ",f14.3)') bexp
	write(9,'("    beta value:       ",f14.3)') bexp
	endif
	else if(denstype.eq.'INFALL') then
	write(*,'("Mass infall rate:     ",e14.3," Msun/yr")') D%Minfall
	write(*,'("Centrifugal radius:   ",f14.3," AU")') D%Rexp
	write(9,'("Mass infall rate:     ",e14.3," Msun/yr")') D%Minfall
	write(9,'("Centrifugal radius:   ",f14.3," AU")') D%Rexp
	else if(denstype.eq.'MEIXNER') then
	write(*,'("Meixner A:            ",f14.3)') MeixA
	write(*,'("Meixner B:            ",f14.3)') MeixB
	write(*,'("Meixner C:            ",f14.3)') MeixC
	write(*,'("Meixner D:            ",f14.3)') MeixD
	write(*,'("Meixner E:            ",f14.3)') MeixE
	write(*,'("Meixner F:            ",f14.3)') MeixF
	write(9,'("Meixner A:            ",f14.3)') MeixA
	write(9,'("Meixner B:            ",f14.3)') MeixB
	write(9,'("Meixner C:            ",f14.3)') MeixC
	write(9,'("Meixner D:            ",f14.3)') MeixD
	write(9,'("Meixner E:            ",f14.3)') MeixE
	write(9,'("Meixner F:            ",f14.3)') MeixF
	write(*,'("Meixner G:            ",f14.3)') MeixG
	write(9,'("Meixner G:            ",f14.3)') MeixG
	write(*,'("Superwind radius:     ",f14.3," AU")') MeixRsw
	write(9,'("Superwind radius:     ",f14.3," AU")') MeixRsw
	if(timeshift.ne.0d0) then
	write(*,'("Time shift:           ",f14.3," years")') timeshift
	write(9,'("Time shift:           ",f14.3," years")') timeshift
	endif
	else if(denstype.eq.'BETAPIC') then
	else if(denstype.eq.'PREVIOUS') then
	else if(denstype.eq.'PRODIMO') then
	write(*,'("File (ProDiMo):       ",a)') densfile(1:len_trim(densfile))
	write(9,'("File (ProDiMo):       ",a)') densfile(1:len_trim(densfile))
	else if(denstype.eq.'ZONES') then
	do i=1,nzones
		write(*,'("Zone      ",i4,":")') i
		write(9,'("Zone      ",i4,":")') i
		write(*,'("Inner radius:         ",f14.3," AU")') Zone(i)%Rin
		write(9,'("Inner radius:         ",f14.3," AU")') Zone(i)%Rin
		write(*,'("Outer radius:         ",f14.3," AU")') Zone(i)%Rout
		write(9,'("Outer radius:         ",f14.3," AU")') Zone(i)%Rout
		write(*,'("Dust mass:            ",e18.3," Msun")') Zone(i)%Mdust
		write(9,'("Dust mass:            ",e18.3," Msun")') Zone(i)%Mdust
		write(*,'("Powerlaw:             ",f14.3)') Zone(i)%denspow
		write(9,'("Powerlaw:             ",f14.3)') Zone(i)%denspow
		write(*,'("Alpha turbulence:     ",f14.3)') Zone(i)%alphaturb
		write(9,'("Alpha turbulence:     ",f14.3)') Zone(i)%alphaturb
		if(Zone(i)%Rexp.lt.1d100) then
			write(*,'("Exponential cutoff:   ",f14.3," AU")') Zone(i)%Rexp
			write(9,'("Exponential cutoff:   ",f14.3," AU")') Zone(i)%Rexp
			write(*,'("Exponent:             ",f14.3," AU")') Zone(i)%gamma_exp
			write(9,'("Exponent:             ",f14.3," AU")') Zone(i)%gamma_exp
		endif
		if(Zone(i)%fix_struct) then
			write(*,'("Reference radius:     ",f14.3," AU")') Zone(i)%Rsh
			write(9,'("Reference radius:     ",f14.3," AU")') Zone(i)%Rsh
			write(*,'("Scaleheight:          ",f14.5," AU")') Zone(i)%sh
			write(9,'("Scaleheight:          ",f14.5," AU")') Zone(i)%sh
			write(*,'("Scaleheight powerlaw: ",f14.3," Msun")') Zone(i)%shpow
			write(9,'("Scaleheight powerlaw: ",f14.3," Msun")') Zone(i)%shpow
		else
			write(*,'("Solving vertical structure")')
			write(9,'("Solving vertical structure")')
		endif
		if(Zone(i)%sizedis) then
			write(*,'(a26,f12.4)') "Lower bound size grid:",Zone(i)%a_min
			write(9,'(a26,f12.4)') "Lower bound size grid:",Zone(i)%a_min
			write(*,'(a26,f12.4)') "Upper bound size grid:",Zone(i)%a_max
			write(9,'(a26,f12.4)') "Upper bound size grid:",Zone(i)%a_max
			write(*,'(a26,f12.2)') "Index for size powerlaw:",Zone(i)%a_pow
			write(9,'(a26,f12.2)') "Index for size powerlaw:",Zone(i)%a_pow
		else
			do j=1,ngrains
				if(Zone(i)%inc_grain(j)) then
					write(*,'("Abundance:            ",f14.3)') Zone(i)%abun(j)
					write(9,'("Abundance:            ",f14.3)') Zone(i)%abun(j)
				endif
			enddo
		endif
		if(use_qhp) then
			write(*,'(a26,f12.4)') "Fraction of QHP:",Zone(i)%fPAH
			write(9,'(a26,f12.4)') "Fraction of QHP:",Zone(i)%fPAH
		endif
		if(abs(Zone(i)%pertA).gt.1d-15) then
			write(*,'(a26,f12.4)') "Pertubation amplitude:",Zone(i)%pertA
			write(*,'(a26,f12.4)') "Pertubation radius:",Zone(i)%pertR
			write(*,'(a26,f12.4)') "Pertubation starting radius:",Zone(i)%pertR0
			write(9,'(a26,f12.4)') "Pertubation amplitude:",Zone(i)%pertA
			write(9,'(a26,f12.4)') "Pertubation radius:",Zone(i)%pertR
			write(9,'(a26,f12.4)') "Pertubation starting radius:",Zone(i)%pertR0
		endif
	enddo
	else
	write(*,'("Error in density type")')
	write(9,'("Error in density type")')
	stop
	endif

	do i=1,ngap
	   write(*,'("Gap in the disk from: ",f14.3," AU")') gap1(i)
	   write(*,'("to:                   ",f14.3," AU")') gap2(i)
	   write(*,'("Density decrease:     ",e14.3)') gap(i)
	   if (gapshape(i).ne.0) then
	      write(*,'("Shape of the gap:     ",f14.3)') gapshape(i)
	   else if (gaproundtype(i).eq.'softedge') then
	      write(*,'("Gap with soft edge ")')
	   else if (gaproundtype(i).eq.'powerlaw') then
	      write(*,'("Round off gap with p: ",f14.3)') gaproundpow(i)
	   else if (gaproundtype(i).eq.'hydro') then
	      write(*,'("Hydro gap with width: ",f14.3)') gaproundpow(i)
	   else
	      write(*,'("Vertical gap")')
	   endif
	   
	   write(9,'("Gap in the disk from: ",f14.3," AU")') gap1(i)
	   write(9,'("to:                   ",f14.3," AU")') gap2(i)
	   write(9,'("Density decrease:     ",e14.3)') gap(i)
	   if (gapshape(i).ne.0) then
	      write(9,'("Shape of the gap:     ",f14.3)') gapshape(i)
	   else if (gaproundtype(i).eq.'softedge') then
	      write(9,'("Gap with soft edge ")')
	   else if (gaproundtype(i).eq.'powerlaw') then
	      write(9,'("Round off gap with p: ",f14.3)') gaproundpow(i)
	   else if (gaproundtype(i).eq.'hydro') then
	      write(9,'("Hydro gap with width: ",f14.3)') gaproundpow(i)
	   else
	      write(9,'("Vertical gap")')
	   endif
	   
	enddo

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	write(*,'("Destroy high T dust:  ",l10)') tdes_iter
	write(9,'("Destroy high T dust:  ",l10)') tdes_iter
	if(tdes_iter) then
	do ii=1,ngrains
	write(*,'("Grain species nr.:    ",i10)') ii
	write(*,'("Maximum temperature:  ",f12.7,f12.7)') TdesA(ii),TdesB(ii)
	write(9,'("Grain species nr.:    ",i10)') ii
	write(9,'("Maximum temperature:  ",f12.7,f12.7)') TdesA(ii),TdesB(ii)
	if(force_vert_gf(ii)) then
	write(*,'("Forcing vertical gasfraction gradient")')
	write(9,'("Forcing vertical gasfraction gradient")')
	else
	write(*,'("Not forcing vertical gasfraction gradient")')
	write(9,'("Not forcing vertical gasfraction gradient")')
	endif
	if(tdes_fast(ii).gt.0d0) then
	write(*,'("Fast sublimation/condensation factor:",f14.3)') tdes_fast(ii)
	write(9,'("Fast sublimation/condensation factor:",f14.3)') tdes_fast(ii)
	endif
	enddo
	write(*,'("Convergence:          ",f14.3)') epsiter
	write(9,'("Convergence:          ",f14.3)') epsiter
	write(*,'("Maximum iterations:   ",i10)') maxiter
	write(9,'("Maximum iterations:   ",i10)') maxiter
	write(*,'("Weight of new iter:   ",f14.3)') f_weight
	write(9,'("Weight of new iter:   ",f14.3)') f_weight
	if(niter0.lt.maxiter) then
	write(*,'("Averaging after iter: ",i10)') niter0
	write(9,'("Averaging after iter: ",i10)') niter0
	endif
	endif
	
	write(*,'("Iterating structure:  ",l10)') struct_iter
	write(9,'("Iterating structure:  ",l10)') struct_iter
	if(struct_iter) then
	write(*,'("Convergence:          ",f14.3)') epsiter
	write(9,'("Convergence:          ",f14.3)') epsiter
	write(*,'("Maximum iterations:   ",i10)') maxiter
	write(9,'("Maximum iterations:   ",i10)') maxiter
	write(*,'("Weight of new iter:   ",f14.3)') f_weight
	write(9,'("Weight of new iter:   ",f14.3)') f_weight
	if(scalesh.eq.'FILE') then
	write(*,'("Scale scaleheight:    FILE")')
	write(*,'("File:                 ",a)') shscalefile(1:len_trim(shscalefile))
	write(9,'("Scale scaleheight:    FILE")')
	write(9,'("File:                 ",a)') shscalefile(1:len_trim(shscalefile))
	else if(scalesh.eq.'VALUE') then
	write(*,'("Scale scaleheight:    VALUE")')
	write(*,'("Scaling value:        ",f14.3)') shscalevalue
	write(9,'("Scale scaleheight:    VALUE")')
	write(9,'("Scaling value:        ",f14.3)') shscalevalue
	else if(scalesh.eq.'POW') then
	write(*,'("Scale scaleheight:    POW")')
	write(*,'("Scaling value:        1+",f6.3,"/(R^",f6.3,")")') shscalevalue-1,shpow
	write(9,'("Scale scaleheight:    POW")')
	write(9,'("Scaling value:        1+",f6.3,"/(R^",f6.3,")")') shscalevalue-1,shpow
	endif
	endif ! Gijsexp: can also use settling if iter=.false.

C       Gijsexp
	!
	!  Isothermal vertical structure (using T=Tmidplane)
	!
	if (mpstr) then
	   write(*,'("Using isothermal vertical structure")')
	   write(9,'("Using isothermal vertical structure")')
	endif
	!
	!  Selfconsistent dust settling
	!
	if (scset) then
	   write(*,'("--------------------------------------------------------")')
	   write(9,'("--------------------------------------------------------")')
	   write(*,'("Using selfconsistent dust settling")')
	   write(9,'("Using selfconsistent dust settling")')
	   if (scsetsave) then
	      write(*,'("Saving convergence times")')
	      write(9,'("Saving convergence times")')
	   endif
	   write(*,'("Alpha (Turbulent)         ",e14.3)') alphaturb
	   write(9,'("Alpha (Turbulent)         ",e14.3)') alphaturb
	   write(*,'("q (Turbulent)             ",f8.1)') qturb
	   write(9,'("q (Turbulent)             ",f8.1)') qturb
	   if (deadzone) then
	      write(*,'("Deadzone column:          ",f9.2)') deadcolumn 
	      write(9,'("Deadzone column:          ",f9.2)') deadcolumn 
	      write(*,'("Deadzone temperature:     ",f9.2)') deadtemp
	      write(9,'("Deadzone temperature:     ",f9.2)') deadtemp
	      write(*,'("Deadzone alpha:           ",e14.3)') deadalpha
	      write(9,'("Deadzone alpha:           ",e14.3)') deadalpha
	   endif
	   if (thinparticle.ne.-2) then
	      write(*,'("Using particle ",i2," in upper layers")') thinparticle
	      write(9,'("Using particle ",i2," in upper layers")') thinparticle
	   endif
	endif
	!
	!  Dust settling using scale height scaling (Dubrulle et al 95)
	!
	if (mpset) then
	   write(*,'("Using isothermal dust settling")')
	   write(9,'("Using isothermal dust settling")')
		if(denstype.eq.'ZONES') then
	   write(*,'("Turbulent alpha set in zones")')
	   write(9,'("Turbulent alpha set in zones")')
		else
	   write(*,'("Alpha (Turbulent)         ",e14.3)') alphaturb
	   write(9,'("Alpha (Turbulent)         ",e14.3)') alphaturb
		endif
	   write(*,'("q (Turbulent)             ",f8.1)') qturb
	   write(9,'("q (Turbulent)             ",f8.1)') qturb
	endif

	if (fixmpset) then
	   write(*,'("Using isothermal dust settling with prescribed soundspeed")')
	   write(9,'("Using isothermal dust settling with prescribed soundspeed")')
		if(denstype.eq.'ZONES') then
	   write(*,'("Turbulent alpha set in zones")')
	   write(9,'("Turbulent alpha set in zones")')
		else
	   write(*,'("Alpha (Turbulent)         ",e14.3)') alphaturb
	   write(9,'("Alpha (Turbulent)         ",e14.3)') alphaturb
		endif
	   write(*,'("q (Turbulent)             ",f8.1)') qturb
	   write(9,'("q (Turbulent)             ",f8.1)') qturb
	endif
	!
	!   Grain size distribution from MRN
	!
	if (mrn) then
	   write(*,'("--------------------------------------------------------")')
	   write(9,'("--------------------------------------------------------")')
	   write(*,'("Calculating grain size distribution from MRN")')
	   write(9,'("Calculating grain size distribution from MRN")')

	   write(*,'(a26,f12.4)') "Lower bound size grid:",mrn_rmin*1d4
	   write(9,'(a26,f12.4)') "Lower bound size grid:",mrn_rmin*1d4
	   write(*,'(a26,f12.4)') "Upper bound size grid:",mrn_rmax*1d4
	   write(9,'(a26,f12.4)') "Upper bound size grid:",mrn_rmax*1d4
	   write(*,'(a26,f12.2)') "Index for size powerlaw:",mrn_index
	   write(9,'(a26,f12.2)') "Index for size powerlaw:",mrn_index
	   if (mrn_ngrains .gt. 0) then
	   	  if (mrn_ngrains+1 .eq. ngrains) then
   	         write(*,'(a16,i2)') "Excluding grain ", ngrains
	         write(9,'(a16,i2)') "Excluding grain ", ngrains
	      else	
   	         write(*,'(a17,i2,a4,i2)') "Excluding grains ",mrn_ngrains+1," to ", ngrains
	         write(9,'(a17,i2,a4,i2)') "Excluding grains ",mrn_ngrains+1," to ", ngrains
	      endif	
	   endif	
	endif
	!
	!   Grain size distribution from coagulation/fragmentation equilibrium
	!
	if (gsd) then
	   write(*,'("--------------------------------------------------------")')
	   write(9,'("--------------------------------------------------------")')
	   write(*,'("Calculating grain size distribution (Experimental)")')
	   write(9,'("Calculating grain size distribution (Experimental)")')
	   write(*,'("Please credit Til Birnstiel if using this part of the code")')
	   write(9,'("Please credit Til Birnstiel if using this part of the code")')
	   if (gsd_full) then
	      write(*,'("Calculating distribution -> not yet")')
	      write(9,'("Calculating distribution -> not yet")')
	      stop 65456
	   else
	      write(*,'("Using the analytical approximation")')
	      write(9,'("Using the analytical approximation")')	      
	   endif
	      write(*,'(a26,f12.4)') "Lower bound size grid:",gsd_rmin*1d4
	      write(9,'(a26,f12.4)') "Lower bound size grid:",gsd_rmin*1d4
	      write(*,'(a26,f12.4)') "Upper bound size grid:",gsd_rmax*1d4
	      write(9,'(a26,f12.4)') "Upper bound size grid:",gsd_rmax*1d4
	      write(*,'(a26,f12.4)') "Fragmentation velocity:",gsd_vfrag
	      write(9,'(a26,f12.4)') "Fragmentation velocity:",gsd_vfrag
	      write(*,'(a26,f12.4)') "Fragmentation slope:",gsd_xi
	      write(9,'(a26,f12.4)') "Fragmentation slope:",gsd_xi
	endif
C       end 

	if(radpress) then
	   write(*,'("Including radiation pressure")')
	   write(9,'("Including radiation pressure")')
	endif	
c	endif ! Gijsexp: can also use settling if iter=.false.
	if(rimscale.ne.0d0) then
	write(*,'("Scaling inner rim:    1+",f6.3,"exp(-((Rin-R)/",f6.3,")^2)")') rimscale,rimwidth
	write(9,'("Scaling inner rim:    1+",f6.3,"exp(-((Rin-R)/",f6.3,")^2)")') rimscale,rimwidth
	endif

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	if(lamtype(1:4).eq.'FILE') then
	write(*,'("Wavelength file:      ",a)') gridfile(1:len_trim(gridfile))
	write(9,'("Wavelength file:      ",a)') gridfile(1:len_trim(gridfile))
	else
	write(*,'("Minimum wavelength:   ",f14.3," micron")') lam1
	write(*,'("Maximum wavelength:   ",f14.3," micron")') lam2
	write(*,'("Number of points:     ",i10)') nlam
	write(9,'("Minimum wavelength:   ",f14.3," micron")') lam1
	write(9,'("Maximum wavelength:   ",f14.3," micron")') lam2
	write(9,'("Number of points:     ",i10)') nlam
	if(nzlam.gt.0) then
	write(*,'("Minimum zoom-in wavel:",f14.3," micron")') zlam1
	write(*,'("Maximum zoom-in wavel:",f14.3," micron")') zlam2
	write(*,'("Number of points:     ",i10)') nzlam
	write(9,'("Minimum zoom-in wavel:",f14.3," micron")') zlam1
	write(9,'("Maximum zoom-in wavel:",f14.3," micron")') zlam2
	write(9,'("Number of points:     ",i10)') nzlam
	endif
	endif

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	if(overflow) then
	write(*,'("Maximum interactions: ",i10)') maxinteract
	write(9,'("Maximum interactions: ",i10)') maxinteract
	endif
	write(*,'("Limit for diffusion:  ",i10," photon packages")') NphotDiffuse
	write(9,'("Limit for diffusion:  ",i10," photon packages")') NphotDiffuse
	write(*,'("Limit for diffusion:  ",f5.3," (T error)")') dTDiffuse
	write(9,'("Limit for diffusion:  ",f5.3," (T error)")') dTDiffuse
	if(FLD) then
	write(*,'("Using Flux Limited Diffusion")')
	write(9,'("Using Flux Limited Diffusion")')
	endif
	if(RNDW) then
	write(*,'("Using Random Walk procedure")')
	write(*,'("Random Walk factor:   ",f14.3)') factRW
	write(9,'("Using Random Walk procedure")')
	write(9,'("Random Walk factor:   ",f14.3)') factRW
	endif
	if(multiwav) then
	write(*,'("Using multi wavelength approach for MC spectra (beta!)")')
	write(9,'("Using multi wavelength approach for MC spectra (beta!)")')
	endif
	if(convection) then
	write(*,'("Correcting the temperature structure for convection")')
	write(9,'("Correcting the temperature structure for convection")')
	endif
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	if(opacity_set) then
c	if(opacity_file_set) then
c	write(*,'("Opacity file:         ",a)') opacityfile(1:len_trim(opacityfile))
c	write(9,'("Opacity file:         ",a)') opacityfile(1:len_trim(opacityfile))
c	else

	!  Set abundances to 1/ngrains if not supplied	
	if(.not.arg_abun .and. .not.(mrn.or.gsd) .and. nzones.eq.0) then
	   write(*,'("Abundances not set, using equal weights")')
	   write(9,'("Abundances not set, using equal weights")')
	   do i=1,ngrains
	      warg(i)=1d0/real(ngrains)
	   enddo
	   arg_abun=.true.
	endif

	!  Set abundances according to mrn distribution
	!  This is also the initial guess for gsd=.true.
	if(mrn.or.gsd) then
		call gsd_MRN(rgrain(1:ngrains),warg(1:ngrains))	! GijsExp
		warg=warg/sum(warg)
	endif

	do i=1,ngrains  !  Start displaying particle information

	!  Display particle file locations
	if(compositionfile.eq.' ') then
	if(parttype(i).eq.1) then
	write(*,'("Particle file:        ",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	write(9,'("Particle file:        ",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	else if(parttype(i).eq.2) then
	write(*,'("File with opacities:  ",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	write(9,'("File with opacities:  ",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	else if(parttype(i).eq.3) then
	write(*,'("Mixed aggregates:     ",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	write(*,'("Mixed aggr. powelaw:  ",f6.2)') powmix(i)
	write(9,'("Mixed aggregates:     ",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	write(9,'("Mixed aggr. powelaw:  ",f6.2)') powmix(i)
	if(Tmix(i).lt.real(TMAX)*dT) then
	write(*,'("Mixed aggr. temp.:    ",f6.1," K")') Tmix(i)
	write(9,'("Mixed aggr. temp.:    ",f6.1," K")') Tmix(i)
	else
	write(*,'("Mixed aggr. radius:   ",f6.2," AU")') radmix(i)
	write(9,'("Mixed aggr. radius:   ",f6.2," AU")') radmix(i)
	endif
	else if(parttype(i).eq.4) then
	write(*,'("Quantum heated part.: ",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	write(9,'("Quantum heated part.: ",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	else if(parttype(i).eq.5) then
	write(*,'("T dependent opacities:",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	write(9,'("T dependent opacities:",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	else if(parttype(i).eq.6) then
	write(*,'("FITS file particle   :",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	write(9,'("FITS file particle   :",a," (",f6.2,"%)")') partarg(i)(1:len_trim(partarg(i))),100d0*warg(i)/sum(warg(1:ngrains))
	endif
	else
	if(parttype(i).eq.1) then
	write(*,'("Particle file:        ",a)') partarg(i)(1:len_trim(partarg(i)))
	write(9,'("Particle file:        ",a)') partarg(i)(1:len_trim(partarg(i)))
	else if(parttype(i).eq.2) then
	write(*,'("File with opacities:  ",a)') partarg(i)(1:len_trim(partarg(i)))
	write(9,'("File with opacities:  ",a)') partarg(i)(1:len_trim(partarg(i)))
	else if(parttype(i).eq.3) then
	write(*,'("Mixed aggregates:     ",a)') partarg(i)(1:len_trim(partarg(i)))
	write(*,'("Mixed aggr. powelaw:  ",f6.2)') powmix(i)
	write(9,'("Mixed aggregates:     ",a)') partarg(i)(1:len_trim(partarg(i)))
	write(9,'("Mixed aggr. powelaw:  ",f6.2)') powmix(i)
	if(Tmix(i).lt.real(TMAX)*dT) then
	write(*,'("Mixed aggr. temp.:    ",f6.1," K")') Tmix(i)
	write(9,'("Mixed aggr. temp.:    ",f6.1," K")') Tmix(i)
	else
	write(*,'("Mixed aggr. radius:   ",f6.2," AU")') radmix(i)
	write(9,'("Mixed aggr. radius:   ",f6.2," AU")') radmix(i)
	endif
	else if(parttype(i).eq.4) then
	write(*,'("Quantum heated part.: ",a)') partarg(i)(1:len_trim(partarg(i)))
	write(9,'("Quantum heated part.: ",a)') partarg(i)(1:len_trim(partarg(i)))
	else if(parttype(i).eq.5) then
	write(*,'("T dependent opacities:",a)') partarg(i)(1:len_trim(partarg(i)))
	write(9,'("T dependent opacities:",a)') partarg(i)(1:len_trim(partarg(i)))
	else if(parttype(i).eq.6) then
	write(*,'("FITS file particle   :",a)') partarg(i)(1:len_trim(partarg(i)))
	write(9,'("FITS file particle   :",a)') partarg(i)(1:len_trim(partarg(i)))
	endif
	endif

	! Display HG scattering properties (if set manually)
	if(asym(i).lt.1d0.and.asym(i).gt.-1d0) then
	   write(*,'("Asymmetry parameter g:",f4.2)') asym(i)
	   write(9,'("Asymmetry parameter g:",f4.2)') asym(i)	   
	   if(wasym2(i).ne.0d0) then
	      write(*,'("Second g:             ",f4.2)') asym2(i)
	      write(9,'("Second g:             ",f4.2)') asym2(i)
	      write(*,'("Weight of second g: ",f6.2)') wasym2(i)	      
	      write(9,'("Weight of second g: ",f6.2)') wasym2(i)	      
	   endif
	endif

	!  Display vertical structure per particle
	if(settle(i)) then
	if(shtype(i).eq.'DISK') then
	if(settlefile(i).eq.' ') then
	write(*,'("Settling to:          ",f5.3," scaleheights")') part_shscale(i)
	write(9,'("Settling to:          ",f5.3," scaleheights")') part_shscale(i)
	else
	write(*,'("Settling file:        ",a)') settlefile(i)(1:len_trim(settlefile(i)))
	write(9,'("Settling file:        ",a)') settlefile(i)(1:len_trim(settlefile(i)))
	endif
	else if(shtype(i).eq.'HALO') then
	write(*,'("Spherical halo")')
	write(9,'("Spherical halo")')
	else if(shtype(i).eq.'LOBE') then
	write(*,'("Lobe wigth opening:   ",f7.3)') part_shscale(i)
	write(9,'("Lobe wigth opening:   ",f7.3)') part_shscale(i)
	else if(shtype(i).eq.'JET') then
	write(*,'("Jet wigth opening:    ",f7.3)') part_shscale(i)
	write(9,'("Jet wigth opening:    ",f7.3)') part_shscale(i)
	else if(shtype(i).eq.'ZODI') then
	write(*,'("Zodiacal dust distribution")')
	write(9,'("Zodiacal dust distribution")')
	else if(shtype(i).eq.'CYLI') then
	write(*,'("Cylinder of height:",f10.3," AU")') part_shscale(i)
	write(9,'("Cylinder of height:",f10.3," AU")') part_shscale(i)
	else if(shtype(i).eq.'CAVITY') then
	write(*,'("Cylindrical cavity radius:",f10.3," AU")') part_shscale(i)
	write(9,'("Cylindrical cavity radius:",f10.3," AU")') part_shscale(i)
	endif
	endif

	!  Display particle properties
	write(*,'("Particle material:    ",a)') material(i)(1:len_trim(material(i)))
	write(9,'("Particle material:    ",a)') material(i)(1:len_trim(material(i)))
	if(scset.or.mpset) then
	   write(*,'("Particle radius(micron):",f12.4)') rgrain(i) * 1d4
	   write(9,'("Particle radius(micron):",f12.4)') rgrain(i) * 1d4
	   write(*,'("Particle density:               ",f4.2)') rhograin(i)
	   write(9,'("Particle density:               ",f4.2)') rhograin(i)
	endif

	!  Display particle inner and outer radii and gap shapes
	if(minrad(i).gt.D%Rin) then
	write(*,'("Destroying inside:    ",f6.2," AU")') minrad(i)
	write(9,'("Destroying inside:    ",f6.2," AU")') minrad(i)
	endif
	if(maxrad(i).lt.D%Rout) then
	write(*,'("Destroying outside:   ",f6.2," AU")') maxrad(i)
	write(9,'("Destroying outside:   ",f6.2," AU")') maxrad(i)
	endif
	if(roundtype(i).ne.' ') then ! Gijsexp
	   write(*,'("  Rounding of the inner wall at ",f6.2," AU")') max(minrad(i),D%Rin)
	   write(9,'("  Rounding of the inner wall at ",f6.2," AU")') max(minrad(i),D%Rin)
	endif
	if(roundtype(i).eq.'powerlaw') then ! Gijsexp
	   write(*,'("  width: ",f6.2," AU, power law index: ",f6.2)') roundwidth(i),roundpow(i)
	   write(9,'("  width: ",f6.2," AU, power law index: ",f6.2)') roundwidth(i),roundpow(i)
	else if(roundtype(i).eq.'hydro') then ! Gijsexp
	   write(*,'("  width: ",f6.2," AU, power law index: ",f6.2)') roundwidth(i),roundpow(i)
	   write(9,'("  width: ",f6.2," AU, power law index: ",f6.2)') roundwidth(i),roundpow(i)
	else if(roundtype(i).eq.'softedge') then
	   write(*,'("Using a soft edge")') 
	   write(9,'("Using a soft edge")')
	endif
	enddo ! loop over ngrains
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	if(compositionfile.ne.' ') then
	write(*,'("Composition file:     ",a)') compositionfile(1:len_trim(compositionfile))
	write(9,'("Composition file:     ",a)') compositionfile(1:len_trim(compositionfile))
	endif
c	endif
	else
	write(*,'("Opacity type:         ",a10)') 'POWERLAW  '
	write(9,'("Opacity type:         ",a10)') 'POWERLAW  '
	endif
	write(*,'("Scattering:           ",a10)') scattype(1:len_trim(scattype))
	write(9,'("Scattering:           ",a10)') scattype(1:len_trim(scattype))
	if(storescatt) then
	if(scat_how.ne.1.and.scat_how.ne.0) then
	write(*,'("Incompatible scatt type for storage of scatt-field")')
	write(9,'("Incompatible scatt type for storage of scatt-field")')
	storescatt=.false.
	else if(scat_how.ne.0) then
	write(*,'("Storing local scattered field")')
	write(9,'("Storing local scattered field")')
	else if(scat_how.eq.0) then
	storescatt=.false.
	endif
	endif
	write(*,'("Thermal contact:      ",l10)') tcontact
	write(9,'("Thermal contact:      ",l10)') tcontact
	if(npow.ne.0) then
	do i=1,npow
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Abundance gradient:")')
	write(*,'("Closeby particle:     ",i10)') powclose(i)
	write(*,'("Far away particle:    ",i10)') powfar(i)
	write(*,'("Powerlaw starting R:  ",f14.3," AU")') powrad0(i)
	write(*,'("Powerlaw slope:       ",f14.3)') powslope(i)
	if (powinner(i).lt.1d0) write(*,'("Abundance decrease factor at small R: ",f14.3)') powinner(i)

	write(9,'("Abundance gradient:")')
	write(9,'("Closeby particle:     ",i10)') powclose(i)
	write(9,'("Far away particle:    ",i10)') powfar(i)
	write(9,'("Powerlaw starting R:  ",f14.3," AU")') powrad0(i)
	write(9,'("Powerlaw slope:       ",f14.3)') powslope(i)
	if (powinner(i).lt.1d0) write(9,'("Abundance decrease factor at small R: ",f14.3)') powinner(i)
	if (warg(powfar(i)) .eq. 0d0) then 
	   write(*,'("Warning: Abundance of Faraway particle should not be zero")')
	   write(9,'("Warning: Abundance of Faraway particle should not be zero")')
	endif
	enddo
	endif

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	if(computeTgas.or.viscous.or.denstype.eq.'PRODIMO') useTgas=.true.

	D%Mstar=D%Mstar*Msun
	D%Rstar=D%Rstar*Rsun
	D%Rstar2=D%Rstar2*Rsun
	D%distance=D%distance*parsec
	D%Mtot=D%Mtot*Msun
	
	D%Mdot=D%Mdot*Msun/(365.25d0*24d0*60d0*60d0)
	D%Minfall=D%Minfall*Msun/(365.25d0*24d0*60d0*60d0)
	if(vexp2.le.0d0) vexp2=vexp1
	vexp1=vexp1*100000d0
	vexp2=vexp2*100000d0

c first setup the wavelength grid
	if(lamtype(1:4).eq.'FILE') then

	i=1
	inquire(file=gridfile,exist=truefalse)
	if(.not.truefalse) then
		write(*,'("Wavelength grid not found")')
		write(9,'("Wavelength grid not found")')
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		stop
	endif
	open(unit=20,file=gridfile,RECL=100)
1	read(20,*,end=2)
	i=i+1
	goto 1
2	close(unit=20)
	nlam=i-1
	if(.not.arraysallocated) then
		allocate(lam(nlam))
		allocate(BB(nlam,0:TMAX))
		allocate(dnu(nlam))
		allocate(nu(nlam))
		allocate(D%Fstar(nlam))
	endif
	open(unit=20,file=gridfile,RECL=100)
	do i=1,nlam
		read(20,*) lam(i)
	enddo
	close(unit=20)

	else
	
	if(.not.arraysallocated) then
		allocate(lam(nlam))
		allocate(BB(nlam,0:TMAX))
		allocate(dnu(nlam))
		allocate(nu(nlam))
		allocate(D%Fstar(nlam))
	endif

	if(exportProDiMo) then

	allocate(spec(0:nlam))

c	NLAM  = 70
	N1UV  = 3
	NUV   = 9
c	lmax  = 4000.0

	spec(0)    = 0.0912                !  912 Angstroem
	spec(N1UV) = 0.1110                ! 1010 Angstroem
	spec(NUV)  = 0.2050                ! 2050 Angstroem
	lam1=spec(0)

	do i=1,N1UV
		spec(i) = exp(log(spec(0))+REAL(i)/REAL(N1UV)*log(spec(N1UV)/spec(0)))
	enddo
	do i=N1UV+1,NUV
		spec(i) = exp(log(spec(N1UV))+REAL(i-N1UV)/REAL(NUV-N1UV)*log(spec(NUV)/spec(N1UV)))
	enddo
	do i=NUV+1,NLAM-nzlam                    ! interval boundaries
		spec(i) = EXP(LOG(spec(NUV))+REAL(i-NUV)/REAL(NLAM-nzlam-NUV)*LOG(lam2/spec(NUV)))
	enddo
	do i=nlam-nzlam+1,NLAM                    ! interval boundaries
		spec(i) = EXP(LOG(zlam1)+REAL(i-(nlam-nzlam))/REAL(nzlam)*LOG(zlam2/zlam1))
	enddo
	call sort(spec(0:nlam),nlam+1)
	do i=1,NLAM
		lam(i) = sqrt(spec(i)*spec(i-1))
	enddo

	deallocate(spec)

	else if(use_qhp.and.Nphot.le.0) then

	nqhp=0
	do ii=1,ngrains
		if(parttype(ii).eq.4.or.parttype(ii).eq.8) then
			nqhp=nqhp+1
			if(outputfits) then
				write(file,'(a,"QHPemis",i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),nqhp/10,nqhp-10*(nqhp/10)
			else
				write(file,'(a,"QHPemis",i1,i1,".dat")') outdir(1:len_trim(outdir)),nqhp/10,nqhp-10*(nqhp/10)
			endif
			write(*,'("Reading wavelength grid from:       ",a)') file(1:len_trim(file))
			write(9,'("Reading wavelength grid from:       ",a)') file(1:len_trim(file))
			call readstruct(file,(/'LAM    '/),1,ii,.true.)
			exit
		endif
	enddo

	else

	if(nzlam.le.0) then
		do i=1,nlam
			lam(i)=10d0**(log10(lam1)+(log10(lam2)-log10(lam1))*real(i-1)/real(nlam-1))
		enddo
	else
		nl=(nlam-nzlam)*(log10(zlam1/lam1)/log10(lam2*zlam1/(zlam2*lam1)))
		do i=1,nl
			lam(i)=10d0**(log10(lam1)+(log10(zlam1)-log10(lam1))*real(i-1)/real(nl-1))
		enddo
		j=i-1
		nl=(nlam-nzlam)-nl
		do i=1,nl
			lam(i+j)=10d0**(log10(zlam2)+(log10(lam2)-log10(zlam2))*real(i-1)/real(nl-1))
		enddo
		j=j+i-1
		do i=1,nzlam
			lam(i+j)=10d0**(log10(zlam1)+(log10(zlam2)-log10(zlam1))*real(i)/real(nzlam+1))
		enddo
		call sort(lam,nlam)
	endif

	endif	

	endif

	clight=2.9979d14
	do i=1,nlam
		nu(i)=clight/lam(i)
	enddo
	do i=2,nlam-1
		dnu(i)=(nu(i-1)-nu(i+1))/2d0
	enddo
	dnu(nlam)=(nu(nlam-1)-nu(nlam))
	dnu(1)=(nu(1)-nu(2))

	do i=1,nlam
		BB(i,0)=0d0
		do j=1,TMAX
			T=real(j)*dT
			BB(i,j)=Planck(T,lam(i))
		enddo
	enddo
	BBint(0)=0d0
	do i=1,TMAX
		call integrate(BB(1:nlam,i),BBint(i))
	enddo

	if(startype.eq.'PLANCK') then
		do i=1,nlam
			D%Fstar(i)=pi*D%Rstar**2*Planck(D%Tstar,lam(i))
		enddo
		call integrate(D%Fstar,D%Lstar)
		write(*,'("Error on the luminosity: ",f7.3," %")') 100d0*(D%Lstar-Luminosity(D%Tstar,D%Rstar))/D%Lstar
		write(9,'("Error on the luminosity: ",f7.3," %")') 100d0*(D%Lstar-Luminosity(D%Tstar,D%Rstar))/D%Lstar
		nlamHR=nlam+9
		allocate(lamHR(nlamHR))
		allocate(FstarHR(nlamHR))
		lamHR(10:nlam+9)=lam(1:nlam)
		lamHR(1)=0.01
		lamHR(2)=0.02
		lamHR(3)=0.03
		lamHR(4)=0.04
		lamHR(5)=0.05
		lamHR(6)=0.06
		lamHR(7)=0.07
		lamHR(8)=0.08
		lamHR(9)=0.09
		call sort(lamHR,nlamHR)
		do i=1,nlamHR
			FstarHR(i)=pi*D%Rstar**2*Planck(D%Tstar,lamHR(i))
		enddo
	else if(startype.eq.'KURUCZ') then
		call ReadKurucz(D%Tstar,D%logg,lam,D%Fstar,nlam)
		call integrate(D%Fstar,tot)
		D%Fstar=D%Lstar*D%Fstar/tot
		D%Rstar=sqrt(D%Lstar/Luminosity(D%Tstar,1d0))
		nlamHR=nlam+9
		allocate(lamHR(nlamHR))
		allocate(FstarHR(nlamHR))
		lamHR(10:nlam+9)=lam(1:nlam)
		lamHR(1)=0.01
		lamHR(2)=0.02
		lamHR(3)=0.03
		lamHR(4)=0.04
		lamHR(5)=0.05
		lamHR(6)=0.06
		lamHR(7)=0.07
		lamHR(8)=0.08
		lamHR(9)=0.09
		call sort(lamHR,nlamHR)
		call ReadKurucz(D%Tstar,D%logg,lamHR,FstarHR,nlamHR)
		FstarHR=D%Lstar*FstarHR/tot
	else if (startype.eq.'FILE') then
		call readstar(starfile,lam,D%Fstar,nlam)
		call readspecHR(starfile)
		call integrate(D%Fstar,tot)
		D%Fstar=D%Lstar*D%Fstar/tot
		FstarHR=D%Lstar*FstarHR/tot
	endif

	if(D%Rstar2.gt.0d0.and.D%Tstar2.gt.0d0) then
		do i=1,nlam
			D%Fstar(i)=D%Fstar(i)+pi*D%Rstar2**2*Planck(D%Tstar2,lam(i))
		enddo
		call integrate(D%Fstar,D%Lstar)
		write(*,'("Total luminosity of the secondary: ",f7.3)') Luminosity(D%Tstar2,D%Rstar2)/Luminosity(5777d0,Rsun)
		write(9,'("Total luminosity of the secondary: ",f7.3)') Luminosity(D%Tstar2,D%Rstar2)/Luminosity(5777d0,Rsun)

		write(*,'("Total luminosity of the binary: ",f7.3)') D%Lstar/Luminosity(5777d0,Rsun)
		write(9,'("Total luminosity of the binary: ",f7.3)') D%Lstar/Luminosity(5777d0,Rsun)
	endif

	! add edges of gaps, minrad and maxrad to Rfix
	do i=1,ngap
	   if(gap1(i).gt.D%Rin.and.gap1(i).lt.D%Rout.and.minval(abs(Rfix(1:nRfix)-gap1(i))).ne.0d0) then
	      nRfix=nRfix+1
	      Rfix(nRfix)=gap1(i)
	   endif
	   if(gap2(i).gt.D%Rin.and.gap2(i).lt.D%Rout.and.minval(abs(Rfix(1:nRfix)-gap2(i))).ne.0d0) then
	      nRfix=nRfix+1
	      Rfix(nRfix)=gap2(i)
	   endif
	enddo
	do ii=1,ngrains
	   if(minrad(ii).gt.D%Rin.and.minrad(ii).lt.D%Rout.and.minval(abs(Rfix(1:nRfix)-minrad(ii))).ne.0d0) then
	      nRfix=nRfix+1
	      Rfix(nRfix)=minrad(ii)
	   endif
	   if(maxrad(ii).gt.D%Rin.and.maxrad(ii).lt.D%Rout.and.minval(abs(Rfix(1:nRfix)-maxrad(ii))).ne.0d0) then
	      nRfix=nRfix+1
	      Rfix(nRfix)=maxrad(ii)
	   endif
	enddo
	
	!  Write Rfix into D%Rfix and sort it 
	D%nRfix=nRfix
	if(.not.arraysallocated) allocate(D%Rfix(D%nRfix))
	D%Rfix(1:D%nRfix)=Rfix(1:nRfix)
	if(D%nRfix.ne.0) call sort(D%Rfix,D%nRfix)
c	write(*,*) "Rfix, ",D%Rfix(1:D%nRfix)
c	write(9,*) "Rfix, ",D%Rfix(1:D%nRfix)

	if(denstype.eq.'FILE'.or.denstype.eq.'PREVIOUS'.or.startiter.ne.' ') then
		if(outputfits) then
			if(startiter.ne.' ') write(densfile,'(a,"denstemp",a,".fits.gz")') 
     &				outdir(1:len_trim(outdir)),startiter(1:len_trim(startiter))
			if(denstype.eq.'PREVIOUS') write(densfile,'(a,"denstemp.fits.gz")') 
     &				outdir(1:len_trim(outdir))
		else
			if(startiter.ne.' ') write(densfile,'(a,"denstemp",a,".dat")') 
     &				outdir(1:len_trim(outdir)),startiter(1:len_trim(startiter))
			if(denstype.eq.'PREVIOUS') write(densfile,'(a,"denstemp.dat")') 
     &				outdir(1:len_trim(outdir))
		endif
		call readstruct(densfile,(/'DENS   '/),1,0,.true.)

		D%nR=D%nR-1
		D%R(1)=D%Rin
		do i=2,D%nR
			D%R(i)=10d0**((2d0*log10(D%R_av(i-1)/AU)-log10(D%R(i-1))))
			if(D%R_av(i)/AU.lt.D%R(i).or.D%R_av(i-1)/AU.gt.D%R(i)) then
				D%R(i)=(D%R_av(i-1)/AU+D%R_av(i)/AU)/2d0
			endif
		enddo
		D%R(D%nR+1)=D%Rout
		D%R(0)=D%Rstar/AU
		call sort(D%R(1:D%nR),D%nR)
		D%nR=D%nR+1

		shscale(0:D%nR)=1d0
		if(scalesh.eq.'FILE') then
			call regrid(shscalefile,D%R_av/AU,shscale,D%nR)
		else if(scalesh.eq.'VALUE') then
			shscale(0:D%nR)=shscalevalue
		else if(scalesh.eq.'POW') then
			shscale(0:D%nR-1)=1d0+(shscalevalue-1d0)*(D%R_av(0:D%nR-1)/AU)**(-shpow)
		else
			shscale(0:D%nR)=1d0
		endif
	else

	if(radfile.ne.' ') then
		open(unit=40,file=radfile,RECL=100)
		D%nR=-1
50		read(40,*,end=60) dummy
		D%nR=D%nR+1
		goto 50
60		close(unit=40)
		D%nR=D%nR+nRfix
	endif
		
	if(.not.arraysallocated) then
		allocate(C(0:D%nR+1,0:D%nTheta+1))
		allocate(D%theta_av(0:D%nTheta+1))
		allocate(D%Theta(0:D%nTheta+1))
		allocate(D%thet(0:D%nTheta+1))
		allocate(D%SinTheta(0:D%nTheta+1))
		allocate(D%R_av(0:D%nR+1))
		allocate(D%R(0:D%nR+1))
		allocate(shscale(0:D%nR+1))
	endif

	D%R(0)=D%Rstar/AU
	D%nR=D%nR+1
	D%nTheta=D%nTheta+1

	i=0
	k=(D%nTheta-ntspan*ntlev)/2
	do j=1,k
		i=i+1
		D%Theta(i)=real(j-1)/real(k-1)
		i=i+1
		D%Theta(i)=cos(0.5d0*pi*real(j)/real(k+1))
	enddo
	call sort(D%Theta(1:i),i)
	do k=1,ntlev
	do j=1,ntspan
		i=i+1
		D%Theta(i)=(D%Theta(j)+D%Theta(j+1))/2d0
	enddo
	call sort(D%Theta(1:i),i)
	enddo
	if(i.lt.D%nTheta) then
		do j=i+1,D%nTheta
			D%Theta(j)=real(j-i+1)/real(D%nTheta-i+1)
		enddo
	endif
	call sort(D%Theta(1:D%nTheta),D%nTheta)
	do i=1,D%nTheta-1
		D%theta_av(i)=acos((D%Theta(D%nTheta-i+1)+D%Theta(D%nTheta-i))/2d0)
	enddo

	if(radfile.eq.' ') then
		nused=0
		if(nzones.gt.0) then
			fpert=8d0
92			npert=0
			do i=1,nzones
				if(abs(Zone(i)%pertA).gt.1d-50) then
					npert=npert+(Zone(i)%Rout-Zone(i)%Rin)*fpert/Zone(i)%pertR
				endif
			enddo
			if(npert.gt.(D%nR-D%nRfix)/2) then
				fpert=fpert*real(D%nR-D%nRfix)/2d0/real(npert)
				goto 92
			endif
			do i=1,nzones
				if(abs(Zone(i)%pertA).gt.1d-50) then
					npert=(Zone(i)%Rout-Zone(i)%Rin)*fpert/Zone(i)%pertR
					do j=1,npert
						if(nused.lt.D%nR) then
							nused=nused+1
							D%R(nused)=Zone(i)%Rin+(Zone(i)%Rout-Zone(i)%Rin)*(real(j)-0.5d0)/real(npert)
						endif
					enddo
				endif
			enddo
		endif
		do i=1,D%nRfix
			if(nused.lt.D%nR) then
				nused=nused+1
				D%R(nused)=D%Rfix(i)
			endif
		enddo
		do i=1,D%nR-nused
			D%R(i+nused)=10d0**(log10(D%Rin)+(log10(D%Rout)-log10(D%Rin))*real(i-1)/real(D%nR-nused-1))
		enddo
		call sort(D%R(1:D%nR),D%nR)
	else
		open(unit=40,file=radfile,RECL=100)
		do i=1,D%nR
			read(40,*,end=40) D%R(i)
		enddo
		close(unit=40)
	endif

	
	endif

	if(scale_R.ne.0) then
		RTmax=D%R(0)*(D%Tstar/Tmax_R)**2
		stretch=(D%R(D%nR)-RTmax)/(D%R(D%nR)-D%R(1))
		do i=D%nR,1,-1
			if(scale_R.eq.1) D%R(i)=D%R(i)*RTmax/D%R(1)
			if(scale_R.eq.2) D%R(i)=RTmax+stretch*(D%R(i)-D%R(1))
		enddo
	endif
	do i=0,D%nR-1
		D%R_av(i)=AU*10d0**((log10(D%R(i))+log10(D%R(i+1)))/2d0)
	enddo
c in the theta grid we actually store cos(theta) for convenience
	D%Theta(1)=1d0
	D%thet(1)=0d0
	D%SinTheta(1)=sin(D%thet(1))
	do i=2,D%nTheta-1
		D%thet(i)=(D%theta_av(i-1)+D%theta_av(i))/2d0
		D%Theta(i)=cos((D%theta_av(i-1)+D%theta_av(i))/2d0)
		D%SinTheta(i)=sin(D%thet(i))
	enddo
	D%Theta(D%nTheta)=0d0
	D%thet(D%nTheta)=pi/2d0
	D%SinTheta(D%nTheta)=sin(D%thet(D%nTheta))

	do i=1,D%nTheta-1
		D%theta_av(i)=acos((D%Theta(i)+D%Theta(i+1))/2d0)
	enddo
	

	write(thetagridfile,'(a,"thetagrid.dat")') outdir(1:len_trim(outdir))
	open(unit=60,file=thetagridfile,RECL=100)
	do i=1,D%nTheta
		write(60,*) acos(D%Theta(i))
	enddo
	close(unit=60)
	
	do j=0,D%nTheta-1
	do i=1,D%nR-1
		C(i,j)%T=0d0
		C(i,j)%Ni=0
		C(i,j)%EJv=0d0
		C(i,j)%iTD=0
		C(i,j)%alphaturb=alphaturb
		C(i,j)%xedge(1)=D%R(i)
		C(i,j)%xedge(2)=D%R(i+1)
		C(i,j)%xedge(3)=D%Theta(j)
		C(i,j)%xedge(4)=D%Theta(j+1)
		C(i,j)%xedge2(1)=D%R(i)**2
		C(i,j)%xedge2(2)=D%R(i+1)**2
		C(i,j)%xedge2(3)=D%Theta(j)**2
		C(i,j)%xedge2(4)=D%Theta(j+1)**2
		C(i,j)%mass=C(i,j)%dens*C(i,j)%V
	enddo
	C(0,j)%T=0d0
	C(0,j)%Ni=0
	C(0,j)%EJv=0d0
	C(0,j)%iTD=0
	C(0,j)%alphaturb=alphaturb
	C(0,j)%xedge(1)=D%R(0)
	C(0,j)%xedge(2)=D%R(1)
	C(0,j)%xedge(3)=D%Theta(j)
	C(0,j)%xedge(4)=D%Theta(j+1)
	C(0,j)%xedge2(1)=D%R(0)**2
	C(0,j)%xedge2(2)=D%R(1)**2
	C(0,j)%xedge2(3)=D%Theta(j)**2
	C(0,j)%xedge2(4)=D%Theta(j+1)**2
	C(0,j)%mass=C(0,j)%dens*C(0,j)%V
	enddo

	do i=1,D%nTheta-1
		C(0,i)%dens=1d-60
		C(0,i)%dens0=1d-60
		C(0,i)%V=(4d0*pi/3d0)*(D%R(1)**3-D%R(0)**3)*
     &			(D%Theta(i)-D%Theta(i+1))*AU**3
		C(0,i)%mass=C(0,i)%V*C(0,i)%dens
	enddo

	do i=0,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%opacity_set=.false.
		enddo
	enddo

	shscale(0:D%nR)=1d0
	if(scalesh.eq.'FILE') then
		call regrid(shscalefile,D%R_av/AU,shscale,D%nR)
	else if(scalesh.eq.'VALUE') then
		shscale(0:D%nR)=shscalevalue
	else if(scalesh.eq.'RADIALPOW') then
		shscale(0:D%nR-1)=shscalevalue*(D%R_av(0:D%nR-1)/D%R_av(1))**(-shpow)
	else if(scalesh.eq.'POW') then
		shscale(0:D%nR-1)=1d0+(shscalevalue-1d0)*(D%R_av(0:D%nR-1)/AU)**(-shpow)
	else
		shscale(0:D%nR)=1d0
	endif

	do i=0,D%nR-1
		shscale(i)=shscale(i)*(1d0+rimscale*exp(-((D%R(1)-D%R(i))/rimwidth)**2))
	enddo

	MassTot=0d0
	D%Vtot=0d0

	if(denstype.eq.'SURFFILE') then
		inquire(file=densfile,exist=truefalse)
		if(.not.truefalse) then
			write(*,'("Surface density file not found")')
			write(9,'("Surface density file not found")')
			write(*,'("--------------------------------------------------------")')
			write(9,'("--------------------------------------------------------")')
			stop
		endif
		allocate(surfacedens(D%nR))
		call regridlog(densfile,D%R_av(1:D%nR)/AU,surfacedens,D%nR)
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
			C(i,j)%dens=surfacedens(i)*(D%R(i+1)**2-D%R(i)**2)/C(i,j)%V
		enddo
		enddo
		deallocate(surfacedens)

		if(gasdensfile.ne.' ') then
			inquire(file=gasdensfile,exist=truefalse)
			if(.not.truefalse) then
				write(*,'("Gas surface density file not found")')
				write(9,'("Gas surface density file not found")')
				write(*,'("--------------------------------------------------------")')
				write(9,'("--------------------------------------------------------")')
				stop
			endif
			allocate(surfacedens(D%nR))
			call regridlog(gasdensfile,D%R_av(1:D%nR)/AU,surfacedens,D%nR)
			do i=1,D%nR-1
			do j=1,D%nTheta-1
				C(i,j)%gasdens=surfacedens(i)*pi*(D%R(i+1)**2-D%R(i)**2)*AU**2/C(i,j)%V/real(D%nTheta-1)
				C(i,j)%gasdens=C(i,j)%gasdens/gas2dust
			enddo
			enddo
			deallocate(surfacedens)
		endif
	endif

	if(denstype.eq.'SHELLFILE') then
		inquire(file=densfile,exist=truefalse)
		if(.not.truefalse) then
			write(*,'("Shell density file not found")')
			write(9,'("Shell density file not found")')
			write(*,'("--------------------------------------------------------")')
			write(9,'("--------------------------------------------------------")')
			stop
		endif
		allocate(surfacedens(D%nR))
		call regridlog(densfile,D%R_av(1:D%nR)/AU,surfacedens,D%nR)
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
			C(i,j)%dens=surfacedens(i)
		enddo
		enddo
		deallocate(surfacedens)
		mdustscale=.false.
	endif
	
	if(denstype.eq.'PRODIMO') then
		inquire(file=densfile,exist=truefalse)
		if(.not.truefalse) then
			write(*,'("ProDiMo file not found")')
			write(9,'("ProDiMo file not found")')
			write(*,'("--------------------------------------------------------")')
			write(9,'("--------------------------------------------------------")')
			stop
		endif
		call ImportProdimo(densfile)
		mdustscale=.false.
	endif

	if(denstype.eq.'ZONES') then
		allocate(zonedens(nzones,ngrains,D%nR-1,D%nTheta-1)) ! first occurence (of 2)
		zonedens=0d0
		do iz=1,nzones
			tot=0d0
			do ii=1,ngrains
				if(Zone(iz)%inc_grain(ii)) then
					tot=tot+Zone(iz)%abun(ii)
				else
					Zone(iz)%abun(ii)=0d0
				endif
			enddo
			Zone(iz)%abun(1:ngrains)=Zone(iz)%abun(1:ngrains)/tot
			tot=0d0
			do i=1,D%nR-1
			do j=1,D%nTheta-1
				C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &					(D%Theta(j)-D%Theta(j+1))*AU**3
				if(D%R_av(i).ge.(Zone(iz)%Rin*AU).and.D%R_av(i).le.(Zone(iz)%Rout*AU)) then
					C(i,j)%alphaturb=Zone(iz)%alphaturb
					zonedens(iz,1:ngrains,i,j)=0d0
					if(Zone(iz)%sh.gt.0d0) then
						njj=10
						do jj=1,njj
							theta=D%thet(j)+(D%thet(j+1)-D%thet(j))*real(jj)/real(njj+1)
							r=D%R_av(i)*sin(theta)/AU
							z=D%R_av(i)*cos(theta)/AU
							hr=Zone(iz)%sh*(r/Zone(iz)%Rsh)**Zone(iz)%shpow
							f1=r**(-Zone(iz)%denspow)*exp(-(D%R_av(i)/(AU*Zone(iz)%Rexp))**(Zone(iz)%gamma_exp))
							if (Zone(iz)%roundtype.ne.'NONE') then
							   f1=f1*RoundOff(D%R_av(i)/AU,Zone(iz)%Rin+Zone(iz)%roundwidth,Zone(iz)%roundtype,
     &								Zone(iz)%roundindex,Zone(iz)%roundscalemin)
							endif
							f2=exp(-(z/hr)**2)
							do ii=1,ngrains
								if(Zone(iz)%inc_grain(ii)) then
									zonedens(iz,ii,i,j)=zonedens(iz,ii,i,j)+Zone(iz)%abun(ii)*f1*f2/hr/real(njj)
								else
									zonedens(iz,ii,i,j)=0d0
								endif
							enddo
						enddo
					else
c	this is a wedge zone!
						theta=D%thet(j)
						r=D%R_av(i)*sin(theta)/AU
						do ii=1,ngrains
							if(Zone(iz)%inc_grain(ii)) then
								if(j.eq.D%nTheta-1.or.(180d0*D%theta_av(j)/pi).gt.(-Zone(iz)%sh)) then
									zonedens(iz,ii,i,j)=zonedens(iz,ii,i,j)+Zone(iz)%abun(ii)*r**(-Zone(iz)%denspow-1d0)
								else
									zonedens(iz,ii,i,j)=0d0
								endif
							else
								zonedens(iz,ii,i,j)=0d0
							endif
						enddo
					endif
					if(D%R_av(i).ge.(Zone(iz)%pertR0*AU)) then
						do ii=1,ngrains
							zonedens(iz,ii,i,j)=zonedens(iz,ii,i,j)*(1d0+Zone(iz)%pertA*sin(2d0*pi*(D%R_av(i)/AU-Zone(iz)%pertR0)/Zone(iz)%pertR))
						enddo
					endif
					do ii=1,ngrains
						tot=tot+zonedens(iz,ii,i,j)*C(i,j)%V
					enddo
				else
					do ii=1,ngrains
						zonedens(iz,ii,i,j)=0d0
					enddo
				endif
			enddo
			enddo
			zonedens(iz,1:ngrains,1:D%nR-1,1:D%nTheta-1)=zonedens(iz,1:ngrains,1:D%nR-1,1:D%nTheta-1)*Zone(iz)%Mdust*Msun/tot
		enddo
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%dens=1d-60
			do iz=1,nzones
				do ii=1,ngrains
					if(Zone(iz)%inc_grain(ii)) C(i,j)%dens=C(i,j)%dens+zonedens(iz,ii,i,j)
				enddo
			enddo
		enddo
		enddo
		deallocate(zonedens)
	endif

	do i=1,D%nR-1
	do j=1,D%nTheta-1
		if(denstype.eq.'PASCUCCI') then
			r=D%R_av(i)*sin(D%theta_av(j))/AU
			z=D%R_av(i)*cos(D%theta_av(j))/AU
			rd=D%Rout/2d0
			zd=D%Rout/8d0
			f1=(r/rd)**(-1d0)
			hr=zd*(r/rd)**1.125d0
			f2=exp(-(pi/4d0)*(z/hr)**2)
			C(i,j)%dens=f1*f2
			scale=D%Mtot/(
     & (128d0/(17d0))*rd**2*zd*((D%Rout/rd)**(17d0/8d0)-(D%Rin/rd)**(17d0/8d0))
     & )/(AU**3)
			C(i,j)%dens=C(i,j)%dens*scale
			mdustscale=.false.
		else if(denstype.eq.'MIN') then
			r=D%R_av(i)*sin(D%theta_av(j))/AU
			z=D%R_av(i)*cos(D%theta_av(j))/AU
			hr=D%sh1AU*r**D%shpow
			f1=r**(-D%denspow)
			f2=exp(-(z/hr)**2)
c			f2=exp(-0.5*(z/hr)**2)
			scale=D%Mtot/(AU**3*2d0*pi*D%sh1AU*sqrt(pi)*(D%Rout-D%Rin))
			C(i,j)%dens=scale*f1*f2*r**(-D%shpow)
		else if(denstype.eq.'PARAMETERIZED') then
			njj=10
			C(i,j)%dens=0d0
			do jj=1,njj
				theta=D%thet(j)+(D%thet(j+1)-D%thet(j))*real(jj)/real(njj+1)
				r=D%R_av(i)*sin(theta)/AU
				z=D%R_av(i)*cos(theta)/AU
				hr=D%sh1AU*r**D%shpow
				hr=hr*shscale(i)
				f1=r**(-D%denspow)
				f2=exp(-(z/hr)**2)
				scale=D%Mtot/(AU**3*2d0*pi*D%sh1AU*sqrt(pi)*(D%Rout-D%Rin))
				C(i,j)%dens=C(i,j)%dens+scale*f1*f2/hr/real(njj)
			enddo
		else if(denstype.eq.'SHELL') then ! gijsexp: shell
			C(i,j)%dens=(AU/D%R_av(i))**(D%denspow) ! density
		else if(denstype.eq.'WEDGE') then ! wedge opening
			njj=10
			C(i,j)%dens=1d-200
			do jj=1,njj
				theta=D%thet(j)+(D%thet(j+1)-D%thet(j))*real(jj)/real(njj+1)
				if((180d0*(pi/2d0-theta)/pi).lt.wedgeopen) then
					C(i,j)%dens=C(i,j)%dens+(AU/D%R_av(i))**(D%denspow+1d0)
				endif
			enddo
		else if(denstype.eq.'POW') then
			C(i,j)%dens=(AU/D%R_av(i))**(D%denspow+1d0)
		else if(denstype.eq.'SIMILARITY') then
			C(i,j)%dens=(AU/D%R_av(i))**(D%denspow+1d0)*exp(-(D%R_av(i)/(AU*D%Rexp))**(D%gamma_exp))
		else if(denstype.eq.'DOUBLEPOW') then
			if(D%R_av(i).lt.(D%Rpow2*AU) .and. D%R_av(i).gt.(D%Rpow3*AU)) then
				C(i,j)%dens=(D%Rpow2*AU/D%R_av(i))**(D%denspow+1d0)
			else
				C(i,j)%dens=(D%Rpow2*AU/D%R_av(i))**(D%denspow2+1d0)
			endif
		else if(denstype.eq.'DOUBLEPOWSIM') then ! gijsexp
		   if(D%R_av(i).lt.(D%Rpow2*AU)) then
		      C(i,j)%dens=(D%Rpow2*AU/D%R_av(i))**(D%denspow+1d0)
		   else
			 C(i,j)%dens=(D%Rpow2*AU/D%R_av(i))**(D%denspow2+1d0)*exp(-(D%R_av(i)/(AU*D%Rexp))**(D%gamma_exp))
		   endif
		else if(denstype.eq.'MASSLOSS') then
			vexp=vexp2*(1.-(1.-(vexp1/vexp2)**(1./bexp))*D%Rstar/D%R_av(i))**bexp
			C(i,j)%dens=D%Mdot/(gas2dust*4d0*pi*D%R_av(i)**2*vexp)
		else if(denstype.eq.'PINTE') then
			r=D%R_av(i)*sin(D%theta_av(j))/AU
			z=D%R_av(i)*cos(D%theta_av(j))/AU
			hr=D%sh1AU*r**D%shpow
			f1=r**(-D%denspow)
			f2=exp(-(z/hr)**2)
			scale=D%Mtot/(AU**3*2d0*pi*D%sh1AU*sqrt(pi)*(D%Rout-D%Rin))
			if(r.lt.D%Rexp) then
				C(i,j)%dens=scale*f1*f2*r**(-D%shpow)
			else
				C(i,j)%dens=1d-60
			endif
		else if(denstype.eq.'INFALL') then
c See Dominik & Dullemond 2008, Eqs. 1 & 2
			mu=cos(D%theta_av(j))
			if(mu.gt.1d0) print*,mu,D%theta_av(j)
			call solve_mu0(D%R_av(i)/(AU*D%Rexp),mu,mu0)
			C(i,j)%dens=((D%Minfall/(4d0*pi*sqrt(6.67300d-8*D%Mstar*D%R_av(i)**3)))
     &					*(1d0+mu/mu0)**(-0.5d0)*(mu/mu0+2d0*mu0**2*D%Rexp*AU/D%R_av(i))**(-1d0))/gas2dust
		else if(denstype.eq.'MEIXNER') then
			r=D%R_av(i)/AU-timeshift*vexp1*0.210944839d-5
			if(r.gt.MeixRin) then
				C(i,j)%dens=((r/MeixRin)**(-MeixB*(1d0+MeixC*sin(D%theta_av(j))**MeixF*
     &					(dexp(-(r/MeixRsw)**MeixD+(MeixRin/MeixRsw)**MeixD)))))*
     &					(1d0+MeixA*(1d0-cos(D%theta_av(j)))**MeixF*
     &					(dexp(-(r/MeixRsw)**MeixE+(MeixRin/MeixRsw)**MeixE)))
				C(i,j)%dens=C(i,j)%dens*(r/(D%R_av(i)/AU))**2
			else
				C(i,j)%dens=1d-60
			endif
		endif
		C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
		C(i,j)%mass=C(i,j)%V*C(i,j)%dens
		MassTot=MassTot+C(i,j)%mass
		D%Vtot=D%Vtot+C(i,j)%V
	enddo
	enddo
	if(.not.mdustscale) D%Mtot=MassTot


	if(denstype.eq.'MIN') then
		scale=D%Mtot/(AU**2*2d0*pi*D%sh1AU*sqrt(pi)*(D%Rout**(-D%denspow+2d0)-D%Rin**(-D%denspow+2d0)))
		write(*,'("Error in the mass:    ",f10.3," %")') 100d0*(MassTot-D%Mtot)/D%Mtot
		write(*,'("Radial optical depth: ",f10.3," kappa")') scale*((D%Rin)**(-D%shpow)-(D%Rout)**(-D%shpow))/D%shpow
		write(9,'("Error in the mass:    ",f10.3," %")') 100d0*(MassTot-D%Mtot)/D%Mtot
		write(9,'("Radial optical depth: ",f10.3," kappa")') scale*((D%Rin)**(-D%shpow)-(D%Rout)**(-D%shpow))/D%shpow
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
	endif
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		C(i,j)%dens=C(i,j)%dens*D%Mtot/MassTot
		C(i,j)%mass=C(i,j)%V*C(i,j)%dens
		C(i,j)%dens0=C(i,j)%dens
	enddo
	enddo

	! Allocate arrays for grains
	if(.not.arraysallocated) then
		allocate(Grain(ngrains))

		ngrains2=0
		do ii=1,ngrains
			if(parttype(ii).eq.3.or.parttype(ii).eq.5) then
				call checkaggregates(partarg(ii),Grain(ii)%nopac)
				if(parttype(ii).eq.5) allocate(Grain(ii)%Topac(Grain(ii)%nopac))
				if(parttype(ii).eq.3) allocate(Grain(ii)%cryst(Grain(ii)%nopac))
			else if(parttype(ii).eq.8) then
				Grain(ii)%nopac=2
			else if(parttype(ii).eq.7.and.computepart_nT(ii).gt.0) then
				Grain(ii)%nopac=computepart_nT(ii)
				allocate(Grain(ii)%Topac(Grain(ii)%nopac))
			else
				Grain(ii)%nopac=1
			endif
			allocate(Grain(ii)%Kabs(Grain(ii)%nopac,nlam))
			allocate(Grain(ii)%Ksca(Grain(ii)%nopac,nlam))
			allocate(Grain(ii)%Kext(Grain(ii)%nopac,nlam))
			allocate(Grain(ii)%F(Grain(ii)%nopac,nlam))
			allocate(Grain(ii)%Kp(Grain(ii)%nopac,0:TMAX))
			allocate(Grain(ii)%Kpstar(Grain(ii)%nopac))
			allocate(Grain(ii)%Kpabsstar(Grain(ii)%nopac))
			allocate(Grain(ii)%shscale(D%nR))
			allocate(Grain(ii)%KabsL(Grain(ii)%nopac))
			allocate(Grain(ii)%KextL(Grain(ii)%nopac))
			allocate(Grain(ii)%rho(Grain(ii)%nopac))
			allocate(Grain(ii)%rscale(Grain(ii)%nopac))
			allocate(Grain(ii)%mscale(Grain(ii)%nopac))
			if(Grain(ii)%nopac.gt.ngrains2) ngrains2=Grain(ii)%nopac
		enddo
		do i=0,D%nR
		do j=0,D%nTheta
			allocate(C(i,j)%w(ngrains))
			allocate(C(i,j)%w0(ngrains))
			allocate(C(i,j)%wopac(ngrains,ngrains2))
			C(i,j)%wopac(1:ngrains,1)=1d0
			if(ngrains2.gt.1) C(i,j)%wopac(1:ngrains,2:ngrains2)=0d0
		enddo
		enddo
	endif

	!  Calculate edges of grain size bins
	if(gsd) call gsd_edges(rgrain(1:ngrains),rgrain_edges(1:ngrains+1))

	! Set grain properties
	do ii=1,ngrains
		Grain(ii)%TdesA=TdesA(ii)
		Grain(ii)%TdesB=TdesB(ii)
		Grain(ii)%force_vert_gf=force_vert_gf(ii)
		Grain(ii)%tdes_fast=tdes_fast(ii)
		Grain(ii)%rv=rgrain(ii)
		Grain(ii)%rho(1:Grain(ii)%nopac)=rhograin(ii)
		Grain(ii)%rscale(1:Grain(ii)%nopac)=1d0
		Grain(ii)%mscale(1:Grain(ii)%nopac)=1d0
		Grain(ii)%Rcryst=radmix(ii)
		Grain(ii)%Tcryst=tmix(ii)
		Grain(ii)%powcryst=powmix(ii)
		if (gsd) Grain(ii)%rvmin=rgrain_edges(ii)
		if (gsd) Grain(ii)%rvmax=rgrain_edges(ii+1)
	enddo
	
	! Read in the particle files / opacities
	nqhp=0
	do ii=1,ngrains
		Grain(ii)%parttype=parttype(ii)
		if(parttype(ii).eq.1) then
			if(scattype.eq.'RAYLEIGH') then
				call readparticle(partarg(ii),Grain(ii),.true.,asym(ii),asym2(ii),wasym2(ii),Pmax(ii),1)
			else
				call readparticle(partarg(ii),Grain(ii),.false.,asym(ii),asym2(ii),wasym2(ii),Pmax(ii),1)
			endif
		else if(parttype(ii).eq.6) then
			if(scattype.eq.'RAYLEIGH') then
				call ReadParticleFits(partarg(ii),Grain(ii),.true.,asym(ii),asym2(ii),wasym2(ii),Pmax(ii),1)
			else
				call ReadParticleFits(partarg(ii),Grain(ii),.false.,asym(ii),asym2(ii),wasym2(ii),Pmax(ii),1)
			endif
		else if(parttype(ii).eq.2) then
			call readopacity(partarg(ii),Grain(ii))
		else if(parttype(ii).eq.3) then
			call MixedAggregates(partarg(ii),Grain(ii))
		else if(parttype(ii).eq.4) then
			nqhp=nqhp+1
			Grain(ii)%qhpnr=nqhp
			call readQHP(partarg(ii),Grain(ii))
		else if(parttype(ii).eq.5) then
			call readTopac(partarg(ii),Grain(ii))
		else if(parttype(ii).eq.7) then
			if(computepart_standard(ii).eq.'DIANA') then
				computepart_abun(ii,1)=1d0-computepart_fcarbon(ii)
				computepart_abun(ii,2)=computepart_fcarbon(ii)
				computepart_norm_abun(ii)=-1d0
				computepart_nT(ii)=0
				partarg(ii)=computepart_standard(ii)
			else if(computepart_standard(ii).eq.'ASTROSIL') then
				computepart_norm_abun(ii)=-1d0
				computepart_nT(ii)=0
				partarg(ii)=computepart_standard(ii)
			endif
			if(computepart_nT(ii).eq.0) then
				call ComputePart(Grain(ii),ii,partarg(ii),computepart_amin(ii),computepart_amax(ii)
     &							,computepart_apow(ii),computepart_fmax(ii),computepart_blend(ii)
     &							,computepart_porosity(ii),computepart_abun(ii,:),1,computepart_norm_abun(ii)
     &							,computepart_standard(ii),computepart_nsubgrains(ii))
			else
				do i=1,Grain(ii)%nopac
					jj=i
					T=1d10
					do j=1,Grain(ii)%nopac
						if(computepart_T(ii,j).lt.T) then
							if(i.eq.1) then
								jj=j
								T=computepart_T(ii,j)
							else if(computepart_T(ii,j).gt.Grain(ii)%Topac(i-1)) then
								jj=j
								T=computepart_T(ii,j)
							endif
						endif
					enddo
					call ComputePart(Grain(ii),ii,computepart_Tfile(ii,jj),computepart_amin(ii),computepart_amax(ii)
     &							,computepart_apow(ii),computepart_fmax(ii),computepart_blend(ii)
     &							,computepart_porosity(ii),computepart_abun(ii,:),i,computepart_norm_abun(ii)
     &							,computepart_standard(ii),computepart_nsubgrains(ii))
					Grain(ii)%Topac(i)=computepart_T(ii,jj)
				enddo
				do i=0,D%nR
				do j=0,D%nTheta
					C(i,j)%wopac(ii,Grain(ii)%nopac)=1d0
					do jj=1,Grain(ii)%nopac-1
						C(i,j)%wopac(ii,jj)=0d0
					enddo
				enddo
				enddo
			endif
		else if(parttype(ii).eq.8) then
			nqhp=nqhp+1
			Grain(ii)%qhpnr=nqhp
			call ComputePAH(Grain(ii),computepart_amin(ii),computepart_amax(ii),computepart_apow(ii))
		endif

		! Set scale height scaling of dust
		Grain(ii)%settle=settle(ii)
		if(settlefile(ii).ne.' ') then
			call regrid(settlefile(ii),D%R_av(1:D%nR)/AU,Grain(ii)%shscale(1:D%nR),D%nR)
		else
			Grain(ii)%shscale(1:D%nR-1)=part_shscale(ii)
		endif
		Grain(ii)%shtype=shtype(ii)
		do i=1,D%nR-1
			if(D%R_av(i).ge.psettR0*AU) then
				Grain(ii)%shscale(i)=Grain(ii)%shscale(i)*
     & (1d0-((D%R_av(i)-psettR0*AU)/D%R_av(i))**psettpow*(1d0-psettphi0**(log10(Grain(ii)%rv)+6d0)))
c				if(Grain(ii)%shscale(i).lt.0.2d0) Grain(ii)%shscale(i)=0.2d0
			endif
		enddo

		! Check for negative opacities and set them to 0
		do i=1,nlam
			do iopac=1,Grain(ii)%nopac
				if(Grain(ii)%Kabs(iopac,i).lt.0d0) Grain(ii)%Kabs(iopac,i)=0d0
				if(Grain(ii)%Kext(iopac,i).lt.0d0) Grain(ii)%Kext(iopac,i)=0d0
				if(Grain(ii)%Ksca(iopac,i).lt.0d0) Grain(ii)%Ksca(iopac,i)=0d0
				Grain(ii)%Kext(iopac,i)=Grain(ii)%Kabs(iopac,i)+Grain(ii)%Ksca(iopac,i)
			enddo
		enddo

		if(parttype(ii).ne.6.and.parttype(ii).ne.7) then
			Grain(ii)%dust_moment1=Grain(ii)%rv*1d4
			Grain(ii)%dust_moment2=(Grain(ii)%rv*1d4)**2
			Grain(ii)%dust_moment3=(Grain(ii)%rv*1d4)**3
		endif
		if(parttype(ii).eq.7.and.computepart_nT(ii).ne.0) Grain(ii)%parttype=5
	enddo
	

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')	

	! Set composition in every grid cell (from keywords, mrn/gsd or file)
	if(denscomposition) then
		write(*,'("Reading composition from density file")')
		write(9,'("Reading composition from density file")')
		call readstruct(compositionfile,(/'SKIP   ','COMP   '/),2,0,.false.)
		do j=1,D%nTheta-1
			C(0,j)%w(1:ngrains)=C(1,j)%w(1:ngrains)
		enddo
	else if(compositionfile.eq.' ') then
		do i=0,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%w(1:ngrains)=warg(1:ngrains)/sum(warg(1:ngrains))
		enddo
		enddo
	else
		call readcomposition(compositionfile)
		do j=1,D%nTheta-1
			C(0,j)%w(1:ngrains)=C(1,j)%w(1:ngrains)
		enddo
	endif

	if(denstype.eq.'MEIXNER'.and.MeixG.ne.0d0) then
		tot=mrn_index
		do i=0,D%nR-1
		do j=1,D%nTheta-1
			mrn_index=MeixG+(tot-MeixG)*(1d0-cos(D%theta_av(j)))**MeixF
			call gsd_MRN(rgrain(1:ngrains),warg(1:ngrains))	! GijsExp
			C(i,j)%w(1:ngrains)=warg(1:ngrains)/sum(warg(1:ngrains))
		enddo
		enddo
	endif
	
	! Radial composition gradient (powerlaw) 
	do ii=1,npow
	   do i=1,D%nR-1
	      r=D%R_av(i)/AU
	      do j=1,D%nTheta-1
		 if(r.gt.powrad0(ii)) then
		    f1=(powrad0(ii)/r)**powslope(ii)
		 else
		    f1=powinner(ii)
		 endif
		 C(i,j)%w(powclose(ii))=C(i,j)%w(powclose(ii))+C(i,j)%w(powfar(ii))*f1
		 C(i,j)%w(powfar(ii))=C(i,j)%w(powfar(ii))*(1d0-f1)
	      enddo
	   enddo
	enddo

	! Set composition for temperature-dependent species
	do ii=1,ngrains
	   if (Grain(ii)%nopac .gt. 1) then
	    do i=0,D%nR-1
		do j=1,D%nTheta-1
		   C(i,j)%wopac(ii,1:Grain(ii)%nopac)=1d0/real(Grain(ii)%nopac)
		enddo
		enddo
	   endif
	enddo

	if(denstype.eq.'ZONES') then
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%w(1:ngrains)=0d0
		enddo
		enddo
		do iz=1,nzones
			if(Zone(iz)%sizedis) then
				mrn_rmin=Zone(iz)%a_min*1d-4
				mrn_rmax=Zone(iz)%a_max*1d-4
				mrn_index=Zone(iz)%a_pow
				allocate(w(ngrains))
				allocate(rtemp(ngrains))
				j=0
				do ii=1,ngrains
					if(Zone(iz)%inc_grain(ii).and..not.Grain(ii)%qhp) then
						j=j+1
						rtemp(j)=Grain(ii)%rv
					endif
				enddo
				mrn_ngrains0=mrn_ngrains
				mrn_ngrains=j
				w(1:ngrains)=0d0
				do ii=j+1,ngrains
					rtemp(ii)=(10d0+real(ii))*mrn_rmax
				enddo
				call gsd_MRN(rtemp(1:ngrains),w(1:ngrains))
				mrn_ngrains=mrn_ngrains0
				j=0
				tot=0d0
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
						tot=tot+warg(ii)
						print*,ii,warg(ii)
					endif
				enddo
				do ii=1,ngrains
					if(Zone(iz)%inc_grain(ii).and..not.Grain(ii)%qhp) then
						j=j+1
						Zone(iz)%abun(ii)=w(j)*(1d0-Zone(iz)%fPAH)
					else if(Grain(ii)%qhp) then
						Zone(iz)%abun(ii)=Zone(iz)%fPAH*warg(ii)/tot
					else
						Zone(iz)%abun(ii)=0d0
					endif
				enddo
				deallocate(w)
				deallocate(rtemp)
			else
				do ii=1,ngrains
					if(Zone(iz)%inc_grain(ii)) then
						Zone(iz)%abun(ii)=ZoneTemp(iz)%abun(ii)
					else
						Zone(iz)%abun(ii)=0d0
					endif
				enddo
			endif
			do i=1,D%nR-1
				if(D%R_av(i).ge.(Zone(iz)%Rin*AU).and.D%R_av(i).le.(Zone(iz)%Rout*AU)) then
					do j=1,D%nTheta-1
						do ii=1,ngrains
							if(Zone(iz)%inc_grain(ii)) C(i,j)%w(ii)=Zone(iz)%abun(ii)
						enddo
					enddo
				endif
			enddo
		enddo
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			tot=sum(C(i,j)%w(1:ngrains))
			if(tot.lt.1d-50) C(i,j)%w(1:ngrains)=1d0/real(ngrains)
		enddo
		enddo

		allocate(zonedens(nzones,ngrains,D%nR-1,D%nTheta-1)) ! second occurence (of 2)
		zonedens=0d0
		D%Mtot=0d0
		do iz=1,nzones
			tot=0d0
			do ii=1,ngrains
				if(Zone(iz)%inc_grain(ii)) tot=tot+Zone(iz)%abun(ii)
			enddo
			Zone(iz)%abun(1:ngrains)=Zone(iz)%abun(1:ngrains)/tot
			tot=0d0
			do i=1,D%nR-1
			do j=1,D%nTheta-1
				C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &					(D%Theta(j)-D%Theta(j+1))*AU**3
				if(D%R_av(i).ge.(Zone(iz)%Rin*AU).and.D%R_av(i).le.(Zone(iz)%Rout*AU)) then
					C(i,j)%alphaturb=Zone(iz)%alphaturb
					zonedens(iz,1:ngrains,i,j)=0d0
					if(Zone(iz)%sh.gt.0d0) then
						njj=10
						do jj=1,njj
							theta=D%thet(j)+(D%thet(j+1)-D%thet(j))*real(jj)/real(njj+1)
							r=D%R_av(i)*sin(theta)/AU
							z=D%R_av(i)*cos(theta)/AU
							hr=Zone(iz)%sh*(r/Zone(iz)%Rsh)**Zone(iz)%shpow
							f1=r**(-Zone(iz)%denspow)*exp(-(D%R_av(i)/(AU*Zone(iz)%Rexp))**(Zone(iz)%gamma_exp))
							if (Zone(iz)%roundtype.ne.'NONE') then
							   f1=f1*RoundOff(D%R_av(i)/AU,Zone(iz)%Rin+Zone(iz)%roundwidth,Zone(iz)%roundtype,
     &								Zone(iz)%roundindex,Zone(iz)%roundscalemin)
							endif
							f2=exp(-(z/hr)**2)
							do ii=1,ngrains
								if(Zone(iz)%inc_grain(ii)) then
									zonedens(iz,ii,i,j)=zonedens(iz,ii,i,j)+Zone(iz)%abun(ii)*f1*f2/hr/real(njj)
								else
									zonedens(iz,ii,i,j)=0d0
								endif
							enddo
						enddo
					else
c	this is a wedge zone!
						r=D%R_av(i)*sin(theta)/AU
						do ii=1,ngrains
							if(Zone(iz)%inc_grain(ii)) then
								if(j.eq.D%nTheta-1.or.(180d0*D%theta_av(j)/pi).gt.(-Zone(iz)%sh)) then
									zonedens(iz,ii,i,j)=zonedens(iz,ii,i,j)+Zone(iz)%abun(ii)*r**(-Zone(iz)%denspow-1d0)
								else
									zonedens(iz,ii,i,j)=0d0
								endif
							else
								zonedens(iz,ii,i,j)=0d0
							endif
						enddo
					endif
					if(D%R_av(i).ge.(Zone(iz)%pertR0*AU)) then
						do ii=1,ngrains
							zonedens(iz,ii,i,j)=zonedens(iz,ii,i,j)*(1d0+Zone(iz)%pertA*sin(2d0*pi*(D%R_av(i)/AU-Zone(iz)%pertR0)/Zone(iz)%pertR))
						enddo
					endif
					do ii=1,ngrains
						tot=tot+zonedens(iz,ii,i,j)*C(i,j)%V
					enddo
				else
					do ii=1,ngrains
						zonedens(iz,ii,i,j)=0d0
					enddo
				endif
			enddo
			enddo
			if(Zone(iz)%maxtauV.lt.0d0) then
				zonedens(iz,1:ngrains,1:D%nR-1,1:D%nTheta-1)=zonedens(iz,1:ngrains,1:D%nR-1,1:D%nTheta-1)*Zone(iz)%Mdust*Msun/tot
			else
				do i=1,nlam-1
					if(lam(i).le.0.55.and.lam(i+1).gt.0.55) then
						jj=i
						wl1=(lam(i+1)-0.55)/(lam(i+1)-lam(i))
						wl2=1d0-wl1
						exit
					endif
				enddo
				maxtauV=0d0
				do j=1,D%nTheta-1
					tau=0d0
					do i=1,D%nR-1
						do ii=1,ngrains
							do iopac=1,Grain(ii)%nopac
								tau=tau+zonedens(iz,ii,i,j)*(wl1*Grain(ii)%Kext(iopac,jj)+wl2*Grain(ii)%Kext(iopac,jj+1))
     &										*(D%R(i+1)-D%R(i))*AU*C(i,j)%wopac(ii,iopac)
							enddo
						enddo
					enddo
					if(tau.gt.maxtauV) maxtauV=tau
				enddo
				zonedens(iz,1:ngrains,1:D%nR-1,1:D%nTheta-1)=
     %	zonedens(iz,1:ngrains,1:D%nR-1,1:D%nTheta-1)*Zone(iz)%maxtauV/maxtauV
				Zone(iz)%Mdust=(tot*Zone(iz)%maxtauV/maxtauV)/Msun
			endif
			D%Mtot=D%Mtot+Zone(iz)%Mdust
		enddo
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%dens=1d-60
			C(i,j)%w(1:ngrains)=0d0
			do iz=1,nzones
				do ii=1,ngrains
					if(Zone(iz)%inc_grain(ii)) then
						C(i,j)%dens=C(i,j)%dens+zonedens(iz,ii,i,j)
						C(i,j)%w(ii)=C(i,j)%w(ii)+zonedens(iz,ii,i,j)
					endif
				enddo
			enddo
			tot=sum(C(i,j)%w(1:ngrains))
			if(tot.lt.1d-50) then
				C(i,j)%w(1:ngrains)=1d0/real(ngrains)
			else
				C(i,j)%w(1:ngrains)=C(i,j)%w(1:ngrains)/tot
			endif
			C(i,j)%dens0=C(i,j)%dens
			C(i,j)%mass=C(i,j)%dens*C(i,j)%V
			C(i,j)%w0(1:ngrains)=C(i,j)%w(1:ngrains)
		enddo
		enddo
		deallocate(zonedens)
		if(compositionfile.ne.' ') then
			write(*,'("Reading composition from file")')
			write(9,'("Reading composition from file")')
			call readcomposition(compositionfile)
		endif
	endif

	MassTot0=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		MassTot0=MassTot0+C(i,j)%dens*C(i,j)%V
	enddo
	enddo

	do ii=1,ngrains
		Grain(ii)%maxrad=maxrad(ii)
		Grain(ii)%minrad=minrad(ii)
		Grain(ii)%maxtheta=maxtheta(ii)
		Grain(ii)%mintheta=mintheta(ii)
		Grain(ii)%shaperad=shaperad(ii)
		Grain(ii)%roundtype=roundtype(ii)
		Grain(ii)%roundpow=roundpow(ii)
		Grain(ii)%roundwidth=roundwidth(ii)
		Grain(ii)%roundtype=roundtype(ii)
	enddo

	MassTot=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		MassTot=MassTot+C(i,j)%dens*C(i,j)%V
	enddo
	enddo
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		C(i,j)%dens0=C(i,j)%dens0*MassTot0/MassTot
		C(i,j)%dens=C(i,j)%dens*MassTot0/MassTot
	enddo
	enddo


	allocate(spec(nlam))
	do ii=1,ngrains
	if(scat_how.eq.0) then
		do i=1,nlam
			Grain(ii)%Ksca(1:Grain(ii)%nopac,i)=0d0
			Grain(ii)%Kext(1:Grain(ii)%nopac,i)=Grain(ii)%Kabs(1:Grain(ii)%nopac,i)
		enddo
	else if(scat_how.eq.1) then
		do i=1,nlam
			do j=1,180
				Grain(ii)%F(1:Grain(ii)%nopac,i)%F11(j)=1d0
			enddo
		enddo
	endif
	do i=1,nlam
		Grain(ii)%F(1:Grain(ii)%nopac,i)%IF11=0d0
		Grain(ii)%F(1:Grain(ii)%nopac,i)%IF12=0d0
		do j=1,180
			thet=pi*(real(j)-0.5d0)/180d0
			Grain(ii)%F(1:Grain(ii)%nopac,i)%IF11=Grain(ii)%F(1:Grain(ii)%nopac,i)%IF11+pi*sin(thet)
     &			*Grain(ii)%F(1:Grain(ii)%nopac,i)%F11(j)/180d0
			Grain(ii)%F(1:Grain(ii)%nopac,i)%IF12=Grain(ii)%F(1:Grain(ii)%nopac,i)%IF12+pi*sin(thet)
     &			*Grain(ii)%F(1:Grain(ii)%nopac,i)%F12(j)/180d0
		enddo
	enddo
	do i=0,360
		phi=pi*real(i-1)/179.5d0
		cos2phi(i)=cos(2d0*phi)
		sin2phi(i)=sin(2d0*phi)
	enddo


	do j=0,TMAX
		do iopac=1,Grain(ii)%nopac
			spec(1:nlam)=BB(1:nlam,j)*Grain(ii)%Kabs(iopac,1:nlam)
			call integrate(spec,tot)
			Grain(ii)%Kp(iopac,j)=tot
		enddo
	enddo
	do iopac=1,Grain(ii)%nopac
		spec(1:nlam)=D%Fstar(1:nlam)*Grain(ii)%Kext(iopac,1:nlam)
		call integrate(spec,Grain(ii)%Kpstar(iopac))
		Grain(ii)%Kpstar(iopac)=Grain(ii)%Kpstar(iopac)/D%Lstar

		spec(1:nlam)=D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)
		call integrate(spec,Grain(ii)%Kpabsstar(iopac))
		Grain(ii)%Kpabsstar(iopac)=Grain(ii)%Kpabsstar(iopac)/D%Lstar
	enddo
	enddo

	deallocate(spec)
	
	scattering=.false.
	if(scat_how.ne.0) then
	do i=1,nlam
		do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			if(Grain(ii)%Ksca(iopac,i).ne.0d0) scattering=.true.
		enddo
		enddo
	enddo
	endif
	

	if(storescatt.and..not.arraysallocated) then
		do i=0,D%nR
			do j=0,D%nTheta
				allocate(C(i,j)%scattfield(1,0:nlam,2))
			enddo
		enddo
		allocate(scattcomputed(nlam))
		allocate(nscattcomputed(nlam))
		do i=1,nlam
			scattcomputed(i)=.false.
			nscattcomputed(i)=0
		enddo
	endif
	do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%w0(1:ngrains)=C(i,j)%w(1:ngrains)
		enddo
	enddo

	if(.not.arraysallocated.and.(.not.tcontact.or.tdes_iter.or.use_qhp)) then
	do i=0,D%nR
	do j=0,D%nTheta
		allocate(C(i,j)%TP(ngrains))
		allocate(C(i,j)%EJvP(ngrains))
		allocate(C(i,j)%EabsP(ngrains))
	enddo
	enddo
	endif
	if(use_qhp) then
		do i=0,D%nR
		do j=0,D%nTheta
			allocate(C(i,j)%QHP(nqhp,nlam))
			allocate(C(i,j)%tdistr(nqhp,NTQHP))
			allocate(C(i,j)%Tqhp(nqhp))
			allocate(C(i,j)%EJvQHP(nqhp))
			allocate(C(i,j)%EabsQHP(nqhp))
		enddo
		enddo
		computeLRF=.true.
	endif
	if(computeLRF) then
		do i=0,D%nR
		do j=0,D%nTheta
			allocate(C(i,j)%LRF(nlam))
			allocate(C(i,j)%nLRF(nlam))
		enddo
		enddo
	endif

	do i=1,D%nR-1
	do j=1,D%nTheta-1
		call CheckMinimumDensity(i,j)
	enddo
	enddo

	! allocate gap arrays
	if(Nphot.ne.0.and.startiter.eq.' ') then
		allocate(D%gap(ngap))
		allocate(D%gap1(ngap))
		allocate(D%gap2(ngap))
		allocate(D%gapshape(ngap))
		allocate(D%gaproundtype(ngap))
		allocate(D%gaproundpow(ngap))
		D%ngap=ngap
		do k=1,ngap
			D%gap(k)=gap(k)
			D%gap1(k)=gap1(k)
			D%gap2(k)=gap2(k)
			D%gapshape(k)=gapshape(k)
			D%gaproundtype(k)=gaproundtype(k)
			D%gaproundpow(k)=gaproundpow(k)
		enddo
		call MakeGaps()
	endif

	! allocate accretion profile array
	if((deadzone.or.gravstable.or.D%Rexp.lt.1d10).and.Nphot.ne.0.and.startiter.eq.' ') then
	   allocate(D%MdotR(D%nR))
	endif

	! allocate deadzone location array
	if(deadzone) allocate(D%MPdead(D%nR))
	
	MassTot=0d0
	do i=0,D%nR-1
	do j=0,D%nTheta-1
		C(i,j)%w0(1:ngrains)=C(i,j)%w(1:ngrains)
		C(i,j)%dens0=C(i,j)%dens
		if(denstype.ne.'PRODIMO'.and.(denstype.ne.'SURFFILE'.or.gasdensfile.eq.' ')) C(i,j)%gasdens=C(i,j)%dens
		C(i,j)%gasfrac=0d0
		MassTot=MassTot+C(i,j)%mass
	enddo
	enddo

	if(viscous) then
		do i=0,D%nR-1
		do j=0,D%nTheta-1
			allocate(C(i,j)%EviscDirect(ngrains))
			C(i,j)%EviscDirect(1:ngrains)=0d0
		enddo
		enddo
	endif

c-----------------------------------------------------------------------
c Calculate the disk density structure?
c-----------------------------------------------------------------------
	
	if(Nphot.gt.0) then
	if(denstype.eq.'POW'.or.denstype.eq.'DOUBLEPOW'.or.
     &     denstype.eq.'SURFFILE'.or.denstype.eq.'SIMILARITY'.or.
     &     denstype.eq.'DOUBLEPOWSIM'.or.
     &     struct_iter.or.mpset.or.scset.or.fixmpset.or.gsd.or.nzones.ne.0) then

           ! Set the disk temperature to the optically thin temperature
	   do i=1,D%nR-1
	      do j=1,D%nTheta-1
		 C(i,j)%dens0=C(i,j)%dens
		 if(denstype.ne.'PRODIMO'.and.(denstype.ne.'SURFFILE'.or.gasdensfile.eq.' ')) C(i,j)%gasdens=C(i,j)%dens0
		 C(i,j)%gasfrac=0d0
		 C(i,j)%T=D%Tstar*sqrt(D%R_av(1)/D%R_av(i))/10d0
		 if(.not.tcontact) then
		    C(i,j)%TP(1:ngrains)=C(i,j)%T
		 endif
	      enddo
	   enddo

		call BackWarming(1d0)
		call OpticallyThin(.false.)
		dosmooth=.false.
		call TempAverage(f_weight)
		if(deadzone) call SetIsoTdensity()
		call MakeDeadZone(.false.)
		if(gsd) call GrainsizeDistribution() ! GijsExp
		call DiskStructure()
		if (thgridrefine) call RegridTheta(D%nTheta/4)
		call MakeDeadZone(.true.)
		if(gsd) call GrainsizeDistribution() ! GijsExp
		call DiskStructure()
		if (thgridrefine) call RegridTheta(D%nTheta/4)
		call DiskStructure()
		dosmooth=.true.
	endif
	endif

	if(use_qhp.and.Nphot.gt.0) then
c		call PAHMCMax(0)
		call Stochastic(0)
		call destroyQHP(TdesQHP)
	endif
	
	if(Nphot.gt.0) then
		MassTot0=0d0
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			MassTot0=MassTot0+C(i,j)%dens*C(i,j)%V
		enddo
		enddo

		call DestroyDustR()

		MassTot=0d0
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			MassTot=MassTot+C(i,j)%dens*C(i,j)%V
		enddo
		enddo
		if(.not.mdustscale) MassTot=MassTot0
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%dens0=C(i,j)%dens0*MassTot0/MassTot
			C(i,j)%dens=C(i,j)%dens*MassTot0/MassTot
			call CheckMinimumDensity(i,j)
		enddo
		enddo
	endif

	if(denstype.eq.'PASCUCCI') then
	write(*,'("Scaling density structure to tau550: ",f10.3)') tau550
	write(9,'("Scaling density structure to tau550: ",f10.3)') tau550
	do i=1,nlam-1
		if(0.55.ge.lam(i).and.0.55.le.lam(i+1)) then
			ilam1=i
			ilam2=i+1
			wl1=(lam(i+1)-0.55)/(lam(i+1)-lam(i))
			wl2=(0.55-lam(i))/(lam(i+1)-lam(i))
		endif
	enddo
	Kext550=0d0
	do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			Kext550=Kext550+Grain(ii)%Kext(iopac,ilam1)*wl1+Grain(ii)%Kext(iopac,ilam2)*wl2
		enddo
	enddo
	rd=D%Rout/2d0
	zd=D%Rout/8d0
	scale=tau550/(rd*log(D%Rout/D%Rin)*Kext550)
	MassTot=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		r=D%R_av(i)*sin(D%theta_av(j))/AU
		z=D%R_av(i)*cos(D%theta_av(j))/AU
		f1=(r/rd)**(-1d0)
		hr=zd*(r/rd)**1.125d0
		f2=exp(-(pi/4d0)*(z/hr)**2)
		C(i,j)%dens=f1*f2*scale/(AU)
		C(i,j)%mass=C(i,j)%V*C(i,j)%dens
		MassTot=MassTot+C(i,j)%mass
		C(i,j)%dens0=C(i,j)%dens
	enddo
	enddo
	D%Mtot=MassTot
	endif
	
	call BackWarming(1d0)

	if(Nphot.gt.0) then
	call OpticallyThin(.true.)
	do i=1,ngrains
		if(Grain(i)%parttype.eq.3) call MakeMixAggregates(i)
	enddo
	endif
	do ii=1,ngrains
		allocate(Grain(ii)%g(Grain(ii)%nopac,nlam))
		do iopac=1,Grain(ii)%nopac
		do i=1,nlam
			Grain(ii)%g(iopac,i)=0d0
			tot=0d0
			do ia=1,180
				Grain(ii)%g(iopac,i)=Grain(ii)%g(iopac,i)+Grain(ii)%F(iopac,i)%F11(ia)*cos(pi*(real(ia)-0.5)/180)
     &						*sin(pi*(real(ia)-0.5)/180d0)
				tot=tot+Grain(ii)%F(iopac,i)%F11(ia)*sin(pi*(real(ia)-0.5)/180d0)
			enddo
			Grain(ii)%g(iopac,i)=Grain(ii)%g(iopac,i)/tot
		enddo
		enddo
		Grain(ii)%material=material(ii)
	enddo
	
	RmaxRefine=D%Rout
	tauincrease=1d30
	do i=0,D%nR-1
	do j=1,D%nTheta-1
		C(i,j)%KextLRF=0d0
		do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			C(i,j)%KextLRF=C(i,j)%KextLRF+C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)*Grain(ii)%Kpstar(iopac)
		enddo
		enddo
		C(i,j)%ILRF=1d0
		if( viscous ) C(i,j)%EviscDirect=0.0d0
	enddo
	enddo

	do ii=1,ngrains
	do iopac=1,Grain(ii)%nopac
	do i=1,nlam
		if(Grain(ii)%Kabs(iopac,i).eq.0d0) then
			Grain(ii)%Kabs(iopac,i)=1d-8
			Grain(ii)%Kext(iopac,i)=Grain(ii)%Kabs(iopac,i)+Grain(ii)%Ksca(iopac,i)
		endif
	enddo
	enddo
	enddo

	iter0=0
	if(startiter.ne.' ') read(startiter,*) iter0

c Add a possible infalling cloud
	if(Nphot.gt.0) then
	do ii=1,ngrains
		if(Grain(ii)%shtype.eq.'INFALL') then
			do i=0,D%nR-1
				do j=1,D%nTheta-1
					mu=cos(D%theta_av(j))
					call solve_mu0(D%R_av(i)/(AU*D%Rexp),mu,mu0)
					if(mu0.lt.D%mu0max) then
						rho=((D%Minfall/(4d0*pi*sqrt(6.67300d-8*D%Mstar*D%R_av(i)**3)))
     &					*(1d0+mu/mu0)**(-0.5d0)*(mu/mu0+2d0*mu0**2*D%Rexp*AU/D%R_av(i))**(-1d0))
					else
						rho=1d-60
					endif
					rho=rho/gas2dust
					C(i,j)%dens=C(i,j)%dens*(1d0-C(i,j)%w(ii))
					C(i,j)%dens0=C(i,j)%dens0*(1d0-C(i,j)%w0(ii))
c					C(i,j)%gasdens=C(i,j)%gasdens*(1d0-C(i,j)%w0(ii))
					C(i,j)%w=C(i,j)%w*C(i,j)%dens
					C(i,j)%w0=C(i,j)%w0*C(i,j)%dens0
					C(i,j)%w(ii)=rho
					C(i,j)%w0(ii)=rho
					C(i,j)%dens=C(i,j)%dens+rho
					C(i,j)%dens0=C(i,j)%dens0+rho
					call CheckMinimumDensity(i,j)
					tot=sum(C(i,j)%w(1:ngrains))
					if(tot.gt.1d-50) then
						C(i,j)%w=C(i,j)%w/tot
					else
						C(i,j)%w(1:ngrains)=warg(1:ngrains)/sum(warg(1:ngrains))
					endif
					tot=sum(C(i,j)%w0(1:ngrains))
					if(tot.gt.1d-50) then
						C(i,j)%w0=C(i,j)%w0/tot
					else
						C(i,j)%w0(1:ngrains)=warg(1:ngrains)/sum(warg(1:ngrains))
					endif
				enddo
			enddo
		endif
	enddo
	endif

	write(file,'(a,"betas.dat")') outdir(1:len_trim(outdir))
	open(unit=90,file=file,RECL=6000)
	write(90,*) ngrains
	allocate(spec(nlam))
	do ii=1,ngrains
	do iopac=1,Grain(ii)%nopac
		spec(1:nlam)=Grain(ii)%Kabs(iopac,1:nlam)+Grain(ii)%Ksca(iopac,1:nlam)*(1d0-Grain(ii)%g(iopac,1:nlam))
		call integrate(spec*D%Fstar,tot)
		call integrate(D%Fstar,tot2)
		tot=tot/tot2
		tot=tot*D%Lstar*Lsun/Luminosity(5777d0,Rsun)
		tot=(tot*1d4/(4d0*pi*clight))/(6.67300d-8*D%Mstar)
		write(90,*) Grain(ii)%rv*Grain(ii)%rscale(iopac),tot,C(D%nR/2,D%nTheta-1)%w0(ii)
	enddo
	enddo
	deallocate(spec)
	close(unit=90)

	do i=1,D%nR-1
		do j=1,D%nTheta-1
			if(computeTgas.or.viscous) then
				C(i,j)%Tgas=D%Tstar*sqrt(D%Rstar/D%R_av(i))
				if(C(i,j)%Tgas.gt.real(TMAX-1)*dT) C(i,j)%Tgas=real(TMAX-1)*dT
				C(i,j)%KappaGas=KappaGas(C(i,j)%gasdens*gas2dust,C(i,j)%Tgas)
			else
				C(i,j)%KappaGas=0d0
			endif
			C(i,j)%useFE=.false.
		enddo
	enddo

	if(multiwav) then
		allocate(column(ngrains,ngrains2))
		allocate(specemit(nlam))
	endif

	if(tdes_iter.and.Nphot.ne.0.and.(iter0.le.nBW.or.nBW.lt.0)) then
		if(RNDW.or.nphotdiffuse.gt.0) call InitRandomWalk()
		do iter=1,6	!15
		write(*,'("Iteration ",i2," of 6")') iter
		RmaxRefine=D%Rout
		if(struct_iter) then
			dosmooth=.false.
			call OpticallyThin(.false.)
			if (thgridrefine) call RegridTheta(D%nTheta/4)
			call TempAverage(1d0)
			call DiskStructure()
			dosmooth=.true.
		endif
		call OpticallyThin(.true.)
		if(iter.lt.7) then
			do i=0,D%nR-1
			do j=1,D%nTheta-1
				C(i,j)%gasfrac=1d0
				C(i,j)%dens=1d-60
			enddo
			enddo
		endif
		if(iter.lt.3) then		! 7) then
			call DestroyDustTDirect(C,Rdes,Rmin,1d10,.true.)
		else
			call DestroyDustT(C,Rdes,Rmin,1d10,.true.)
		endif
		call DestroyDustR()
		MassTot=0d0
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			MassTot=MassTot+C(i,j)%dens*C(i,j)%V
			if(number_invalid(C(i,j)%dens).ne.0) then
				print*,i,j,C(i,j)%dens,C(i,j)%w(1:ngrains),C(i,j)%TP(1:ngrains)
			endif
		enddo
		enddo
		print*,'MassTot:',MassTot/Msun
		write(*,'("Optically thin destruction R: ",f13.6," AU")') Rmin
		write(9,'("Optically thin destruction R: ",f13.6," AU")') Rmin
		if(Rdes.gt.D%R(1)) then
			write(*,'("No dust within:               ",f13.6," AU")') Rdes
			write(9,'("No dust within:               ",f13.6," AU")') Rdes
		endif
		if(iter.lt.4) then
			call RegridR(Rdes,RmaxRefine)
		else
			call RegridR(Rdes,RmaxRefine)
		endif
		enddo
		if(struct_iter) then
			dosmooth=.false.
			call OpticallyThin(.false.)
			call TempAverage(1d0)
			call DiskStructure()
			dosmooth=.true.
		endif
		call OpticallyThin(.true.)
		call DestroyDustT(C,Rdes,Rmin,1d10,.true.)
		write(*,'("Optically thin destruction R: ",f13.6," AU")') Rmin
		write(9,'("Optically thin destruction R: ",f13.6," AU")') Rmin
		if(Rdes.gt.D%R(1)) then
			write(*,'("No dust within:               ",f13.6," AU")') Rdes
			write(9,'("No dust within:               ",f13.6," AU")') Rdes
		endif
c		f_weight=f_weight_backup
	else if(denstype.ne.'FILE'.and.denstype.ne.'PREVIOUS'.and.gridrefine.and.startiter.eq.' ') then
c		f_weight_backup=f_weight
c		f_weight=1d0
		do iter=1,5
			RmaxRefine=D%Rout
			if(iter.lt.3) then
				call RegridR(D%R(2),RmaxRefine)
			else
				call RegridR(D%R(2),RmaxRefine)
			endif
			call DestroyDustR()
			if(denstype.eq.'PINTE') then
				MassTot=0d0
				D%nRfix=20
				if(allocated(D%Rfix)) deallocate(D%Rfix)
				allocate(D%Rfix(D%nRfix))
				do j=1,D%nRfix
					D%Rfix(j)=D%Rexp+(D%Rout-D%Rexp)*real(j-1)/real(D%nRfix)
				enddo
				do i=1,D%nR-1
				do j=1,D%nTheta-1
					r=D%R_av(i)*sin(D%theta_av(j))/AU
					z=D%R_av(i)*cos(D%theta_av(j))/AU
					hr=D%sh1AU*r**D%shpow
					f1=r**(-D%denspow)
					f2=exp(-(z/hr)**2)
					scale=D%Mtot/(AU**3*2d0*pi*D%sh1AU*sqrt(pi)*(D%Rout-D%Rin))
					if(r.lt.D%Rexp) then
						C(i,j)%dens=scale*f1*f2*r**(-D%shpow)
					else
						C(i,j)%dens=1d-60
					endif
					if(C(i,j)%dens.lt.1d-50) C(i,j)%dens=1d-50
					C(i,j)%mass=C(i,j)%dens*C(i,j)%V
					C(i,j)%dens0=C(i,j)%dens
					MassTot=MassTot+C(i,j)%mass
				enddo
				enddo
				do i=1,D%nR-1
				do j=1,D%nTheta-1
					C(i,j)%mass=D%Mtot*C(i,j)%mass/MassTot
					C(i,j)%dens=D%Mtot*C(i,j)%dens/MassTot
					C(i,j)%dens0=D%Mtot*C(i,j)%dens0/MassTot
				enddo
				enddo
				call DestroyDustR()
			endif
		enddo
c		f_weight=f_weight_backup
	endif
	if (use_topac) call Topac(0)
	tauincrease=1.25d0

	iter0=0
	if(startiter.ne.' '.or.denstype.eq.'PREVIOUS') then
		if(outputfits) then
			write(file,'(a,"denstemp",a,".fits.gz")') outdir(1:len_trim(outdir)),startiter(1:len_trim(startiter))
		else
			write(file,'(a,"denstemp",a,".dat")') outdir(1:len_trim(outdir)),startiter(1:len_trim(startiter))
		endif
		write(*,'("Reading starting iteration from: ",a)') file(1:len_trim(file))
		write(9,'("Reading starting iteration from: ",a)') file(1:len_trim(file))
		call readstruct(file,(/'DENS   ','TEMP   ','COMP   ','GASDENS','DENS0  '/),5,0,.false.)
		do i=1,D%nR-1
			do j=1,D%nTheta-1
				C(i,j)%TMC=C(i,j)%T
			enddo
		enddo

		if(computeTgas) then
			if(outputfits) then
				write(file,'(a,"denstempGas",a,".fits.gz")') outdir(1:len_trim(outdir)),startiter(1:len_trim(startiter))
			else
				write(file,'(a,"denstempGas",a,".dat")') outdir(1:len_trim(outdir)),startiter(1:len_trim(startiter))
			endif
			write(*,'("Reading gas temperatures from:   ",a)') file(1:len_trim(file))
			write(9,'("Reading gas temperatures from:   ",a)') file(1:len_trim(file))
			call readstruct(file,(/'GASDENS','GASTEMP'/),2,0,.false.)
			do i=1,D%nR-1
				do j=1,D%nTheta-1
					C(i,j)%KappaGas=KappaGas(C(i,j)%gasdens*gas2dust,C(i,j)%Tgas)
				enddo
			enddo
		endif

		D%nR=D%nR-1
		D%R(1)=D%Rin
		do i=2,D%nR
			D%R(i)=10d0**((2d0*log10(D%R_av(i-1)/AU)-log10(D%R(i-1))))
		enddo
		D%R(D%nR+1)=D%Rout
		D%R(0)=D%Rstar/AU
		call sort(D%R(1:D%nR),D%nR)
		D%nR=D%nR+1

		do j=0,D%nTheta-1
		do i=1,D%nR-1
			C(i,j)%T=0d0
			C(i,j)%Ni=0
			C(i,j)%EJv=0d0
			C(i,j)%iTD=0
			C(i,j)%xedge(1)=D%R(i)
			C(i,j)%xedge(2)=D%R(i+1)
			C(i,j)%xedge(3)=D%Theta(j)
			C(i,j)%xedge(4)=D%Theta(j+1)
			C(i,j)%xedge2(1)=D%R(i)**2
			C(i,j)%xedge2(2)=D%R(i+1)**2
			C(i,j)%xedge2(3)=D%Theta(j)**2
			C(i,j)%xedge2(4)=D%Theta(j+1)**2
		enddo
		C(0,j)%T=0d0
		C(0,j)%Ni=0
		C(0,j)%EJv=0d0
		C(0,j)%iTD=0
		C(0,j)%xedge(1)=D%R(0)
		C(0,j)%xedge(2)=D%R(1)
		C(0,j)%xedge(3)=D%Theta(j)
		C(0,j)%xedge(4)=D%Theta(j+1)
		C(0,j)%xedge2(1)=D%R(0)**2
		C(0,j)%xedge2(2)=D%R(1)**2
		C(0,j)%xedge2(3)=D%Theta(j)**2
		C(0,j)%xedge2(4)=D%Theta(j+1)**2
		enddo

		do i=1,D%nTheta-1
			C(0,i)%dens=1d-60
			C(0,i)%dens0=1d-60
			C(0,i)%V=(4d0*pi/3d0)*(D%R(1)**3-D%R(0)**3)*
     &			(D%Theta(i)-D%Theta(i+1))*AU**3
			C(0,i)%mass=C(0,i)%V*C(0,i)%dens
		enddo

		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
			C(i,j)%mass=C(i,j)%V*C(i,j)%dens
		enddo
		enddo

		shscale(0:D%nR)=1d0
		if(scalesh.eq.'FILE') then
			call regrid(shscalefile,D%R_av/AU,shscale,D%nR)
		else if(scalesh.eq.'VALUE') then
			shscale(0:D%nR)=shscalevalue
		else if(scalesh.eq.'POW') then
			shscale(0:D%nR-1)=1d0+(shscalevalue-1d0)*(D%R_av(0:D%nR-1)/AU)**(-shpow)
		else
			shscale(0:D%nR)=1d0
		endif
		if(startiter.ne.' ') read(startiter,*) iter0

	endif

	if(Nphot.le.0) then
		call reset2d()
		call ReadTemperatures()
	endif

	do ii=1,ngrains
		Grain(ii)%trace=trace(ii)
	enddo

	if((RNDW.or.nphotdiffuse.gt.0).and.Nphot.ne.0) call InitRandomWalk()

	if(reprocess.gt.0d0.and.Nphot.gt.0) then
	write(*,'("Scaling density structure to reprocessing: ",f10.3)') reprocess
	write(9,'("Scaling density structure to reprocessing: ",f10.3)') reprocess
	do iter=1,10
		tot=0d0
		tot2=0d0
		do j=1,D%nTheta-1
			radtau=0d0
			do i=1,D%nR-1
				do ii=1,ngrains
					do iopac=1,Grain(ii)%nopac
						radtau=radtau+(D%R(i+1)-D%R(i))*AU*C(i,j)%dens*C(i,j)%wopac(ii,iopac)*C(i,j)%w(ii)*Grain(ii)%Kpabsstar(iopac)
					enddo
				enddo
			enddo
			tot=tot+(D%Theta(j)-D%Theta(j+1))*(1d0-exp(-radtau))
			tot2=tot2+(D%Theta(j)-D%Theta(j+1))
		enddo
		print*,(tot/tot2)
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%dens=C(i,j)%dens*reprocess/(tot/tot2)
			C(i,j)%dens0=C(i,j)%dens0*reprocess/(tot/tot2)
			C(i,j)%gasdens=C(i,j)%gasdens*reprocess/(tot/tot2)
			C(i,j)%mass=C(i,j)%mass*reprocess/(tot/tot2)
		enddo
		enddo
		D%Mtot=D%Mtot*reprocess/(tot/tot2)
	enddo

	endif

c Set the spectrum and energy for the interstellar radiation field
	E_IRF=0d0
	if(use_IRF) then
		allocate(IRF(nlam))
		do i=1,nlam
			IRF(i)=pi*(D%R(D%nR)*AU)**2*(9.85357d-17*1.71*Planck(T_IRF,lam(i))*F_IRF
     &                                               +Planck(T_BG,lam(i)))
		enddo
		call integrate(IRF,E_IRF)
		if(T_ISMdust.gt.0d0) then
			allocate(spec(nlam))
			do i=1,nlam
				spec(i)=pi*(D%R(D%nR)*AU)**2*Planck(T_ISMdust,lam(i))
			enddo
			call integrate(spec,tot)
			call integrate(IRF,E_IRF)
			call ReadISRF(lam,spec,nlam)
			call integrate(spec,tot2)
			if(tot.gt.E_IRF) then
				spec=spec*(tot-E_IRF)/tot2
				IRF=IRF+spec
				call integrate(IRF,E_IRF)
			endif
			deallocate(spec)
		endif
		open(unit=45,file='ISRF.dat',RECL=500)
		do i=1,nlam
			write(45,*) lam(i),IRF(i)
		enddo
		close(unit=45)
	endif

	if(nplanets.ne.0) then
	allocate(Planets(nplanets))
	ia=1
	i=0
600	call getarg(2+ia,line)
	if(line.eq.' ') goto 601
	if(line(1:2).eq.'-p') then
		i=i+1
		call getarg(3+ia,Planets(i)%file)
	endif
	ia=ia+1
	goto 600
601	continue
	endif

	allocate(ncoolingtime(0:D%nR))
	allocate(coolingtime(0:D%nR))

	if(Nphot.gt.0) then
		do i=0,D%nR-1
			do j=1,D%nTheta-1
				call CheckMinimumDensity(i,j)
			enddo
		enddo
	endif

	write(surfdensfile,'(a,"surfacedens.dat")') outdir(1:len_trim(outdir))
	open(unit=90,file=surfdensfile,RECL=100)
	do i=1,D%nR-1
	Vtot=0d0
	MassTot=0d0
	tot=0d0
	MassTot0=0d0
	do j=1,D%nTheta
		MassTot=MassTot+C(i,j)%mass
		tot=tot+C(i,j)%dens0*C(i,j)%V
		MassTot0=MassTot0+C(i,j)%gasdens*C(i,j)%V
	enddo
	write(90,*) D%R_av(i)/AU,MassTot/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2),tot/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)
     &						,MassTot0*gas2dust/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)
	enddo
	close(unit=90)

	call KappaAverage()

	arraysallocated=.true.
	
	NsigDiskstructure=1

	deallocate(computepart_T)
	deallocate(computepart_Tfile)

	return
	end

c-----------------------------------------------------------------------
c This subroutine deallocates all arrays set in the subroutine 
c initialize.
c-----------------------------------------------------------------------

	subroutine FreeAllMem()
	use Parameters
	use Diskstruct	
	IMPLICIT NONE

	do i=1,ngrains
		deallocate(Grain(i)%Kabs)
		deallocate(Grain(i)%Ksca)
		deallocate(Grain(i)%Kext)
		deallocate(Grain(i)%g)
		deallocate(Grain(i)%F)
		if (allocated(Grain(i)%Topac)) deallocate(Grain(i)%Topac)
		if (allocated(Grain(i)%cryst)) deallocate(Grain(i)%cryst)
		deallocate(Grain(i)%KabsL)
		deallocate(Grain(i)%KextL)
	enddo
	do i=0,D%nR
	do j=0,D%nTheta
		deallocate(C(i,j)%w)
		if(storescatt) deallocate(C(i,j)%scattfield)
		if(.not.tcontact.or.tdes_iter.or.use_qhp) then
			deallocate(C(i,j)%TP)
			deallocate(C(i,j)%EJvP)
			deallocate(C(i,j)%EabsP)
		endif
		if(allocated(C(i,j)%QHP)) deallocate(C(i,j)%QHP)
		if(allocated(C(i,j)%LRF)) deallocate(C(i,j)%LRF)
		if(allocated(C(i,j)%nLRF)) deallocate(C(i,j)%nLRF)
		if(allocated(C(i,j)%tdistr)) deallocate(C(i,j)%tdistr)
		if(allocated(C(i,j)%EJvQHP)) deallocate(C(i,j)%EJvQHP)
		if(allocated(C(i,j)%EabsQHP)) deallocate(C(i,j)%EabsQHP)
		if(allocated(C(i,j)%EviscDirect)) deallocate(C(i,j)%EviscDirect)
		if(allocated(C(i,j)%FE)) deallocate(C(i,j)%FE)
	enddo
	enddo
	if(allocated(IRF)) deallocate(IRF)
	if(storescatt) then
		deallocate(scattcomputed)
		deallocate(nscattcomputed)
	endif
	deallocate(Grain)

	deallocate(C)
	deallocate(D%theta_av)
	deallocate(D%Theta)
	deallocate(D%thet)
	deallocate(D%R_av)
	deallocate(D%R)
	if(allocated(D%gap)) then
		deallocate(D%gap)
		deallocate(D%gap1)
		deallocate(D%gap2)
		deallocate(D%gapshape)
		deallocate(D%gaproundtype)
		deallocate(D%gaproundpow)
	endif

	if (allocated(D%MdotR)) deallocate(D%MdotR)
	if (allocated(D%MPdead)) deallocate(D%MPdead)

	deallocate(lam)
	deallocate(BB)
	deallocate(dnu)
	deallocate(nu)
	deallocate(D%Fstar)
	deallocate(lamHR)
	deallocate(FstarHR)

	arraysallocated=.false.

	if(allocated(Temp)) deallocate(Temp)
	if(allocated(fBW)) deallocate(fBW)

	if(allocated(D%Rfix)) deallocate(D%Rfix)
	if(allocated(muRad)) deallocate(muRad)

	if(multiwav) then
		deallocate(column)
		deallocate(specemit)
	endif
	if(allocated(D%SinTheta)) deallocate(D%SinTheta)
	if(allocated(shscale)) deallocate(shscale)
	if(allocated(Rthindes)) deallocate(Rthindes)
	if(allocated(Rbwdes)) deallocate(Rbwdes)
	if(allocated(iRfirst)) deallocate(iRfirst)
	if(allocated(Planets)) deallocate(Planets)

	return
	end


c-----------------------------------------------------------------------
c This subroutine resets the temperature and photon statistics in the
c disk.
c-----------------------------------------------------------------------

	subroutine reset2d()
	use Parameters
	IMPLICIT NONE
	integer i,j,ii
		
	do j=1,D%nTheta-1
	do i=0,D%nR-1
		C(i,j)%T=0d0
		C(i,j)%TMC=0d0
		C(i,j)%dT=0d0
		C(i,j)%Ni=0
		C(i,j)%EJv=0d0
		C(i,j)%EJv2=0d0
		C(i,j)%Egas=0d0
		C(i,j)%lastphotnr=0
		C(i,j)%randomwalk=.false.
		if(storescatt) C(i,j)%scattfield(1,0:nlam,1:2)=0d0
		if(.not.tcontact.or.tdes_iter) then
			do ii=1,ngrains
				C(i,j)%TP(ii)=0d0
				C(i,j)%EJvP(ii)=0d0
			enddo
		endif
		if(use_qhp) then
			if(qhp_solver.eq.2) C(i,j)%Tqhp(1:nqhp)=0d0
			C(i,j)%EJvQHP(1:nqhp)=0d0
			C(i,j)%EabsQHP(1:nqhp)=0d0
		endif
		if(computeLRF) then
			C(i,j)%LRF(1:nlam)=0d0
			C(i,j)%nLRF(1:nlam)=0
		endif
		C(i,j)%FradR=0d0
		C(i,j)%FradZ=0d0
	enddo
	enddo

	return
	end

c-----------------------------------------------------------------------
c This subroutine reads a 'particle' file. This is a specially setup
c file containing the opacities and scattering phase functions of 
c and ensemble of particles at each wavelength.
c The interpolation is done on the set wavelength grid.
c-----------------------------------------------------------------------
	subroutine readparticle(input,p,Rayleigh,asym,asym2,wasym2,Pmax,iopac)
	use Parameters
	IMPLICIT NONE
	type(particle) p,p0,p1
	integer i,j,ia,iopac
	character*500 input
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
	
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
	   write(*,'("Particle file not found: ",a)') input(1:len_trim(input))
	   write(9,'("Particle file not found: ",a)') input(1:len_trim(input))
	   write(*,'("--------------------------------------------------------")')
	   write(9,'("--------------------------------------------------------")')
	   stop
	endif
	write(*,'("Reading particle file: ",a)') input(1:len_trim(input))
	write(9,'("Reading particle file: ",a)') input(1:len_trim(input))
	open(unit=20,file=input,RECL=6000)
	i=1
	read(20,*,end=102) l0,p0%Kext(1,1),p0%Kabs(1,1),p0%Ksca(1,1)
	read(20,*,end=102) (p0%F(1,1)%F11(j),j=1,180)
	read(20,*,end=102) (p0%F(1,1)%F12(j),j=1,180)
	read(20,*,end=102) (p0%F(1,1)%F22(j),j=1,180)
	read(20,*,end=102) (p0%F(1,1)%F33(j),j=1,180)
	read(20,*,end=102) (p0%F(1,1)%F34(j),j=1,180)
	read(20,*,end=102) (p0%F(1,1)%F44(j),j=1,180)
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
100	read(20,*,end=102) l1,p1%Kext(1,1),p1%Kabs(1,1),p1%Ksca(1,1)
	read(20,*,end=102) (p1%F(1,1)%F11(j),j=1,180)
	read(20,*,end=102) (p1%F(1,1)%F12(j),j=1,180)
	read(20,*,end=102) (p1%F(1,1)%F22(j),j=1,180)
	read(20,*,end=102) (p1%F(1,1)%F33(j),j=1,180)
	read(20,*,end=102) (p1%F(1,1)%F34(j),j=1,180)
	read(20,*,end=102) (p1%F(1,1)%F44(j),j=1,180)
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
	close(unit=20)


	if(nspike.gt.0.and.nspike.lt.180.and.scat_how.gt.1) then
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

	if(nspike.gt.0.and.nspike.lt.180.and.scat_how.gt.1) then
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


	return
	end

c-----------------------------------------------------------------------
c This subroutine reads an 'opcaity' file. This is a file containing 
c the opacities of scattering and absorption. The scattering matrix is
c assumed Rayleigh. The file is assumed to be as follows:
c
c 1st column: wavelenght
c 2nd column: Kext (cm^2/g)
c 3th column: Kabs (cm^2/g)
c 4th column: Ksca (cm^2/g)
c
c The interpolation is done on the set wavelength grid.
c-----------------------------------------------------------------------
	subroutine readopacity(input,p)
	use Parameters
	IMPLICIT NONE
	type(particle) p,p0,p1
	integer i,j,ia
	character*500 input
	logical truefalse
	real*8 l0,l1,tot,tot2,theta

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
	
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,'("Opacity file not found")')
		write(9,'("Opacity file not found")')
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		stop
	endif
	write(*,'("Reading particle file: ",a)') input(1:len_trim(input))
	write(9,'("Reading particle file: ",a)') input(1:len_trim(input))
	open(unit=20,file=input,RECL=6000)
	i=1
	read(20,*,end=102) l0,p0%Kext(1,1),p0%Kabs(1,1),p0%Ksca(1,1)
	do ia=1,180
		theta=(real(ia)-0.5d0)*pi/180d0
		p0%F(1,1)%F11(ia)=(1d0+cos(theta)**2)/2d0
		p0%F(1,1)%F12(ia)=-(1d0-cos(theta)**2)/2d0
		p0%F(1,1)%F22(ia)=(1d0+cos(theta)**2)/2d0
		p0%F(1,1)%F33(ia)=cos(theta)
		p0%F(1,1)%F34(ia)=0d0
		p0%F(1,1)%F44(ia)=cos(theta)
	enddo
103	if(l0.ge.lam(i)) then
		p%Kext(1,i)=p0%Kext(1,1)
		p%Ksca(1,i)=p0%Ksca(1,1)
		p%Kabs(1,i)=p0%Kabs(1,1)
		p%F(1,i)=p0%F(1,1)
		call tellertje(i,nlam)
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) l1,p1%Kext(1,1),p1%Kabs(1,1),p1%Ksca(1,1)
	do ia=1,180
		theta=(real(ia)-0.5d0)*pi/180d0
		p1%F(1,1)%F11(ia)=(1d0+cos(theta)**2)/2d0
		p1%F(1,1)%F12(ia)=-(1d0-cos(theta)**2)/2d0
		p1%F(1,1)%F22(ia)=(1d0+cos(theta)**2)/2d0
		p1%F(1,1)%F33(ia)=cos(theta)
		p1%F(1,1)%F34(ia)=0d0
		p1%F(1,1)%F44(ia)=cos(theta)
	enddo
101	if(lam(i).le.l1.and.lam(i).ge.l0) then
		p%Kext(1,i)=p1%Kext(1,1)+(lam(i)-l1)*(p0%Kext(1,1)-p1%Kext(1,1))/(l0-l1)
		p%Ksca(1,i)=p1%Ksca(1,1)+(lam(i)-l1)*(p0%Ksca(1,1)-p1%Ksca(1,1))/(l0-l1)
		p%Kabs(1,i)=p1%Kabs(1,1)+(lam(i)-l1)*(p0%Kabs(1,1)-p1%Kabs(1,1))/(l0-l1)
		p%F(1,i)%F11(1:180)=p1%F(1,1)%F11(1:180)+(lam(i)-l1)*(p0%F(1,1)%F11(1:180)-p1%F(1,1)%F11(1:180))/(l0-l1)
		p%F(1,i)%F12(1:180)=p1%F(1,1)%F12(1:180)+(lam(i)-l1)*(p0%F(1,1)%F12(1:180)-p1%F(1,1)%F12(1:180))/(l0-l1)
		p%F(1,i)%F22(1:180)=p1%F(1,1)%F22(1:180)+(lam(i)-l1)*(p0%F(1,1)%F22(1:180)-p1%F(1,1)%F22(1:180))/(l0-l1)
		p%F(1,i)%F33(1:180)=p1%F(1,1)%F33(1:180)+(lam(i)-l1)*(p0%F(1,1)%F33(1:180)-p1%F(1,1)%F33(1:180))/(l0-l1)
		p%F(1,i)%F34(1:180)=p1%F(1,1)%F34(1:180)+(lam(i)-l1)*(p0%F(1,1)%F34(1:180)-p1%F(1,1)%F34(1:180))/(l0-l1)
		p%F(1,i)%F44(1:180)=p1%F(1,1)%F44(1:180)+(lam(i)-l1)*(p0%F(1,1)%F44(1:180)-p1%F(1,1)%F44(1:180))/(l0-l1)
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
		p%Ksca(1,j)=p%Ksca(1,i-1)*(lam(i-1)/lam(j))**4
		p%Kabs(1,j)=p%Kabs(1,i-1)*(lam(i-1)/lam(j))**2
		p%Kext(1,j)=p%Kabs(1,j)+p%Ksca(1,j)
		p%F(1,j)=p%F(1,i-1)
	enddo
	close(unit=20)
	do j=1,nlam
	tot=0d0
	tot2=0d0
	do i=1,180
		tot=tot+p%F(1,j)%F11(i)*sin(pi*(real(i)-0.5)/180d0)
		tot2=tot2+sin(pi*(real(i)-0.5)/180d0)
	enddo
	do i=1,180
		p%F(1,j)%F11(i)=tot2*p%F(1,j)%F11(i)/tot
		p%F(1,j)%F12(i)=tot2*p%F(1,j)%F12(i)/tot
		p%F(1,j)%F22(i)=tot2*p%F(1,j)%F22(i)/tot
		p%F(1,j)%F33(i)=tot2*p%F(1,j)%F33(i)/tot
		p%F(1,j)%F34(i)=tot2*p%F(1,j)%F34(i)/tot
		p%F(1,j)%F44(i)=tot2*p%F(1,j)%F44(i)/tot
	enddo
	enddo

	deallocate(p0%Kabs)
	deallocate(p0%Ksca)
	deallocate(p0%Kext)
	deallocate(p0%F)
	deallocate(p1%Kabs)
	deallocate(p1%Ksca)
	deallocate(p1%Kext)
	deallocate(p1%F)


	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	recursive subroutine MakeMatrixCoupled(abunA,coupledabun,i,k,ngrains)
	IMPLICIT NONE
	integer ngrains,i,j,k,coupledabun(*)
	real*8 abunA(ngrains,ngrains)

	do j=1,ngrains
		if(coupledabun(j).eq.i) then
			abunA(k,j)=1d0
			call MakeMatrixCoupled(abunA,coupledabun,j,k,ngrains)
		endif
	enddo

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine readcomposition(input)
	use Parameters
	IMPLICIT NONE
	integer i,j,jj
	real*8 x0,y0(ngrains),x1,y1(ngrains),tot
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,RECL=6000)
	i=1
	read(20,*,end=102) x0,y0(1:ngrains)
103	if(x0.ge.(D%R_av(i)/AU)) then
		do jj=1,D%nTheta-1
			C(i,jj)%w(1:ngrains)=y0(1:ngrains)
			tot=sum(C(i,jj)%w(1:ngrains))
			if(tot.eq.0d0) then
				C(i,jj)%w(1:ngrains)=1d0/real(ngrains)
				C(i,jj)%w0(1:ngrains)=1d0/real(ngrains)
			else
				C(i,jj)%w(1:ngrains)=C(i,jj)%w(1:ngrains)/tot
				C(i,jj)%w0(1:ngrains)=C(i,jj)%w(1:ngrains)
			endif
		enddo
		i=i+1
		if(i.gt.D%nR) goto 102
		goto 103
	endif
100	read(20,*,end=102) x1,y1(1:ngrains)
101	if((D%R_av(i)/AU).le.x1.and.(D%R_av(i)/AU).gt.x0) then
		do jj=1,D%nTheta-1
			C(i,jj)%w(1:ngrains)=y1(1:ngrains)+((D%R_av(i)/AU)-x1)*(y0(1:ngrains)-y1(1:ngrains))/(x0-x1)
			tot=sum(C(i,jj)%w(1:ngrains))
			if(tot.eq.0d0) then
				C(i,jj)%w(1:ngrains)=1d0/real(ngrains)
				C(i,jj)%w0(1:ngrains)=1d0/real(ngrains)
			else
				C(i,jj)%w(1:ngrains)=C(i,jj)%w(1:ngrains)/tot
				C(i,jj)%w0(1:ngrains)=C(i,jj)%w(1:ngrains)
			endif
		enddo
		i=i+1
		if(i.gt.D%nR) goto 102
		goto 101
	endif
	x0=x1
	y0(1:ngrains)=y1(1:ngrains)
	goto 100
102	continue
	do j=i,D%nR
		do jj=1,D%nTheta-1
			C(j,jj)%w(1:ngrains)=C(i-1,jj)%w(1:ngrains)*D%R_av(i-1)/D%R_av(j)
			tot=sum(C(i,jj)%w(1:ngrains))
			if(tot.eq.0d0) then
				C(i,jj)%w(1:ngrains)=1d0/real(ngrains)
				C(i,jj)%w0(1:ngrains)=1d0/real(ngrains)
			else
				C(i,jj)%w(1:ngrains)=C(i,jj)%w(1:ngrains)/tot
				C(i,jj)%w0(1:ngrains)=C(i,jj)%w(1:ngrains)
			endif
		enddo
	enddo
	close(unit=20)
	return
	end


	subroutine checkaggregates(file,nopac)
	use Parameters
	IMPLICIT NONE
	character*500 file,name
	real*8 fract
	integer nopac
	
	open(unit=50,file=file,RECL=6000)
	
	nopac=0
1	read(50,*,end=2) fract
	read(50,*) name
	nopac=nopac+1
	goto 1
2	close(unit=50)
	
	return
	end
	
	

c-----------------------------------------------------------------------
c This subroutine reads an 'qhp' file. This is a file containing 
c the opacities of scattering and absorption and the single photon
c emission for each wavelength. The phase function is
c assumed isotropic. The file is assumed to be as follows:
c
c 1st row contains Na Ma Td: The # of atoms, the weight of one atom (in proton mass)
c								and the Td (450 for C, 175 for silicates)
c
c 1st column: wavelenght
c 2nd column: Kext (cm^2/g)
c 3th column: Kabs (cm^2/g)
c 4th column: Ksca (cm^2/g)
c
c The interpolation is done on the set wavelength grid.
c-----------------------------------------------------------------------
	subroutine readQHP(input,p)
	use Parameters
	IMPLICIT NONE
	type(particle) p,p0,p1
	integer i,j,ia
	character*500 input,line
	logical truefalse
	real*8 l0,l1,tot,tot2,theta,Mc
	parameter(Mc=12d0*1.66d-24) !mass of a carbon atom in gram

	p%qhp=.true.
	p%gascoupled=.true.

	allocate(p0%Kabs(1,nlam))
	allocate(p0%Ksca(1,nlam))
	allocate(p0%Kext(1,nlam))
	allocate(p0%F(1,nlam))
	allocate(p1%Kabs(1,nlam))
	allocate(p1%Ksca(1,nlam))
	allocate(p1%Kext(1,nlam))
	allocate(p1%F(1,nlam))
	
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,'("Opacity file not found")')
		write(9,'("Opacity file not found")')
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		stop
	endif
	write(*,'("Reading particle file: ",a)') input(1:len_trim(input))
	write(9,'("Reading particle file: ",a)') input(1:len_trim(input))
	open(unit=20,file=input,RECL=6000)
	read(20,'(a500)',end=1,err=1) line
	read(line,*,end=1,err=1) p%Nc,p%Mc,p%Td_qhp
	goto 2
1	continue
	read(line,*) p%Nc
	p%Mc=12d0
	p%Td_qhp=450d0
2	continue
	print*,p%Nc,p%Mc,p%Td_qhp

	p%rv=(p%Nc/468d0)**(1d0/3d0)*1d-7
	p%rho(1:p%nopac)=p%Nc*(p%Mc*Mc/12d0)/(4d0*pi*p%rv**3/3d0)
	print*,p%Nc,p%rv*1d4,p%rho(1:p%nopac)

	i=1
	read(20,*,end=102) l0,p0%Kext(1,1),p0%Kabs(1,1),p0%Ksca(1,1)
	do ia=1,180
		theta=(real(ia)-0.5d0)*pi/180d0
		p0%F(1,1)%F11(ia)=(1d0+cos(theta)**2)/2d0
		p0%F(1,1)%F12(ia)=-(1d0-cos(theta)**2)/2d0
		p0%F(1,1)%F22(ia)=(1d0+cos(theta)**2)/2d0
		p0%F(1,1)%F33(ia)=cos(theta)
		p0%F(1,1)%F34(ia)=0d0
		p0%F(1,1)%F44(ia)=cos(theta)
	enddo
103	if(l0.ge.lam(i)) then
		p%Kext(1,i)=p0%Kext(1,1)
		p%Ksca(1,i)=p0%Ksca(1,1)
		p%Kabs(1,i)=p0%Kabs(1,1)
		p%F(1,i)=p0%F(1,1)
		call tellertje(i,nlam)
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) l1,p1%Kext(1,1),p1%Kabs(1,1),p1%Ksca(1,1)
	do ia=1,180
		theta=(real(ia)-0.5d0)*pi/180d0
		p1%F(1,1)%F11(ia)=(1d0+cos(theta)**2)/2d0
		p1%F(1,1)%F12(ia)=-(1d0-cos(theta)**2)/2d0
		p1%F(1,1)%F22(ia)=(1d0+cos(theta)**2)/2d0
		p1%F(1,1)%F33(ia)=cos(theta)
		p1%F(1,1)%F34(ia)=0d0
		p1%F(1,1)%F44(ia)=cos(theta)
	enddo
101	if(lam(i).le.l1.and.lam(i).ge.l0) then
		p%Kext(1,i)=p1%Kext(1,1)+(lam(i)-l1)*(p0%Kext(1,1)-p1%Kext(1,1))/(l0-l1)
		p%Ksca(1,i)=p1%Ksca(1,1)+(lam(i)-l1)*(p0%Ksca(1,1)-p1%Ksca(1,1))/(l0-l1)
		p%Kabs(1,i)=p1%Kabs(1,1)+(lam(i)-l1)*(p0%Kabs(1,1)-p1%Kabs(1,1))/(l0-l1)
		p%F(1,i)%F11(1:180)=p1%F(1,1)%F11(1:180)+(lam(i)-l1)*(p0%F(1,1)%F11(1:180)-p1%F(1,1)%F11(1:180))/(l0-l1)
		p%F(1,i)%F12(1:180)=p1%F(1,1)%F12(1:180)+(lam(i)-l1)*(p0%F(1,1)%F12(1:180)-p1%F(1,1)%F12(1:180))/(l0-l1)
		p%F(1,i)%F22(1:180)=p1%F(1,1)%F22(1:180)+(lam(i)-l1)*(p0%F(1,1)%F22(1:180)-p1%F(1,1)%F22(1:180))/(l0-l1)
		p%F(1,i)%F33(1:180)=p1%F(1,1)%F33(1:180)+(lam(i)-l1)*(p0%F(1,1)%F33(1:180)-p1%F(1,1)%F33(1:180))/(l0-l1)
		p%F(1,i)%F34(1:180)=p1%F(1,1)%F34(1:180)+(lam(i)-l1)*(p0%F(1,1)%F34(1:180)-p1%F(1,1)%F34(1:180))/(l0-l1)
		p%F(1,i)%F44(1:180)=p1%F(1,1)%F44(1:180)+(lam(i)-l1)*(p0%F(1,1)%F44(1:180)-p1%F(1,1)%F44(1:180))/(l0-l1)
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
		p%Ksca(1,j)=p%Ksca(1,i-1)*(lam(i-1)/lam(j))**4
		p%Kabs(1,j)=p%Kabs(1,i-1)*(lam(i-1)/lam(j))**2
		p%Kext(1,j)=p%Kabs(1,j)+p%Ksca(1,j)
		p%F(1,j)=p%F(1,i-1)
	enddo
	close(unit=20)
	do j=1,nlam
	tot=0d0
	tot2=0d0
	do i=1,180
		tot=tot+p%F(1,j)%F11(i)*sin(pi*(real(i)-0.5)/180d0)
		tot2=tot2+sin(pi*(real(i)-0.5)/180d0)
	enddo
	do i=1,180
		p%F(1,j)%F11(i)=tot2*p%F(1,j)%F11(i)/tot
		p%F(1,j)%F12(i)=tot2*p%F(1,j)%F12(i)/tot
		p%F(1,j)%F22(i)=tot2*p%F(1,j)%F22(i)/tot
		p%F(1,j)%F33(i)=tot2*p%F(1,j)%F33(i)/tot
		p%F(1,j)%F34(i)=tot2*p%F(1,j)%F34(i)/tot
		p%F(1,j)%F44(i)=tot2*p%F(1,j)%F44(i)/tot
	enddo
	enddo

	deallocate(p0%Kabs)
	deallocate(p0%Ksca)
	deallocate(p0%Kext)
	deallocate(p0%F)
	deallocate(p1%Kabs)
	deallocate(p1%Ksca)
	deallocate(p1%Kext)
	deallocate(p1%F)


	
	return
	end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine MakeGaps()
	use Parameters
	IMPLICIT NONE
	integer i,j,k
	real*8 MassTot,MassTot0,r,scale,RoundOff

	if(D%ngap.eq.0) return

	MassTot0=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		MassTot0=MassTot0+C(i,j)%mass
	enddo
	enddo
	
	do i=1,D%nR-1
	   do j=1,D%nTheta-1
	      C(i,j)%densscale=1d0
	      do k=1,D%ngap
		 r=D%R_av(i)*sin(D%theta_av(j))**D%gapshape(k)/AU
		 if(r.gt.D%gap1(k).and.r.lt.D%gap2(k)) then
		    if (D%gaproundtype(k).ne.' ') then 
		       scale=RoundOff(r,D%gap2(k),D%gaproundtype(k),
     1                                D%gaproundpow(k),D%gap(k))
		    else
		       scale=D%gap(k)
		    endif
		    C(i,j)%dens=C(i,j)%dens*scale
		    C(i,j)%dens0=C(i,j)%dens0*scale
		    C(i,j)%mass=C(i,j)%dens*C(i,j)%V
		    C(i,j)%densscale=C(i,j)%densscale*scale
		 endif
	      enddo
	   enddo
	enddo
	MassTot=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		MassTot=MassTot+C(i,j)%mass
	enddo
	enddo
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		C(i,j)%mass=C(i,j)%mass*MassTot0/MassTot
		C(i,j)%dens=C(i,j)%dens*MassTot0/MassTot
		C(i,j)%dens0=C(i,j)%dens0*MassTot0/MassTot
		call CheckMinimumDensity(i,j)
	enddo
	enddo

	return
	end



	subroutine UnMakeGaps()
	use Parameters
	IMPLICIT NONE
	integer i,j,k
	real*8 MassTot,MassTot0

	if(D%ngap.eq.0) return

	MassTot0=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		MassTot0=MassTot0+C(i,j)%mass
	enddo
	enddo
	
	do i=1,D%nR-1
		do j=1,D%nTheta-1
			if(C(i,j)%densscale.ne.0d0) then
				C(i,j)%dens=C(i,j)%dens/C(i,j)%densscale
				C(i,j)%dens0=C(i,j)%dens0/C(i,j)%densscale
				C(i,j)%mass=C(i,j)%dens*C(i,j)%V
				C(i,j)%densscale=1d0
			endif
		enddo
	enddo

	MassTot=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		MassTot=MassTot+C(i,j)%mass
	enddo
	enddo
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		C(i,j)%mass=C(i,j)%mass*MassTot0/MassTot
		C(i,j)%dens=C(i,j)%dens*MassTot0/MassTot
		C(i,j)%dens0=C(i,j)%dens0*MassTot0/MassTot
		call CheckMinimumDensity(i,j)
	enddo
	enddo

	return
	end

c-----------------------------------------------------------------------
c
c   This subroutine calculates the shape for a rounded of gap edge/rim
c   It returns only a factor which is the decrease in density, so it
c   can be used with dust destruction [DestroyDustR()] or for the total
c   surface density [MakeGaps()] 
c
c   + roundtype='softedge': a soft edge inside of rmax, based on angular momentum 
c                           conservation. Described in Woitke++ 2009.
c                           (only a 10% effect, needs grid refinement and benchmark)
c
c   + roundtype='powerlaw': a powerlaw surface density between up to rmax,
c                           proportional to r^-roundpow.
c
c   + roundtype='hydro':    a gaussian like shape that fits hydro simulations
c                           of embedded planets well (Lubow & D'Angelo 2006)
c
c
c   OUPUT:  scaling factor at radius r
c
c-----------------------------------------------------------------------

	function RoundOff(r,rmax,type,pow,scalemin)
	use Parameters
	IMPLICIT NONE
	real*8 RoundOff

	character*20 type
	integer i,j,k
	real*8 scale,r,rmax,pow,scalemin
	doubleprecision, PARAMETER :: GG=6.6720000d-08
	doubleprecision, PARAMETER :: amu=1.66053886d-24
	
	scale=scalemin

	!  do not apply a scale factor if outside rmax
	if(r.ge.rmax) then
	   scale= 1d0

        !  use a soft edge
	else if(type.eq.'softedge') then
	   scale=       (2.3 * amu)/(kb*C(i,D%nTheta-1)%T)* (GG * D%Mstar)/D%R_av(i)
	   scale=scale* (1d0 - 0.5d0*(rmax/r - r/rmax))
	   scale=exp(scale) * ((r/rmax)**(D%denspow)) 
	
        !  use a shape from Hydrodynamic simulations
	else if (type.eq.'hydro') then
	   scale=exp( -(((1-r/rmax)/pow)**3d0))

	!  use a powerlaw
	else if (type.eq.'powerlaw') then
	   scale=(rmax/r)**(-pow)

	   ! scale from Sigma=r^-p (not robust for p =!= 1)
	   scale=scale* ((r/rmax)**(D%denspow)) 
	endif

	!  return gap depth
	RoundOff=max(scale,scalemin)

	return
	end

c-----------------------------------------------------------------------


	

	subroutine ReadPlanet(P,lam0)
	use Parameters
	IMPLICIT NONE
	type(ExoPlanet) P
	character*500 key,value,line,Ffile
	integer i,ia
	real*8 theta,ran2,lam0,tot,tot2,lam1,lam2,A1,A2
	real*8 F11_1(180),F11_2(180),F12_1(180),F12_2(180)
	logical phiset,Tset
	
	P%R=1d0
	P%d=5d0
	P%phi=360d0*ran2(idum)
	P%theta=pi/2d0
	p%A=1d0
	P%P=0.3d0

	P%name=P%file(1:index(P%file,'.')-1)

	phiset=.false.
	Tset=.false.

	Ffile=' '

	open(unit=50,file=P%file)

	write(*,'("Reading planet file: ",a)') P%file(1:len_trim(P%file))
	write(9,'("Reading planet file: ",a)') P%file(1:len_trim(P%file))

1	call ignorestar(50)
	read(50,'(a500)',end=2,err=2) line

	key=line(1:index(line,'=')-1)
	value=line(index(line,'=')+1:len_trim(line))

	do i=1,len_trim(key)
		if(iachar(key(i:i)).ge.65.and.iachar(key(i:i)).le.90) then
			key(i:i)=achar(iachar(key(i:i))+32)
		endif
	enddo
	if(value(1:1).eq.'"'.or.value(1:1).eq."'") then
		value=value(2:len_trim(value)-1)
	endif

	if(key.eq.'r') read(value,*) P%R ! should be in km
	if(key.eq.'phi') then
		read(value,*) P%phi
		phiset=.true.
	endif
	if(key.eq.'t') then
		read(value,*) P%T
		Tset=.true.
	endif
	if(key.eq.'albedo') read(value,*) P%A
	if(key.eq.'pmax') read(value,*) P%P
	if(key.eq.'distance') read(value,*) P%d ! should be in AU
	if(key.eq.'name') P%name=value
	if(key.eq.'file') Ffile=value
	if(key.eq.'rantheta') P%theta=acos(1d0-2d0*ran2(idum))

	goto 1

2	continue
	close(unit=50)
	
	if(Ffile.eq.' ') then
		do i=1,180
			theta=(real(i)-0.5d0)*pi/180d0
			P%F%F11(i)=(8d0/3d0)*(sin(theta)-theta*cos(theta))/pi
			P%F%F12(i)=-P%P*P%F%F11(i)*(1d0-cos(theta)**2)/(1d0+cos(theta)**2)
			P%F%F22(i)=P%F%F11(i)
			P%F%F33(i)=cos(theta)
			P%F%F34(i)=0d0
			P%F%F44(i)=cos(theta)
		enddo
	else
		write(*,'("  reading file: ",a)') Ffile(1:len_trim(Ffile))
		write(9,'("  reading file: ",a)') Ffile(1:len_trim(Ffile))
		open(unit=50,file=Ffile,RECL=6000)
		read(50,*) lam1,A1
		read(50,*) F11_1(1:180)
		read(50,*) F12_1(1:180)
3		read(50,*,end=4) lam2,A2
		read(50,*) F11_2(1:180)
		read(50,*) F12_2(1:180)
		if(lam2.ge.lam0) goto 4
		lam1=lam2
		F11_1=F11_2
		F12_1=F12_2
		A1=A2
		goto 3
4		continue
		close(unit=50)
		P%A=A1+(A2-A1)*(lam0-lam1)/(lam2-lam1)
		P%F%F11=F11_1+(F11_2-F11_1)*(lam0-lam1)/(lam2-lam1)
		P%F%F12=F12_1+(F12_2-F12_1)*(lam0-lam1)/(lam2-lam1)
		do i=1,180
			theta=(real(i)-0.5d0)*pi/180d0
			P%F%F22(i)=P%F%F11(i)
			P%F%F33(i)=cos(theta)
			P%F%F34(i)=0d0
			P%F%F44(i)=cos(theta)
		enddo
	endif
	
	write(*,'("  Albedo: ",f15.3)') P%A
	write(*,'("  Polarization at 90 degrees: ",f7.2,"%")') -100d0*P%F%F12(90)/P%F%F11(90)
	write(*,'("  Phi:     ",f15.1," degrees")') P%phi
	write(*,'("  Theta:   ",f15.1," degrees")') P%theta*180d0/pi

	tot=0d0
	tot2=0d0
	do ia=1,180
		theta=pi*(real(ia)-0.5)/180d0
		tot=tot+p%F%F11(ia)*sin(theta)
		tot2=tot2+sin(theta)
	enddo
	p%F%F11=p%F%F11*tot2/tot
	p%F%F12=p%F%F12*tot2/tot
	p%F%F22=p%F%F22*tot2/tot
	p%F%F33=p%F%F33*tot2/tot
	p%F%F34=p%F%F34*tot2/tot
	p%F%F44=p%F%F44*tot2/tot

	P%R=P%R*1d5

	if(.not.phiset) then
		P%phi=360d0*lifetime/sqrt(P%d**3*Msun/D%Mstar)
	endif
	
5	if(P%phi.gt.360d0) then
		P%phi=P%phi-360d0
		goto 5
	endif
	if(P%phi.lt.0d0) then
		P%phi=P%phi+360d0
		goto 5
	endif

	if(.not.Tset) then
		P%T=D%Tstar*(P%d*AU/D%Rstar)**(-0.5)
	endif
	
	return
	end
	
	

	
	real*8 function HG(g,theta)
	IMPLICIT NONE
	real*8 g,theta,t,tot,f
	integer i,n

c	HG=(1d0-g**2)/((1d0-2d0*g*cos(theta)+g**2)**(3.0/2.0))/2d0

	n=1000
	tot=0d0
	HG=0d0
	do i=1,n
		t=theta+(real(i-1)/real(n-1)-0.5)*3.1415926536/180d0
		f=(1d0-g**2)/((1d0-2d0*g*cos(t)+g**2)**(3.0/2.0))/2d0
		HG=HG+f*sin(t)
		tot=tot+sin(t)
	enddo
	HG=HG/tot

	return
	end

c-----------------------------------------------------------------------
c  This subroutine reads in a file for the temperate-dependent opacities.
c  The files contains 2*nopac lines. The first line contains the temperature,
c  the second line the location of the .particle files, and so on.
c-----------------------------------------------------------------------
	
	subroutine readTopac(file,p)
	use Parameters
	IMPLICIT NONE
	character*500 file,partfile
	type(Particle) p
	integer ii,iopac
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	write(*,'("Read in temperature dependent opacities")')
	write(9,'("Read in temperature dependent opacities")')


	open(unit=50,file=file,RECL=6000)

	do iopac=1,p%nopac
 1	   read(50,*) p%Topac(iopac)
	   read(50,*) partfile

	   if (iopac .gt. 1) then
	      if (p%Topac(iopac) .le. p%Topac(iopac-1)) then
		 write(*,'("Temperatures should be in increasing order:")')
		 write(*,'("T[",i2,"]=",i5," <!< T[",i2,"]=",i5)')
     &                 iopac-1,p%Topac(iopac-1),iopac,p%Topac(iopac)
		 stop
	      endif
	   endif

	   call readparticle(partfile,p,.false.,2d0,2d0,0d0,1d0,iopac)
	   
	enddo
2	close(unit=50)

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	return
	end
	
	


	subroutine MixedAggregates(file,p)
	use Parameters
	IMPLICIT NONE
	character*500 file,partfile
	type(Particle) p
	integer ii,iopac
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	write(*,'("Creating mixed aggregate")')
	write(9,'("Creating mixed aggregate")')


	open(unit=50,file=file,RECL=6000)

	do iopac=1,p%nopac
 1	   read(50,*) p%cryst(iopac)
	   read(50,*) partfile

	   if (iopac .gt. 1) then
	      if (p%cryst(iopac) .gt. p%cryst(iopac-1)) then
		 write(*,'("crystallinity should be in decreasing order:")')
		 write(*,'("cryst[",i2,"]=",i5," >!> cryst[",i2,"]=",i5)')
     &                 iopac-1,p%cryst(iopac-1),iopac,p%cryst(iopac)
		 stop
	      endif
	   endif

	   call readparticle(partfile,p,.false.,2d0,2d0,0d0,1d0,iopac)
	enddo
2	close(unit=50)


	return
	end
	

	subroutine MakeMixAggregates(jj)
	use Parameters
	IMPLICIT NONE
	real*8 w0,w,f
	integer k,i,j,jj
	
	if(Grain(jj)%Tcryst.lt.real(TMAX)*dT) then
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		if(C(i,j)%T.lt.Grain(jj)%Tcryst) then
			Grain(jj)%Rcryst=D%R(i)
			goto 1
		endif
	enddo
	enddo
	endif

1	continue

	write(*,'("Mixing powerlaw start: ",f10.3," AU")') Grain(jj)%Rcryst
	write(9,'("Mixing powerlaw start: ",f10.3," AU")') Grain(jj)%Rcryst

	do i=0,D%nR-1
	do j=1,D%nTheta-1
c		if(C(i,j)%T.ge.Grain(jj)%Tcryst) then
c			f=1d0
c		else
			f=(Grain(jj)%Rcryst*AU/D%R_av(i))**Grain(jj)%powcryst
c		endif
		if(f.gt.1d0) f=1d0
		C(i,j)%wopac(jj,1:Grain(jj)%nopac)=0d0
		do k=1,Grain(jj)%nopac-1
			if(Grain(jj)%cryst(k).ge.f.and.Grain(jj)%cryst(k+1).lt.f) goto 3
		enddo
		k=Grain(jj)%nopac-1
3		continue
		C(i,j)%wopac(jj,k)=(Grain(jj)%cryst(k+1)-f)/(Grain(jj)%cryst(k+1)-Grain(jj)%cryst(k))
		C(i,j)%wopac(jj,k+1)=(f-Grain(jj)%cryst(k))/(Grain(jj)%cryst(k+1)-Grain(jj)%cryst(k))
	enddo
	enddo

	return
	end


	subroutine solve_mu0(r,mu,mu0)
	IMPLICIT NONE
	real*8 r,mu,mu0,min,mu1,d
	integer n,i
	n=1000
	
	min=1d200
	do i=1,n
		mu1=real(i-1)/real(n-1)
		d=(1d0-mu1**2)-r*(1d0-mu/mu1)
		if(dabs(d).lt.min) then
			mu0=mu1
			min=dabs(d)
		endif
	enddo
	
	return
	end
	


c-----------------------------------------------------------------------
c subroutine that only reads in the entire file, no regridding
c-----------------------------------------------------------------------
	subroutine readspecHR(starfile)
	use Parameters
	IMPLICIT NONE
	character*500 starfile
	real*8 x,y
	integer i
	
	open(unit=20,file=starfile,RECL=6000)
	nlamHR=0
1	read(20,*,end=2,err=1) x,y
	nlamHR=nlamHR+1
	goto 1
2	close(unit=20)

	allocate(lamHR(nlamHR))
	allocate(FstarHR(nlamHR))

	open(unit=20,file=starfile,RECL=6000)
	do i=1,nlamHR
3		read(20,*,err=3) lamHR(i),FstarHR(i)
	enddo
	close(unit=20)
	
	return
	end
		
	
	
	subroutine ConnectZones()
	use Parameters
	IMPLICIT NONE
	integer izone,iter,i1,i2,ir,nr
	real*8 Mtot,M1,M2,R,S1,S2,SdensZone,RR1,RR2,RR

	do izone=1,nzones
		if(Zone(izone)%Rconnect.lt.0d0) Zone(izone)%Rconnect=Zone(izone)%Rout
	enddo

	nr=250
	do iter=1,100
		do izone=1,nzones
			if(Zone(izone)%Mconnect.gt.0) then
				i1=izone
				i2=Zone(izone)%Mconnect
				R=Zone(izone)%Rconnect
				M1=0d0
				M2=0d0
				do ir=1,nr
					RR1=Zone(i1)%Rin+(Zone(i1)%Rout-Zone(i1)%Rin)*real(ir-1)/real(nr)
					RR2=Zone(i1)%Rin+(Zone(i1)%Rout-Zone(i1)%Rin)*real(ir)/real(nr)
					RR=sqrt(RR1*RR2)
					S1=SdensZone(RR,Zone(i1)%denspow,Zone(i1)%Rexp,Zone(i1)%gamma_exp)
					M1=M1+S1*(RR2-RR1)*RR

					RR1=Zone(i2)%Rin+(Zone(i2)%Rout-Zone(i2)%Rin)*real(ir-1)/real(nr)
					RR2=Zone(i2)%Rin+(Zone(i2)%Rout-Zone(i2)%Rin)*real(ir)/real(nr)
					RR=sqrt(RR1*RR2)
					S2=SdensZone(RR,Zone(i2)%denspow,Zone(i2)%Rexp,Zone(i2)%gamma_exp)
					M2=M2+S2*(RR2-RR1)*RR
				enddo
				S1=SdensZone(R,Zone(i1)%denspow,Zone(i1)%Rexp,Zone(i1)%gamma_exp)
				S2=SdensZone(R,Zone(i2)%denspow,Zone(i2)%Rexp,Zone(i2)%gamma_exp)
				M1=M1/S1
				M2=M2/S2
				Mtot=Zone(i1)%Mdust+Zone(i2)%Mdust
				Zone(i1)%Mdust=Mtot*M1/(M1+M2)
				Zone(i2)%Mdust=Mtot*M2/(M1+M2)
			endif
			if(Zone(izone)%Sconnect.gt.0) then
				i1=izone
				i2=Zone(izone)%Mconnect
				R=Zone(izone)%Rconnect
				Zone(i2)%Rsh=R
				Zone(i2)%sh=Zone(i1)%sh*(R/Zone(i1)%Rsh)**Zone(i1)%shpow
			endif
		enddo
	enddo
	
	return
	end


	real*8 function SdensZone(R,denspow,Rexp,gamma)
	IMPLICIT NONE
	real*8 R,denspow,Rexp,gamma
	
	SdensZone=R**(-denspow)*exp(-(R/Rexp)**gamma)
	
	return
	end
	

