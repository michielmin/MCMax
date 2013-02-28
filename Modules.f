	module Parameters
	IMPLICIT NONE
	integer nlam,idum,TMAX,MAXOBS,NPHISCATT,NTQHP
	real*8 pi,AU,Rsun,Msun,Lsun,sigma,parsec,kb,dT,RWmax,gas2dust,totaldist
	parameter(TMAX=10000,dT=1d0,NTQHP=100)
	parameter(MAXOBS=100,NPHISCATT=180)
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(AU=1.49598e13)
	parameter(parsec=3.08568025e18)
	parameter(Rsun=6.955e10,Msun=1.98892e33,Lsun=3.827e33)
	parameter(kb=1.3806503d-16,sigma=5.6704d-5)
c	parameter(gas2dust=100d0) ! Gijsexp, need it to be variable
	real*8,allocatable :: lam(:),BB(:,:),dnu(:),nu(:),IRF(:)
	real*8,allocatable :: shscale(:),column(:,:),specemit(:),muRad(:),Kext_column(:)
	real*8 xsf(NPHISCATT),ysf(NPHISCATT),zsf(NPHISCATT),tautot,epsiter,tau_max,nEJv
	real*8 xsn(NPHISCATT),ysn(NPHISCATT),zsn(NPHISCATT),dTDiffuse,factRW,tgridqhp(NTQHP)
	real*8 xin(NPHISCATT),yin(NPHISCATT),zin(NPHISCATT),f_weight,tauincrease,TdesQHP
	real*8 sin2phi(0:360),cos2phi(0:360),RmaxRefine,KDext(0:TMAX),KDabs(0:TMAX),BBint(0:TMAX)
	real*8 alphavis,dimstar,alphavispow,E_IRF
	real*8 lifetime,alphaturb,qturb !Gijsexp
	real*8 deadalpha,deadheight !Gijsexp: deadzone in midplane
	real*8 gsd_rmin,gsd_rmax,gsd_xi,gsd_vfrag !Gijsexp
	real*8 mrn_rmin,mrn_rmax,mrn_index !Gijsexp
	real*8 tau1_lam(100) !Gijsexp
	integer*8 nemit,ninteract,maxinteract,nmaxinteract
	integer ngrains,nobs,scat_how,NPhotDiffuse,maxiter,nRWinter,iRWinter,NphotFinal
	integer nspan,nlev,ntspan,ntlev,nexits,nruns,iTD(0:TMAX),NsigDiskstructure,nBW,nqhp
	integer nplanets,niter0,nspike,ngrains2,nzones,maxruntime
	integer gsd_diag !Gijsexp
	integer ntau1_lam !Gijsexp
	integer mrn_ngrains,thinparticle !Gijsexp
	logical struct_iter,scattering,arraysallocated,RNDW,dosmooth,use_obs_TMC,exportProDiMo
	logical FLD,storescatt,overflow,tcontact,tdes_iter,shell1D,forcediff,multiwav,outputfits
	logical useobspol,readmcscat,makeangledependence,gridrefine,etrace,use_qhp,use_topac,computeLRF
	logical tracestar,traceemis,tracescat,tracegas,radpress,haloswitch,raditer,viscous,computeTgas,getalpha
	logical fastviscous,convection,outfluxcontr,forcefirst,use_IRF,useTgas,g2d_heat,Tsmooth,emptylower
	logical scset,scsetsave,scseteq,mpset,mpstr ! Gijsexp
	logical fixmpset
	logical gsd,gsd_full,gsd_plot		!Gijsexp
	logical mrn		!Gijsexp
	logical deadzone	!Gijsexp
	logical topac_interpol	!Gijsexp
	logical,allocatable :: scattcomputed(:)
	integer,allocatable :: nscattcomputed(:)
	character*500 outdir,particledir
	real*8,allocatable :: coolingtime(:)
	integer,allocatable :: ncoolingtime(:)
	integer icoolingtime
	real*8 HotGasMinRad,HotGasMaxRad,HotGasT
	
	type Mueller
		real*8 F11(180),F12(180),F22(180)
		real*8 F33(180),F44(180),F34(180)
		real*8 IF11,IF12
	end type Mueller
	
	type Particle
		real*8 m,rho,Nc,rv,rvmin,rvmax
		real*8,allocatable :: Kabs(:,:),Ksca(:,:),Kext(:,:),shscale(:),g(:,:)
		real*8,allocatable :: Kp(:,:),Kpstar(:),Kpabsstar(:),Topac(:),cryst(:)
		real*8,allocatable :: KabsL(:),KextL(:)
		real*8 TdesA,TdesB
		real*8 Td_qhp,Mc	! Needed for QHP particles. Mc is the weight of one atom of the grain (in proton masses)
		type(Mueller),allocatable :: F(:,:)
		character*20 name
		logical settle,qhp,gascoupled,trace,force_vert_gf
		real*8 tdes_fast
		integer qhpnr
		character*20 material,shtype,roundtype
		integer nopac,parttype
		real*8 Rcryst,Tcryst,powcryst,maxtau,maxrad,minrad,shaperad
		real dust_moment1,dust_moment2,dust_moment3
		real*8 roundwidth,roundpow !Gijsexp, roundoff
! parttype=1	Normal particle
! parttype=2	Opacity file
! parttype=3	Mixed aggregates (still in beta)
! parttype=4	Quantum heated particle
! parttype=5	T dependent opacities
! parttype=6	Fits file particle
	end type Particle

	type Cell
		real*8 T,dens,mass,xedge(6),EJv,V,xedge2(6),Eabs
		integer Ni,iedge(6),lastphotnr,iTD,nrg
		real*8 Kabs,Ksca,Kext,Albedo,dT,ddens,gasfrac,EJv2,dEJv,TMC
		real*8 tauexit,dens0,KDext,KDabs,gasdens,HighTQHP,KDQHP,Tav
		real*8,allocatable :: TP(:),EJvP(:),EJvQHP(:),EviscDirect(:),KabsTot(:),KscaTot(:)
		type(Mueller) F
		real*8,allocatable :: w(:),scattfield(:,:,:),w0(:),QHP(:,:),LRF(:)
		integer,allocatable :: nLRF(:)
		real*8,allocatable :: tdistr(:,:),Tqhp(:),EabsP(:),thetrg(:),Tphi(:)
		real*8,allocatable :: scattQ(:,:,:),scattU(:,:,:),scattV(:,:,:)
		logical thick,diff,randomwalk
		real*8 KextLRF,ILRF
		real*8 FradR,FradZ
		real*8 densscale
		real*8,allocatable :: wopac(:,:),FE(:)
		real*8 KappaGas
		real*8 Tgas,Egas
		logical useFE,opacity_set
		real*8,allocatable :: line_emis(:),line_abs(:),velo_T(:)
	end type Cell

	type Disk
		real*8 Tstar,Rstar,Lstar,Mstar,Tstar2,Rstar2
		real*8 Mtot,Vtot,distance,PA,IA,Mdot,mu0max,Av
		real*8 denspow,Rin,Rout,shpow,sh1AU,denspow2,Rpow2,Rpow3,Rexp
		real*8,allocatable :: Fstar(:),theta_av(:),R_av(:),Rfix(:)
		real*8,allocatable :: R(:),Theta(:),thet(:),SinTheta(:)
		integer nR,nTheta,nRfix,ngap
		real*8,allocatable :: gap(:),gap1(:),gap2(:),gapshape(:)
		real*8,allocatable :: gaproundpow(:)
		character*20,allocatable :: gaproundtype(:)
	end type Disk

	type DiskZone
		real*8 Rin,Rout,Mdust,denspow,shpow,sh,Rsh
		real*8 a_min,a_max,a_pow,Rexp,gamma_exp
		logical fix_struct,sizedis
		real*8,allocatable :: abun(:)
		logical,allocatable :: inc_grain(:)
	end type DiskZone
	
	type(DiskZone),allocatable :: Zone(:)

	type(Cell),allocatable :: C(:,:)

	type(Disk) D

	type(Particle),allocatable :: Grain(:)
	
	type Photon
		real*8 x,y,z,vx,vy,vz,lam,nu,E
		integer ilam1,ilam2,edgeNr,i,j
		real*8 wl1,wl2
		logical onEdge,scatt,pol,viscous
		integer nr,irg
		real*8 fnr,Q,U,V,Sx,Sy,Sz
	end type Photon

	type path
		real*8 x,y,z,vx,vy,vz
		real*8,allocatable :: v(:),phi1(:),phi2(:)
		real*8,allocatable :: velo1(:),velo2(:)
		integer,allocatable :: i(:),j(:),jphi1(:),jphi2(:),k(:),irg(:)
		integer n
		logical hitstar
	end type path

	type RPhiImage
		real*8,allocatable :: image(:,:),R(:),Phi(:)
		real*8,allocatable :: imageQ(:,:),imageU(:,:),imageV(:,:)
		real*8 angle,lam,flux,rscale,zscale
		integer nr,nphi,scaletype
		type(path),allocatable :: p(:,:)
	end type RPhiImage

	type Telescope
		character*100 kind,flag
		integer nphi,Nphot,NphotAngle,nr,nexits,nt,nstar
		integer nbaseline,nangle,npixel,nint,nfov,scaletype,nlam ! Gijsexp
		real*8 lam1,lam2,angle,D,dlam,texp,D2,spider
		real*8 angle1,angle2,width,opening,mask,wmask,iwa,owa,strehl
		real*8 Ptelescope,APtelescope,snoise,RON,iprad
		real*8,allocatable :: b(:),theta(:),fov(:),lam(:) ! Gijsexp
		logical usepol,readmcscat,traceinverse,fluxcontr
		logical,allocatable :: trace(:)
		logical tracestar,traceemis,tracescat,tracegas
		character*500 psffile,linefile,popfile
		real*8 dvelo,abun
		integer nvelo,trans_nr1,trans_nr2
	end type Telescope

	type ExoPlanet
		character*500 name,file
		real*8 R,d,phi,A,P,T,theta
		real*8 E,Q,U,V
		integer i,j,k
		type(Mueller) F
	end type ExoPlanet

	type(ExoPlanet),allocatable :: Planets(:)

	end module Parameters
	
	
	module NAG
!c	use f90_unix_io
!c	use f90_unix_proc
	end module NAG
	
	subroutine get_command_argument(i,in)
	IMPLICIT NONE
	integer i
	character*500 in
	call getarg(i,in)
	return
	end
	
