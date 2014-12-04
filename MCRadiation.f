!c-----------------------------------------------------------------------
!c This subroutine computes the temperature structure of the disk
!c using Monte Carlo radiative transfer with Nphot photons.
!c The disk has to be initialized using the subroutine initialize
!c-----------------------------------------------------------------------
	subroutine MCRadiation(NphotTot,niter,dofastvisc)
	use NAG
	use Parameters
	use InputOutput
	IMPLICIT NONE
	integer i,j,Nphot,NphotTot,nextTout,ii,iruns,niter,iT,iopac
	real*8 ran2,tau,ct,inp,x,y,z,r,theta,increaseT
	logical escape,hitstar
	type(Photon) phot,phot2
!c-----------------------------------------
!c this is for the output MC spectrum
!c-----------------------------------------
	integer iangle,nangle
	real*8 clight
	parameter(nangle=30,clight=2.9979d8)
	real*8 angle(0:nangle),wangle(0:nangle),spec(nlam,0:nangle),cosangle(0:nangle)
	real*8 wlam(nlam),scatspec(nlam,nangle),spectemp(nlam),tot,EUV
	real*8 dspec(nlam),spec2(nlam,nangle),Mdot,FracVis,starttime,checktime
	integer ispec(nlam,nangle),iscatspec(nlam,nangle)
	character*500 specfile,pressurefile,timefile
	real*8,allocatable :: EJvTot(:,:),EJv2Tot(:,:),EJvTotP(:,:,:),betadisk(:,:,:)
	logical scatbackup,emitted,dofastvisc,phot_irf
	integer nsplit,isplit,iter,l
	real*8 split(2),determineT,determineTP,T,ShakuraSunyaevIJ,where_emit,FracIRF

	real*8 mu,G,Rad,phi,Evis,Efrac(0:D%nR,0:D%nTheta),Er,Sig,Einner,R1,R2
	real*8 F(ngrains),maxT,minT,determinegasfrac,Tevap,A(ngrains),ExtISM(nlam)
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
	parameter(G=6.67300d-8) ! in cm^3/g/s
	real*8 Reddening,compute_dlam
	real*8 WeightedAlpha,obs_opening

	real*8 radtau,radialtau,tau1,Av
	integer ntau1

6	continue

	Av=D%Av
	if(adjustAv) then
		radtau=radialtau(0.55d0,tau1,ntau1,1)
		Av=Av+2.5*log10(exp(-radtau))
		if(Av.lt.0d0) Av=0d0
		write(*,'("Av adjusted to: ",f5.2)') Av
	endif
	do i=1,nlam
		ExtISM(i)=Reddening(lam(i),compute_dlam(lam(i)),Av)
	enddo

	do i=1,nlam
	spec(i,0)=0d0
	do j=1,nangle
		spec(i,j)=0d0
		spec2(i,j)=0d0
		ispec(i,j)=0
		scatspec(i,j)=0d0
		iscatspec(i,j)=0
	enddo
	if(i.eq.1) then
		wlam(i)=-1d23/(nu(2)/2d0)
	else if(i.eq.nlam) then
		wlam(i)=-1d23/((nu(nlam)-nu(nlam-1))/2d0)
	else
		wlam(i)=-1d23/((nu(i+1)-nu(i-1))/2d0)
	endif
	enddo

	do i=0,nangle
		cosangle(i)=dfloat(i)/dfloat(nangle)
		angle(i)=acos(cosangle(i))
	enddo
	do i=1,nangle
		wangle(i)=1d0/(cosangle(i)-cosangle(i-1))
	enddo
	obs_opening=1d0+4d0*cos(D%IA*pi/180d0)
	wangle(0)=1d0/(cos(pi*(D%IA-obs_opening)/180d0)-cos(pi*(D%IA+obs_opening)/180d0))

c====================================================================================
c include the inner disk like in the Pringle 1981, Akeson 2005 papers
c====================================================================================

	if (deadzone.or.gravstable.or.D%Rexp.lt.1d10) then
	   Mdot=D%MdotR(1) ! accretion at inner disk edge
	else 
	   Mdot=D%Mdot
	endif

	Einner=0d0
	if(inner_gas) then
		Einner=G*D%Mstar*Mdot/(4d0*D%Rstar)
		Einner=Einner*(1d0+2d0*(D%Rstar/(D%R(1)*AU))**(3d0/2d0)-3d0*(D%Rstar/(D%R(1)*AU)))
		Einner=Einner/(4d0*pi)

c EUV is the energy from the innermost region, probably comes out in the UV
		EUV=G*D%Mstar*D%Mdot/(4d0*D%Rstar)
		EUV=EUV*(1d0+2d0*(1d0/(Rinner_gas))**(3d0/2d0)-3d0*(1d0/(Rinner_gas)))
		EUV=EUV/(4d0*pi)

		Einner=Einner-EUV
	endif

c====================================================================================
c  Compute the viscous heating per cell
c====================================================================================
	
	Evis=0d0
	if(viscous) then
	do i=1,D%nR-1
		if(.not.allocated(C(i,j)%EviscDirect)) allocate(C(i,j)%EviscDirect(ngrains))
		if(.not.allocated(C(i,j)%FE)) then
			allocate(C(i,j)%FE(0:ngrains))
			C(i,j)%useFE=.false.
		endif
		if (deadzone.or.gravstable.or.D%Rexp.lt.1d10) then
		   Mdot=D%MdotR(i) 
		else 
		   Mdot=D%Mdot
		endif

		! Viscous heating in cells i
		Er=(2d0*sqrt(D%Rstar/(AU*D%R(i+1)))-3d0)/(3d0*D%R(i+1)*AU)
		Er=Er-(2d0*sqrt(D%Rstar/(AU*D%R(i)))-3d0)/(3d0*D%R(i)*AU)
		Er=2d0*pi*Er*3d0*G*D%Mstar*Mdot/(4d0*pi)
		Er=Er/(4d0*pi)

                !  Calculate surface density and density-weighted viscosity (in case of deadzones)
		Sig=0d0
		WeightedAlpha=0d0
		do j=1,D%nTheta-1
			Sig=Sig+C(i,j)%gasdens*C(i,j)%V
			WeightedAlpha=WeightedAlpha+C(i,j)%gasdens*C(i,j)%V*C(i,j)%alphavis
		enddo
		WeightedAlpha=WeightedAlpha/Sig

		!  Verical distribution of viscous energy proportional to:
                !   density or dens*alpha (in case of deadzone)
		!  Temperature is ignored! (for stability?)
		do j=1,D%nTheta-1
			Efrac(i,j)=Er*C(i,j)%gasdens*C(i,j)%V/Sig *C(i,j)%alphavis/WeightedAlpha
			do ii=1,ngrains
				C(i,j)%EviscDirect(ii)=0d0
			enddo
			Evis=Evis+Efrac(i,j)
		enddo
	enddo

c	Er=(2d0*sqrt(D%Rstar/(AU*D%R(D%nR)))-3d0)/(3d0*D%R(D%nR)*AU)
c	Er=Er-(2d0*sqrt(D%Rstar/(AU*D%R(1)))-3d0)/(3d0*D%R(1)*AU)
c	Er=2d0*pi*Er*3d0*G*D%Mstar*Mdot/(4d0*pi)
c	print*,100d0*G*D%Mstar*D%Mdot/(2d0*D%Rstar)/(4d0*sigma*D%Tstar**4*pi*D%Rstar**2)
c	print*,100d0*(Er/(4d0*pi))/(D%Lstar+Er/(4d0*pi))

	if(dofastvisc) Evis=0d0
	endif

	nemit=0
	nmaxinteract=0

	call reset2d()
	tautot=0d0
	totaldist=0d0

	tau=0d0
	do j=1,D%nTheta-1
	do i=1,D%nR-1
		do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			tau=tau+C(i,j)%dens*(D%R(i+1)-D%R(i))*AU*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
     &				*Grain(ii)%Kpstar(iopac)*(D%thet(j+1)-D%thet(j))*2d0/pi
		enddo
		enddo
	enddo
	enddo
	if(forcefirst) then
		write(*,'("Forcing first interaction")')
		write(9,'("Forcing first interaction")')
	endif
	
	Nphot=NphotTot/nruns
	
	allocate(EJvTot(0:D%nR-1,1:D%nTheta-1))
	allocate(EJv2Tot(0:D%nR-1,1:D%nTheta-1))
	if(.not.tcontact.or.tdes_iter) then
		allocate(EJvTotP(1:ngrains,0:D%nR-1,1:D%nTheta-1))
	endif
	do j=1,D%nTheta-1
	do i=0,D%nR-1
		EJvTot(i,j)=0d0
		EJv2Tot(i,j)=0d0
		if(.not.tcontact.or.tdes_iter) then
			do ii=1,ngrains
				EJvTotP(ii,i,j)=0d0
			enddo
		endif
		C(i,j)%FradZ=0d0
		C(i,j)%FradR=0d0
		if(.not.allocated(C(i,j)%KabsTot)) then
			allocate(C(i,j)%KabsTot(nlam))
			allocate(C(i,j)%KscaTot(nlam))
		endif
		do l=1,nlam
			C(i,j)%KabsTot(l)=0d0
			C(i,j)%KscaTot(l)=0d0
			do ii=1,ngrains
				if(.not.Grain(ii)%qhp) then
					do iopac=1,Grain(ii)%nopac
						C(i,j)%KabsTot(l)=C(i,j)%KabsTot(l)
     &	+Grain(ii)%Kabs(iopac,l)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
					enddo
				endif
				do iopac=1,Grain(ii)%nopac
					C(i,j)%KscaTot(l)=C(i,j)%KscaTot(l)
     &	+Grain(ii)%Ksca(iopac,l)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
				enddo
			enddo
		enddo
		C(i,j)%opacity_set=.true.
		C(i,j)%fQHP=0d0
		do ii=1,ngrains
			if(Grain(ii)%qhp) C(i,j)%fQHP=C(i,j)%fQHP+C(i,j)%w(ii)
		enddo
	enddo
	enddo


	FracVis=(Evis+Einner)/(D%Lstar+Evis+E_IRF+Einner)
	if(FracVis.gt.0.5d0) then
		FracVis=0.5d0
		write(*,'("Increasing statistics from the star to ",f6.2,"%")') 100d0*(1d0-FracVis)
		write(9,'("Increasing statistics from the star to ",f6.2,"%")') 100d0*(1d0-FracVis)
	endif
	if(use_IRF) then
		FracIRF=E_IRF/(D%Lstar+Evis+E_IRF)
		write(*,'("Energy: ",f7.3,"% from the IRF")') 100d0*E_IRF/(D%Lstar+Evis+E_IRF+Einner)
		write(9,'("Energy: ",f7.3,"% from the IRF")') 100d0*E_IRF/(D%Lstar+Evis+E_IRF+Einner)
		if(FracIRF.lt.0.01) then
			FracIRF=0.01
		endif
	else
		FracIRF=0d0
	endif

	if(viscous.or.inner_gas) then
		write(*,'("Energy: ",f7.3,"% from the star")') 100d0*D%Lstar/(D%Lstar+Evis+E_IRF+Einner)
		write(*,'("Energy: ",f7.3,"% from the disk")') 100d0*(Einner+Evis)/(D%Lstar+Evis+E_IRF+Einner)
		write(9,'("Energy: ",f7.3,"% from the star")') 100d0*D%Lstar/(D%Lstar+Evis+E_IRF+Einner)
		write(9,'("Energy: ",f7.3,"% from the disk")') 100d0*(Einner+Evis)/(D%Lstar+Evis+E_IRF+Einner)
		write(*,'("Energy: ",f7.3,"% not included (innermost region)")') 100d0*(EUV)/(D%Lstar+Evis+E_IRF+Einner+EUV)
		write(9,'("Energy: ",f7.3,"% not included (innermost region)")') 100d0*(EUV)/(D%Lstar+Evis+E_IRF+Einner+EUV)
	endif
	
	write(*,'("Emitting ",i10," photon packages")') NphotTot
	write(9,'("Emitting ",i10," photon packages")') NphotTot
	call cpu_time(starttime)

	do iruns=1,nruns
	
	if(nruns.gt.1) then
		write(*,'("Run number",i4," of",i4)') iruns,nruns
		write(9,'("Run number",i4," of",i4)') iruns,nruns
	endif

	do j=1,D%nTheta-1
	do i=0,D%nR-1
		C(i,j)%T=0d0
		C(i,j)%EJv=0d0
		C(i,j)%EJv2=0d0
		C(i,j)%Eabs=0d0
		C(i,j)%Egas=0d0
		C(i,j)%lastphotnr=0
		if(.not.tcontact.or.tdes_iter) then
			do ii=1,ngrains
				C(i,j)%TP(ii)=0d0
				C(i,j)%EJvP(ii)=0d0
				C(i,j)%EabsP(ii)=0d0
			enddo
		endif
		if(use_qhp) then
			do ii=1,ngrains
				if(Grain(ii)%qhp) C(i,j)%EJvQHP(Grain(ii)%qhpnr)=0d0
			enddo
		endif
		if(computeLRF) then
			C(i,j)%LRF(1:nlam)=0d0
			C(i,j)%nLRF(1:nlam)=0
		endif
		C(i,j)%KextLRF=0d0
		C(i,j)%ILRF=0d0
		call CheckMinimumDensity(i,j)
	enddo
	enddo
	ncoolingtime=0
	coolingtime=0d0
	
	do i=1,Nphot
		call tellertje(i,Nphot)

		if(maxruntime.gt.0.and.Nphot.gt.10) then
			call cpu_time(checktime)
			checktime=checktime-starttime
			if((checktime.gt.real(maxruntime)
     &	.or.(i.gt.(Nphot/10).and.(checktime*real(Nphot)/real(i)).gt.real(5*maxruntime)))
     &	.and.(real(i)/real(Nphot).lt.0.9)) then
				write(*,'("STOPPING DUE TO TIME CONSTRAINT!!")')
				write(*,'("REDUCING NUMBER OF PHOTON PACKAGES!!")')
				write(9,'("STOPPING DUE TO TIME CONSTRAINT!!")')
				write(9,'("REDUCING NUMBER OF PHOTON PACKAGES!!")')
				NphotTot=NphotTot*(real(maxruntime)/checktime)*(0.5d0*real(i)/real(Nphot))
				if(NphotTot.lt.10) NphotTot=10
				deallocate(EJvTot)
				deallocate(EJv2Tot)
				if(.not.tcontact.or.tdes_iter) then
					deallocate(EJvTotP)
				endif
				goto 6
			endif
		endif
				
		phot%nr=i
		phot%fnr=real(i)/real(Nphot)

3		continue
		ninteract=0
		iRWinter=0
		nEJv=0d0


		where_emit=ran2(idum)
		phot_irf=.false.
		if(where_emit.gt.(FracVis+FracIRF)) then
c emit from the star
			phot%viscous=.false.
			call randomdirection(phot%x,phot%y,phot%z)
			phot%x=D%R(0)*phot%x
			phot%y=D%R(0)*phot%y
			phot%z=D%R(0)*phot%z

			call emit(phot,D%Fstar,D%Lstar)

			phot%E=D%Lstar/(real(Nphot)*(1d0-FracVis-FracIRF))

			call randomdirection(phot%vx,phot%vy,phot%vz)
			if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.lt.0d0) then
				phot%vx=-phot%vx
				phot%vy=-phot%vy
				phot%vz=-phot%vz
			endif

			phot%edgeNr=1
			phot%onEdge=.true.
			ct=abs(phot%z)/D%R(0)
			do j=1,D%nTheta-1
				if(ct.lt.D%Theta(j).and.ct.gt.D%Theta(j+1)) then
					phot%i=0
					phot%j=j
				endif
			enddo

		else if(where_emit.gt.FracIRF) then
!c emit from viscous heating
			phot%viscous=.true.
			Er=(Evis+Einner)*ran2(idum)
			if(Er.lt.Evis) then
!c emit from the viscous heating inside the dust disk
9			continue
			do ii=1,D%nR-1
				do j=1,D%nTheta-1
					Er=Er-Efrac(ii,j)
					if(Er.lt.0d0) goto 8
				enddo
			enddo
			goto 9
		
8			continue

			Rad=ran2(idum)
			Rad=sqrt(D%R(ii)**2*Rad+D%R(ii+1)**2*(1d0-Rad))
			Theta=ran2(idum)
			Theta=D%Theta(j)*Theta+D%Theta(j+1)*(1d0-Theta)
			phot%z=Rad*Theta
			if(ran2(idum).lt.0.5) phot%z=-phot%z
			phi=ran2(idum)*pi*2d0
			phot%x=Rad*sqrt(1d0-Theta**2)*sin(phi)
			phot%y=Rad*sqrt(1d0-Theta**2)*cos(phi)
			call randomdirection(phot%vx,phot%vy,phot%vz)
			phot%i=ii
			phot%j=j
			phot%onEdge=.false.
			icoolingtime=ii
			ncoolingtime(ii)=ncoolingtime(ii)+1

			call randomdirection(phot%vx,phot%vy,phot%vz)

			phot%E=(Evis+Einner)/(real(Nphot)*FracVis)

c emit the viscous photon
			call EmitViscous(phot)

			else
c emit from the inner gas disk (Pringle (1981), Akeson (2005)
			phot%viscous=.false.

			Er=ran2(idum)*((2d0*(D%Rstar/(D%R(1)*AU))**(3d0/2d0)-3d0*(D%Rstar/(D%R(1)*AU)))-
     &					  (2d0*(1d0/(Rinner_gas))**(3d0/2d0)-3d0*(1d0/(Rinner_gas))))
			R1=D%Rstar*Rinner_gas/AU
			R2=D%R(1)

			Rad=R1+ran2(idum)*(R2-R1)
			do iter=1,10
				tot=(2d0*(D%Rstar/(Rad*AU))**(3d0/2d0)-3d0*(D%Rstar/(Rad*AU)))-
     &					  (2d0*(1d0/(Rinner_gas))**(3d0/2d0)-3d0*(1d0/(Rinner_gas)))
				if(tot.gt.Er) then
					R2=Rad
				else
					R1=Rad
				endif
				Rad=(R2+R1)/2d0
			enddo

			j=D%nTheta-1
			Theta=0.5
			Theta=D%Theta(j)*Theta+D%Theta(j+1)*(1d0-Theta)
			phot%z=Rad*Theta
			if(ran2(idum).lt.0.5) phot%z=-phot%z
			phi=ran2(idum)*pi*2d0
			phot%x=Rad*sqrt(1d0-Theta**2)*sin(phi)
			phot%y=Rad*sqrt(1d0-Theta**2)*cos(phi)
			call randomdirection(phot%vx,phot%vy,phot%vz)
			phot%i=0
			phot%j=j
			phot%onEdge=.false.
			icoolingtime=0
			ncoolingtime(0)=ncoolingtime(0)+1

			call randomdirection(phot%vx,phot%vy,phot%vz)

			phot%E=(Evis+Einner)/(real(Nphot)*FracVis)

c emit the viscous photon
			call EmitViscous(phot)

			endif
		else
c emit from the interstellar radiation field
			phot%viscous=.false.
			phot_irf=.true.
			call randomdirection(phot%x,phot%y,phot%z)
			phot%x=D%R(D%nR)*phot%x
			phot%y=D%R(D%nR)*phot%y
			phot%z=D%R(D%nR)*phot%z

			call emit(phot,IRF,E_IRF)

			phot%E=E_IRF/(real(Nphot)*FracIRF)

7			call randomdirection(phot%vx,phot%vy,phot%vz)
			if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.gt.0d0) then
				phot%vx=-phot%vx
				phot%vy=-phot%vy
				phot%vz=-phot%vz
			endif
			if(abs(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz)/D%R(D%nR).lt.ran2(idum)) goto 7

			phot%edgeNr=2
			phot%onEdge=.true.
			ct=abs(phot%z)/D%R(D%nR)
			do j=1,D%nTheta-1
				if(ct.lt.D%Theta(j).and.ct.gt.D%Theta(j+1)) then
					phot%i=D%nR-1
					phot%j=j
				endif
			enddo
		endif

		if(scattering) then
			phot%pol=.false.
			phot%scatt=.false.
			phot%Q=0d0
			phot%U=0d0
			phot%V=0d0
			x=-phot%vy
			y=phot%vx
			z=0d0
			r=sqrt(x**2+y**2+z**2)
			phot%Sx=x/r
			phot%Sy=y/r
			phot%Sz=z/r
		endif
		phot2=phot
		nsplit=1
		split(1)=1d0
		if(forcefirst) then
			call tau2exit(phot,tau)
			if(tau.lt.1d-9) then
				nsplit=1
				split(1)=1d0
			else 
			if(tau.gt.1d-6) then
				nsplit=2
				split(1)=exp(-tau)
				split(2)=1d0-split(1)
			else
				nsplit=2
				split(1)=1d0-tau
				split(2)=tau
			endif
			endif
		endif

		do isplit=1,nsplit

		ninteract=0
		iRWinter=0
		nEJv=0d0
		if(forcefirst) call emit(phot,D%Fstar,D%Lstar)
		phot=phot2

		phot%E=phot2%E*split(isplit)
		if(forcefirst) then
			if(isplit.eq.1) then
				tau=1d200
				do ii=1,D%nR-1
					C(ii,phot%j)%Ni=C(ii,phot%j)%Ni-1
				enddo
			else
				tau=-log(ran2(idum)*split(2)+split(1))
			endif
		else
			tau=-log(ran2(idum))
		endif
1		continue

		call trace2d(phot,tau,escape,hitstar,.true.)
		if(hitstar) then
			if(forcefirst) goto 4
			goto 3
		endif
		if(escape) goto 2
		ninteract=ninteract+1
		phot_irf=.false.
		call interact(phot)
		if(overflow) then
			if(ran2(idum).lt.(1d0/real(maxinteract))) goto 4
			phot%E=phot%E*real(maxinteract)/real(maxinteract-1)
		endif

		tau=-log(ran2(idum))
		goto 1

2		continue
		if(phot_irf) goto 4
!c-----------------------------------------
!c this is for the output MC spectrum
!c-----------------------------------------
		if(phot%vz.gt.0d0) then
			iangle=phot%vz*real(nangle)+1
		else
			iangle=-phot%vz*real(nangle)+1
		endif
		if(multiwav) then
			spectemp(1:nlam)=0d0
			do ii=1,ngrains
				do iopac=1,Grain(ii)%nopac
					spectemp(1:nlam)=spectemp(1:nlam)+Grain(ii)%Kext(iopac,1:nlam)*column(ii,iopac)
				enddo
			enddo
			do j=1,nlam
				if(spectemp(j).lt.1000d0) then
					spectemp(j)=specemit(j)*exp(-spectemp(j))
				else
					spectemp(j)=0d0
				endif
			enddo
			call integrate(spectemp,tot)
			if(tot.gt.1d-100) then
				spectemp=spectemp/tot
				spectemp=1d23*wangle(iangle)*phot%E*spectemp
				spec(1:nlam,iangle)=spec(1:nlam,iangle)+spectemp(1:nlam)/real(nruns)
				ispec(1:nlam,iangle)=ispec(1:nlam,iangle)+1
				spec2(1:nlam,iangle)=spec2(1:nlam,iangle)+spectemp(1:nlam)**2
				if(phot%scatt) then
					scatspec(1:nlam,iangle)=scatspec(1:nlam,iangle)+spectemp(1:nlam)/real(nruns)
				endif
				if(abs(acos(abs(phot%vz))*180d0/pi-D%IA).lt.obs_opening) then
					spec(1:nlam,0)=spec(1:nlam,0)+spectemp(1:nlam)*wangle(0)/real(nruns)/wangle(iangle)
				endif
			else
				if(phot%wl1.gt.phot%wl2) then
					spec(phot%ilam1,iangle)=spec(phot%ilam1,iangle)+wangle(iangle)*phot%E*wlam(phot%ilam1)/real(nruns)
					ispec(phot%ilam1,iangle)=ispec(phot%ilam1,iangle)+1
					if(phot%scatt) then
						scatspec(phot%ilam1,iangle)=scatspec(phot%ilam1,iangle)+wangle(iangle)*phot%E*wlam(phot%ilam1)/real(nruns)
						iscatspec(phot%ilam1,iangle)=iscatspec(phot%ilam1,iangle)+1
					endif
				else
					spec(phot%ilam2,iangle)=spec(phot%ilam2,iangle)+wangle(iangle)*phot%E*wlam(phot%ilam2)/real(nruns)
					ispec(phot%ilam2,iangle)=ispec(phot%ilam2,iangle)+1
					if(phot%scatt) then
						scatspec(phot%ilam2,iangle)=scatspec(phot%ilam2,iangle)+wangle(iangle)*phot%E*wlam(phot%ilam2)/real(nruns)
						iscatspec(phot%ilam2,iangle)=iscatspec(phot%ilam2,iangle)+1
					endif
				endif
			endif
		else
			if(phot%wl1.gt.phot%wl2) then
				spec(phot%ilam1,iangle)=spec(phot%ilam1,iangle)+wangle(iangle)*phot%E*wlam(phot%ilam1)/real(nruns)
				ispec(phot%ilam1,iangle)=ispec(phot%ilam1,iangle)+1
				if(phot%scatt) then
					scatspec(phot%ilam1,iangle)=scatspec(phot%ilam1,iangle)+wangle(iangle)*phot%E*wlam(phot%ilam1)/real(nruns)
					iscatspec(phot%ilam1,iangle)=iscatspec(phot%ilam1,iangle)+1
				endif
			else
				spec(phot%ilam2,iangle)=spec(phot%ilam2,iangle)+wangle(iangle)*phot%E*wlam(phot%ilam2)/real(nruns)
				ispec(phot%ilam2,iangle)=ispec(phot%ilam2,iangle)+1
				if(phot%scatt) then
					scatspec(phot%ilam2,iangle)=scatspec(phot%ilam2,iangle)+wangle(iangle)*phot%E*wlam(phot%ilam2)/real(nruns)
					iscatspec(phot%ilam2,iangle)=iscatspec(phot%ilam2,iangle)+1
				endif
			endif
		endif
4		continue
	enddo
!c-----------------------------------------
!c-----------------------------------------
5	continue
	enddo

	print*,'Average distance travelled per photon:',totaldist/real(Nphot)

	do j=1,D%nTheta-1
	do i=0,D%nR-1
		EJvTot(i,j)=EJvTot(i,j)+C(i,j)%EJv/real(nruns)
		EJv2Tot(i,j)=EJv2Tot(i,j)+C(i,j)%EJv2/real(nruns)
		if(.not.tcontact.or.tdes_iter) then
			do ii=1,ngrains
				EJvTotP(ii,i,j)=EJvTotP(ii,i,j)+C(i,j)%EJvP(ii)/real(nruns)
			enddo
		endif
	enddo
	enddo

	enddo

	do j=1,D%nTheta-1
	do i=0,D%nR-1
		C(i,j)%EJv=EJvTot(i,j)
		C(i,j)%EJv2=EJv2Tot(i,j)
		if(.not.tcontact.or.tdes_iter) then
			do ii=1,ngrains
				C(i,j)%EJvP(ii)=EJvTotP(ii,i,j)
			enddo
		endif
		if(computeLRF) then
			C(i,j)%LRF(1:nlam)=C(i,j)%LRF(1:nlam)/C(i,j)%V
		endif
		if(.not.tcontact.or.tdes_iter) then
			C(i,j)%KextLRF=C(i,j)%KextLRF/C(i,j)%V
			C(i,j)%ILRF=C(i,j)%ILRF/C(i,j)%V
		endif
		if(C(i,j)%useFE.and.computeTgas.and.g2d_heat) then
			C(i,j)%EJv=C(i,j)%EJv*(1d0-C(i,j)%FE(0))
			if(.not.tcontact.or.tdes_iter) then
				do ii=1,ngrains
					C(i,j)%EJvP(ii)=C(i,j)%EJvP(ii)*(1d0-C(i,j)%FE(ii))
				enddo
			endif
			if(use_qhp) then
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
					do iopac=1,Grain(ii)%nopac
						C(i,j)%EJvQHP(Grain(ii)%qhpnr)=C(i,j)%EJvQHP(Grain(ii)%qhpnr)
     &	*(1d0-C(i,j)%FE(0))
					enddo
					endif
				enddo
			endif
		endif
	enddo
	enddo

	write(*,'("Number of photons emitted:                ",i12)') NphotTot
	write(*,'("Numer of interaction overflows:           ",i12)') nmaxinteract
	write(*,'("Average number of reemissions per photon: ",f12.4)') real(nemit-NphotTot)/real(NphotTot)
	write(*,'("Average optical depth:                    ",f12.4)') tautot/real(NphotTot)

	write(9,'("Number of photons emitted:                ",i12)') NphotTot
	write(9,'("Numer of interaction overflows:           ",i12)') nmaxinteract
	write(9,'("Average number of reemissions per photon: ",f12.4)') real(nemit-NphotTot)/real(NphotTot)
	write(9,'("Average optical depth:                    ",f12.4)') tautot/real(NphotTot)


	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	call flush(9)


	if(emptylower) then
		spec=(spec+1d23*D%Fstar(i))/2d0
		scatspec=scatspec/2d0
	endif

!c-----------------------------------------
!c this is for the output MC spectrum
!c-----------------------------------------
	write(*,'("Writing MC spectra")')
	write(*,'("Number of angles:     ",i10)') nangle
	write(9,'("Writing MC spectra")')
	write(9,'("Number of angles:     ",i10)') nangle
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	do j=1,nangle
	angle(j)=180d0*acos((dfloat(j)-0.5d0)/dfloat(nangle))/pi
	write(specfile,'(a,"MCSpec",i1,f3.1,".dat")') outdir(1:len_trim(outdir))
     &			,int((angle(j))/10d0),angle(j)-10d0*int((angle(j)/10d0))
	open(unit=20,file=specfile,RECL=1000)
	if(multiwav) then
		dspec(1:nlam)=(spec(1:nlam,j)
     &				+sqrt(abs(spec2(1:nlam,j)-spec(1:nlam,j)**2/real(ispec(1:nlam,j)))))
     &				/sqrt(real(ispec(1:nlam,j)))
	else
		dspec(1:nlam)=spec(1:nlam,j)/sqrt(real(ispec(1:nlam,j)))
	endif
	do i=1,nlam
		write(20,*) lam(i),spec(i,j)*ExtISM(i)/D%distance**2,dspec(i)*ExtISM(i)/D%distance**2
     &			,scatspec(i,j)*ExtISM(i)/D%distance**2,iscatspec(i,j),1d23*D%Fstar(i)/D%distance**2
	enddo
	close(unit=20)
	enddo

	write(specfile,'(a,"MCSpec.dat")') outdir(1:len_trim(outdir))
	open(unit=20,file=specfile,RECL=1000)
	do i=1,nlam
		write(20,*) lam(i),spec(i,0)*ExtISM(i)/D%distance**2
	enddo
	close(unit=20)


C	if(scattering) then
C	write(*,'("Writing scattered spectra")')
C	write(*,'("Number of angles:     ",i10)') nangle
C	write(9,'("Writing scattered spectra")')
C	write(9,'("Number of angles:     ",i10)') nangle
C	write(*,'("--------------------------------------------------------")')
C	write(9,'("--------------------------------------------------------")')
C	write(specfile,'(a,"MCScattered.dat")') outdir(1:len_trim(outdir))
C	open(unit=20,file=specfile,RECL=1000)
C	write(20,'(i10,"    : number of angles")') nangle
C	write(20,'(i10,"    : number of wavelengths")') nlam
C	do j=nangle,1,-1
C		angle(j)=180d0*acos((dfloat(j)-0.5d0)/dfloat(nangle))/pi
C		write(20,'(f10.4)') angle(j)
C	enddo
C	do i=1,nlam
C		write(20,*) lam(i),(scatspec(i,j)/D%distance**2,j=nangle,1,-1),(iscatspec(i,j),j=nangle,1,-1)
C	enddo
C	close(unit=20)
C	endif


	if(viscous) then
		write(timefile,'(a,"coolingtime.dat")') outdir(1:len_trim(outdir))
		open(unit=20,file=timefile,RECL=6000)
		do i=1,D%nR-1
			if(ncoolingtime(i).gt.0) then
				write(20,*) D%R_av(i)/AU,(coolingtime(i)*AU/real(ncoolingtime(i)))/9.46e17
			else
				write(20,*) D%R_av(i)/AU,0d0
			endif
		enddo
		close(unit=20)
	endif

!c-----------------------------------------
!c-----------------------------------------		
	
	

c	write(pressurefile,'(a,"RadPressure.dat")') outdir(1:len_trim(outdir))
c	open(unit=20,file=pressurefile,RECL=6000)
c	write(20,'("# Format number")')
c	write(20,'(i6)') 1
c	write(20,'("# NR, NT")')
c	write(20,'(i6,i6)') D%nR-1,D%nTheta-1
c	write(20,'("# Spherical radius grid [cm] (middle of cell)")')
c	do i=1,D%nR-1
c		write(20,*) D%R_av(i)
c	enddo
c	write(20,'("# Theta grid [rad, from pole] (middle of cell)")')
c	do i=1,D%nTheta-1
c		write(20,*) D%theta_av(i)
c	enddo
c	write(20,'("# Radial component (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
c	do i=1,D%nR-1
c		do j=1,D%nTheta-1
c			write(20,*) C(i,j)%FradR
c		enddo
c	enddo
c	write(20,'("# Vertical component (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
c	do i=1,D%nR-1
c		do j=1,D%nTheta-1
c			write(20,*) C(i,j)%FradZ
c		enddo
c	enddo
c	close(unit=20)

	do i=1,D%nR-1
	do j=1,D%nTheta-1
		C(i,j)%FradZ=C(i,j)%FradZ/(-(cos(D%theta_av(j))*C(i,j)%dens*C(i,j)%V*6.67300d-8*D%Mstar/D%R_av(i)**2))
		C(i,j)%FradR=C(i,j)%FradR/(-(sin(D%theta_av(j))*C(i,j)%dens*C(i,j)%V*6.67300d-8*D%Mstar/D%R_av(i)**2))
		C(i,j)%FradZ=C(i,j)%FradZ/gas2dust
		C(i,j)%FradR=C(i,j)%FradR/gas2dust
	enddo
	enddo

	if(radpress) then
		allocate(betadisk(D%nR-1,D%nTheta-1,2))
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			betadisk(i,j,1)=C(i,j)%FradZ*gas2dust
			betadisk(i,j,2)=C(i,j)%FradR*gas2dust
		enddo
		enddo
		if(outputfits) then
			write(pressurefile,'(a,"RadPressure.fits.gz")') outdir(1:len_trim(outdir))
		else
			write(pressurefile,'(a,"RadPressure.dat")') outdir(1:len_trim(outdir))
		endif
		call outputstruct(pressurefile,(/'ARRAY  '/),1,0,betadisk(1:D%nR-1,1:D%nTheta-1,1:2),2)
		deallocate(betadisk)
	endif
	
	if(dofastvisc.and.viscous) then
	do i=1,D%nR-1
		do j=1,D%nTheta-1
			T=ShakuraSunyaevIJ(i,j)
			iT=T/dT
			if(iT.gt.TMAX-1)iT=TMAX-1
			if(iT.lt.1)iT=1
			do ii=1,ngrains
				do iopac=1,Grain(ii)%nopac
					C(i,j)%EJv=C(i,j)%EJv+Grain(ii)%Kp(iopac,iT)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)*C(i,j)%V
				enddo
			enddo
			C(i,j)%T=determineT(phot)
			C(i,j)%Ni=1000
		enddo
	enddo

	endif

	do j=1,D%nTheta-1
	do i=0,D%nR-1
		deallocate(C(i,j)%KabsTot)
		deallocate(C(i,j)%KscaTot)
		C(i,j)%opacity_set=.false.
	enddo
	enddo

	return
	end
	
	

!c-----------------------------------------------------------------------
!c This subroutine traces a photon through the disk over an
!c optical depth tau0. The photon has to be initialized fully
!c On output the escape flag is set to .true. if the photon escaped
!c the disk. The flag hitstar is set to .true. if the photon hits
!c the central star. The flag locfield is an input flag to indicate
!c if the local field should be constructed in order to compute the
!c temperature structure afterwards.
!c 2008-07-11:	Bugfix with the KabsBBGrains parameter (which was wrong)
!c-----------------------------------------------------------------------
	subroutine trace2d(phot,tau0,escape,hitstar,locfield)
	use Parameters
	IMPLICIT NONE
	type(Photon) phot
	real*8 tau0,tau,tau1,v,Kext,ct,EJv,Kabs,MinDist,dmin,wlam,KabsQHP
	real*8 KabsBBGrains,nRWinteract,vx,vy,vz,inp,vrPress,vzPress
	logical escape,hitstar,locfield,RandomWalk
	integer inext,jnext,ntrace,j,ii,iopac
c	character*500 s

	escape=.false.
	hitstar=.false.
	tau=0d0
	ntrace=0

	if(RNDW.and.C(phot%i,phot%j)%thick.and..not.phot%onEdge) then
2		dmin=MinDist(phot)
		if(RandomWalk(phot,dmin,nRWinteract)) then
c			if(nRWinteract.gt.1d0) then
3				call randomdirection(vx,vy,vz)
				inp=phot%vx*vx+phot%vy*vy+phot%vz*vz
				if(inp.lt.0d0) goto 3
				phot%vx=vx
				phot%vy=vy
				phot%vz=vz
c			endif
			goto 2
		endif
	endif

	if(etrace) call DepositEnergy(phot)

	vrPress=sqrt(phot%vx**2+phot%vy**2)
	if((phot%vx*phot%x+phot%vy*phot%y).lt.0d0) vrPress=-vrPress
	if(phot%z.ge.0d0) then
		vzPress=phot%vz
	else
		vzPress=-phot%vz
	endif

1	continue

	if(C(phot%i,phot%j)%lastphotnr.ne.phot%nr.and..not.etrace) then
		C(phot%i,phot%j)%Ni=C(phot%i,phot%j)%Ni+1
		C(phot%i,phot%j)%lastphotnr=phot%nr
	endif

	Kabs=0d0
	Kext=0d0
	KabsQHP=0d0
	if(.not.scattering) then
	do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			Kabs=Kabs+Grain(ii)%KabsL(iopac)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
			if(Grain(ii)%qhp) KabsQHP=KabsQHP+Grain(ii)%KabsL(iopac)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
		enddo
	enddo
	Kext=Kabs
	else
	do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			Kabs=Kabs+Grain(ii)%KabsL(iopac)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
			Kext=Kext+Grain(ii)%KextL(iopac)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
			if(Grain(ii)%qhp) KabsQHP=KabsQHP+Grain(ii)%KabsL(iopac)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
		enddo
	enddo
	endif
	KabsBBGrains=Kabs-KabsQHP
	if(KabsBBGrains.le.0d0) KabsBBGrains=Kabs

	call Trace2edge(phot,v,inext,jnext)

	tau1=Kext*v*C(phot%i,phot%j)%dens*AU
	if((tau+tau1).gt.tau0) then
		v=v*(tau0-tau)/tau1
		call RadiationPressure(phot,v,vrPress,vzPress,Kext)
		tau1=tau0-tau
		tautot=tautot+tau1
		totaldist=totaldist+v
		if(phot%viscous) coolingtime(icoolingtime)=coolingtime(icoolingtime)+v
c		write(s,*) tautot
c		if(s.eq.' NaN'.or.s.eq.' -NaN') then
c			print*,'dens:     ',C(phot%i,phot%j)%dens
c			print*,'tau,tau0: ',tau,tau0
c			print*,'lam:      ',phot%lam,phot%wl1,phot%wl2
c			print*,'Kext,Kabs:',Kext,Kabs,KabsBBGrains
c			stop
c		endif
		phot%x=phot%x+phot%vx*v
		phot%y=phot%y+phot%vy*v
		phot%z=phot%z+phot%vz*v
		tautot=tautot+(tau0-tau)
		phot%onEdge=.false.
		if(locfield.and..not.etrace) then
			EJv=phot%E*v*AU  !*C(phot%i,phot%j)%dens
			if(computeLRF) then
				C(phot%i,phot%j)%LRF(phot%ilam1)=C(phot%i,phot%j)%LRF(phot%ilam1)+EJv*phot%wl1/dnu(phot%ilam1)
				C(phot%i,phot%j)%LRF(phot%ilam2)=C(phot%i,phot%j)%LRF(phot%ilam2)+EJv*phot%wl2/dnu(phot%ilam2)
				C(phot%i,phot%j)%nLRF(phot%ilam1)=C(phot%i,phot%j)%nLRF(phot%ilam1)+1
				C(phot%i,phot%j)%nLRF(phot%ilam2)=C(phot%i,phot%j)%nLRF(phot%ilam2)+1
				if(use_qhp) then
					do ii=1,ngrains
						if(Grain(ii)%qhp) then
						do iopac=1,Grain(ii)%nopac
							C(phot%i,phot%j)%EJvQHP(Grain(ii)%qhpnr)=C(phot%i,phot%j)%EJvQHP(Grain(ii)%qhpnr)
     &	+C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)*EJv*Grain(ii)%KabsL(iopac)
						enddo
						endif
					enddo
				endif
			endif
			C(phot%i,phot%j)%EJv=C(phot%i,phot%j)%EJv+EJv*KabsBBGrains
c ------------------------------------------------
c give some of the UV photons directly to the gas
c			if(phot%lam.lt.0.3) C(phot%i,phot%j)%Egas=C(phot%i,phot%j)%Egas+EJv*(KabsBBGrains*0.05+KabsQHP*0.5)
c ------------------------------------------------
c ------------------------------------------------
			nEJv=nEJv+EJv*KabsBBGrains
			if(.not.tcontact.or.tdes_iter) then
				do ii=1,ngrains
				do iopac=1,Grain(ii)%nopac
					C(phot%i,phot%j)%EJvP(ii)=C(phot%i,phot%j)%EJvP(ii)
     &	+EJv*C(phot%i,phot%j)%wopac(ii,iopac)*Grain(ii)%KabsL(iopac)
				enddo
     			enddo
				C(phot%i,phot%j)%KextLRF=C(phot%i,phot%j)%KextLRF+EJv*Kext
				C(phot%i,phot%j)%ILRF=C(phot%i,phot%j)%ILRF+EJv
     		endif
		endif
		if(multiwav) then
			do ii=1,ngrains
				do iopac=1,Grain(ii)%nopac
					column(ii,iopac)=column(ii,iopac)+v*C(phot%i,phot%j)%dens
     &					*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)*AU
				enddo
			enddo
		endif
		iRWinter=iRWinter+1
		return
	endif

	iRWinter=0

	call RadiationPressure(phot,v,vrPress,vzPress,Kext)

	phot%x=phot%x+phot%vx*v
	phot%y=phot%y+phot%vy*v
	phot%z=phot%z+phot%vz*v
	if(locfield.and..not.etrace) then
		EJv=phot%E*v*AU		!*C(phot%i,phot%j)%dens
		if(computeLRF) then
			C(phot%i,phot%j)%LRF(phot%ilam1)=C(phot%i,phot%j)%LRF(phot%ilam1)+EJv*phot%wl1/dnu(phot%ilam1)
			C(phot%i,phot%j)%LRF(phot%ilam2)=C(phot%i,phot%j)%LRF(phot%ilam2)+EJv*phot%wl2/dnu(phot%ilam2)
			C(phot%i,phot%j)%nLRF(phot%ilam1)=C(phot%i,phot%j)%nLRF(phot%ilam1)+1
			C(phot%i,phot%j)%nLRF(phot%ilam2)=C(phot%i,phot%j)%nLRF(phot%ilam2)+1
			if(use_qhp) then
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
					do iopac=1,Grain(ii)%nopac
						C(phot%i,phot%j)%EJvQHP(Grain(ii)%qhpnr)=C(phot%i,phot%j)%EJvQHP(Grain(ii)%qhpnr)
     &	+C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)*EJv*Grain(ii)%KabsL(iopac)
					enddo
					endif
				enddo
			endif
		endif
		nEJv=nEJv+EJv*KabsBBGrains
		C(phot%i,phot%j)%EJv2=C(phot%i,phot%j)%EJv2+nEJv**2
		nEJv=0d0
		C(phot%i,phot%j)%EJv=C(phot%i,phot%j)%EJv+EJv*KabsBBGrains
c ------------------------------------------------
c give some of the UV photons directly to the gas
c		if(phot%lam.lt.0.3) C(phot%i,phot%j)%Egas=C(phot%i,phot%j)%Egas+EJv*(KabsBBGrains*0.05+KabsQHP*0.5)
c ------------------------------------------------
c ------------------------------------------------
		if(.not.tcontact.or.tdes_iter) then
			do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				C(phot%i,phot%j)%EJvP(ii)=C(phot%i,phot%j)%EJvP(ii)
     &	+EJv*C(phot%i,phot%j)%wopac(ii,iopac)*Grain(ii)%KabsL(iopac)
			enddo
			enddo
			C(phot%i,phot%j)%KextLRF=C(phot%i,phot%j)%KextLRF+EJv*Kext
			C(phot%i,phot%j)%ILRF=C(phot%i,phot%j)%ILRF+EJv
     	endif
	endif
	if(multiwav) then
		do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				column(ii,iopac)=column(ii,iopac)+v*C(phot%i,phot%j)%dens
     &						*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)*AU
			enddo
		enddo
	endif
	phot%onEdge=.true.
	tau=tau+tau1
	tautot=tautot+tau1
	totaldist=totaldist+v
	if(phot%viscous) coolingtime(icoolingtime)=coolingtime(icoolingtime)+v

	if(inext.ge.D%nR) then
		escape=.true.
		return
	endif
	if(inext.lt.0) then
		escape=.true.
		hitstar=.true.
		return
	endif

	if(jnext.eq.D%nTheta-1..and.phot%j.eq.D%nTheta-1.and.inext.eq.phot%i) vzPress=-vzPress

	if(emptylower) then
		if(jnext.eq.D%nTheta-1..and.phot%j.eq.D%nTheta-1.and.inext.eq.phot%i) then
			escape=.true.
			return
		endif
	endif

	phot%i=inext
	phot%j=jnext

	ntrace=ntrace+1

	if(tau.lt.0d0) then
		escape=.true.
		hitstar=.true.
		return
	endif

	if(ntrace.gt.D%nR*D%nTheta*2) then
		print*,'raar!',phot%i,phot%j,phot%x,phot%y,phot%z,sqrt(phot%vx**2+phot%vy**2+phot%vz**2),v,tau,tau0
		escape=.true.
		hitstar=.true.
		return
	endif

	goto 1

	return
	end

	
!c-----------------------------------------------------------------------
!c This subroutine determines the distance a photon has to travel to
!c the next border of the cell. The output variables are v, the distance
!c inext, the indicator for the next radial cell, and jnext, for the
!c next theta cell.
!c-----------------------------------------------------------------------
	subroutine Trace2edge(phot,v,inext,jnext)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 b,v,r,R1,R2,T1,T2,vR1,vR2,vT1,vT2
	integer inext,jnext
	logical hitR1,hitR2,hitR,hitT1,hitT2,hitT,hitTsame

	r=phot%x**2+phot%y**2+phot%z**2
	R1=C(phot%i,phot%j)%xedge2(1)
	R2=C(phot%i,phot%j)%xedge2(2)
	T1=C(phot%i,phot%j)%xedge2(3)
	T2=C(phot%i,phot%j)%xedge2(4)

	b=2d0*(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz)

	if(.not.phot%onEdge) then
		hitR1=hitR(phot,R1,r,b,vR1)
		hitR2=hitR(phot,R2,r,b,vR2)
		hitT1=hitT(phot,T1,r,b,vT1)
		hitT2=hitT(phot,T2,r,b,vT2)
	else
	if(phot%edgeNr.eq.1) then
		hitR1=.false.
		vR1=1d200
		hitR2=hitR(phot,R2,r,b,vR2)
		hitT1=hitT(phot,T1,r,b,vT1)
		hitT2=hitT(phot,T2,r,b,vT2)
	else if(phot%edgeNr.eq.2) then
		hitR1=hitR(phot,R1,r,b,vR1)
		hitR2=.true.
		vR2=-b
		hitT1=hitT(phot,T1,r,b,vT1)
		hitT2=hitT(phot,T2,r,b,vT2)
	else if(phot%edgeNr.eq.3) then
		hitR1=hitR(phot,R1,r,b,vR1)
		hitR2=hitR(phot,R2,r,b,vR2)
		hitT1=.false.
		vT1=1d200
		if(phot%j.ne.(D%nTheta-1)) then
			hitT2=hitT(phot,T2,r,b,vT2)
		else
			vT2=-phot%z/phot%vz
			if(vT2.gt.0d0) then
				hitT2=.true.
			else
				hitT2=.false.
			endif
		endif
	else if(phot%edgeNr.eq.4) then
		hitR1=hitR(phot,R1,r,b,vR1)
		hitR2=hitR(phot,R2,r,b,vR2)
		if(phot%j.ne.1) then
			hitT1=hitT(phot,T1,r,b,vT1)
		else
			hitT1=.false.
			vT1=1d200
		endif
		if(phot%j.ne.(D%nTheta-1)) then
			hitT2=hitTsame(phot,T2,r,b,vT2)
		else
			hitT2=.false.
			vT2=1d200
		endif
	endif
	endif

	if(.not.hitR2) then
		print*,'Cannot hit outer boundary...'
		print*,sqrt(phot%x**2+phot%y**2+phot%z**2),D%R(phot%i),D%R(phot%i+1)
		print*,sqrt(phot%vx**2+phot%vy**2+phot%vz**2)
	endif

	if(.not.hitR1.and..not.hitR2.and..not.hitT1.and..not.hitT2) then
		print*,'nothing to hit!',vT1,vT2,vR1,vR2
		hitR1=hitR(phot,R1,r,b,vR1)
		hitR2=hitR(phot,R2,r,b,vR2)
		hitT1=hitT(phot,T1,r,b,vT1)
		hitT2=hitT(phot,T2,r,b,vT2)
		if(.not.hitR1.and..not.hitR2.and..not.hitT1.and..not.hitT2) print*,'still nothing'
		print*,phot%z,phot%vz,T1,T2
		stop
	endif

	v=1d200
	if(hitR1.and.vR1.lt.v) then
		v=vR1
		inext=phot%i-1
		jnext=phot%j
		phot%edgeNr=2
	endif
	if(hitR2.and.vR2.lt.v) then
		v=vR2
		inext=phot%i+1
		jnext=phot%j
		phot%edgeNr=1
	endif
	if(hitT1.and.vT1.lt.v) then
		v=vT1
		inext=phot%i
		jnext=phot%j-1
		phot%edgeNr=4
	endif
	if(hitT2.and.vT2.lt.v) then
		v=vT2
		inext=phot%i
		jnext=phot%j+1
		phot%edgeNr=3
	endif


	if(jnext.eq.D%nTheta) then
		jnext=D%nTheta-1
		phot%edgeNr=4
	endif
	if(jnext.eq.0) then
		jnext=1
		phot%edgeNr=3
	endif


	return
	end	
	
	logical function hitR(phot,Rad,r,b,v)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 Rad,r,b,cc,discr,vr1,vr2,v,q
	
	hitR=.false.
	v=1d200

	cc=r-Rad
	discr=(b**2-4d0*cc)
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(b.gt.0d0) then
			q=-0.5d0*(b+discr)
		else
			q=-0.5d0*(b-discr)
		endif
		vr1=q
		vr2=cc/q
		if(vr1.gt.0d0) then
			v=vr1
			hitR=.true.
		endif
		if(vr2.gt.0d0.and.vr2.lt.v) then
			v=vr2
			hitR=.true.
		endif
	endif
	return
	end

	logical function hitT(phot,Thet,r,b,v)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 Thet,r,b,at,bt,ct,discr,vt1,vt2,v,q

	hitT=.false.
	v=1d200

	at=Thet-phot%vz*phot%vz
	bt=Thet*b-2d0*phot%z*phot%vz
	ct=Thet*r-phot%z*phot%z
	discr=bt*bt-4d0*at*ct
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(bt.gt.0d0) then
			q=-0.5d0*(bt+discr)
		else
			q=-0.5d0*(bt-discr)
		endif
		vt1=q/at
		vt2=ct/q
		if(vt1.gt.0d0) then
			v=vt1
			hitT=.true.
		endif
		if(vt2.gt.0d0.and.vt2.lt.v) then
			v=vt2
			hitT=.true.
		endif
	endif
	return
	end

	logical function hitTsame(phot,Thet,r,b,v)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 Thet,r,b,at,bt,v

	hitTsame=.true.
	v=1d200

	bt=Thet*b-2d0*phot%z*phot%vz
	at=Thet-phot%vz**2
	v=-bt/at
	if(v.le.0d0) hitTsame=.false.

	return
	end

!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------

	
	subroutine checkbox(phot)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 theta,r,R1,R2,T1,T2
	
	r=phot%x**2+phot%y**2+phot%z**2
	theta=phot%z**2/r
	R1=C(phot%i,phot%j)%xedge(1)**2*0.999
	R2=C(phot%i,phot%j)%xedge(2)**2*1.001
	T1=C(phot%i,phot%j)%xedge(3)**2*1.001
	T2=C(phot%i,phot%j)%xedge(4)**2*0.999
	if(theta.gt.T1) print*,'oeps! T1',theta,T1,T2,phot%i,phot%j
	if(theta.lt.T2) print*,'oeps! T2',theta,T1,T2,phot%i,phot%j
	if(r.lt.R1) print*,'oeps! R1',r,R1,R2,phot%i,phot%j
	if(r.gt.R2) print*,'oeps! R2',r,R1,R2,phot%i,phot%j
	
	if(theta.gt.T1.or.theta.lt.T2.or.r.lt.R1.or.r.gt.R2) then
		print*,sqrt(r),sqrt(theta)
		print*,C(phot%i,phot%j)%xedge(1)
		print*,C(phot%i,phot%j)%xedge(2)
		print*,C(phot%i,phot%j)%xedge(3)
		print*,C(phot%i,phot%j)%xedge(4)
		stop
	endif
	return
	end


	subroutine allinfo(phot,CC)
	use Parameters
	type(photon) phot
	type(Cell) CC
	
	print*,phot%x,phot%y,phot%z
	print*,sqrt(phot%x**2+phot%y**2+phot%z**2)
	print*,phot%vx,phot%vy,phot%vz
	print*,CC%xedge(1)
	print*,CC%xedge(2)
	print*,CC%xedge(3)
	print*,CC%xedge(4)
	print*,phot%i,phot%j,D%nTheta
	print*,phot%edgenr
	stop
	return
	end
	

	subroutine DepositEnergy(phot0)
	use Parameters
	IMPLICIT NONE
	type(photon) phot0,phot
	real*8 frac,exptau,Kabs,Kext,v,tau,EJv
	integer inext,jnext,ii,ntrace,iopac
	
	phot=phot0

	frac=1d0
	ntrace=0d0
	
1	continue

	if(C(phot%i,phot%j)%lastphotnr.ne.phot%nr) then
		C(phot%i,phot%j)%Ni=C(phot%i,phot%j)%Ni+1
		C(phot%i,phot%j)%lastphotnr=phot%nr
	endif

	Kabs=0d0
	Kext=0d0
	if(.not.scattering) then
	do ii=1,ngrains
	do iopac=1,Grain(ii)%nopac
		Kabs=Kabs+(Grain(ii)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kabs(iopac,phot%ilam2)*phot%wl2)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
	enddo
	enddo
	Kext=Kabs
	else
	do ii=1,ngrains
	do iopac=1,Grain(ii)%nopac
		Kabs=Kabs+(Grain(ii)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kabs(iopac,phot%ilam2)*phot%wl2)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
		Kext=Kext+(Grain(ii)%Kext(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kext(iopac,phot%ilam2)*phot%wl2)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
	enddo
	enddo
	endif

	call Trace2edge(phot,v,inext,jnext)

	tau=Kext*v*C(phot%i,phot%j)%dens*AU

	phot%x=phot%x+phot%vx*v
	phot%y=phot%y+phot%vy*v
	phot%z=phot%z+phot%vz*v

	exptau=exp(-tau)
	if(tau.gt.1d-6) then
		EJv=phot%E*frac*(1d0-exptau)/(C(phot%i,phot%j)%dens*Kext)
	else
		EJv=phot%E*frac*v*AU		!*C(phot%i,phot%j)%dens
	endif

	nEJv=nEJv+EJv*Kabs
	C(phot%i,phot%j)%EJv2=C(phot%i,phot%j)%EJv2+nEJv**2
	nEJv=0d0
	C(phot%i,phot%j)%EJv=C(phot%i,phot%j)%EJv+EJv*Kabs
	if(.not.tcontact.or.tdes_iter) then
		do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			C(phot%i,phot%j)%EJvP(ii)=C(phot%i,phot%j)%EJvP(ii)
     &	+EJv*C(phot%i,phot%j)%wopac(ii,iopac)*(Grain(ii)%Kabs(iopac,phot%ilam1)*phot%wl1
     &	+Grain(ii)%Kabs(iopac,phot%ilam2)*phot%wl2)
		enddo
		enddo
	endif

	phot%onEdge=.true.

	if(inext.ge.D%nR.or.inext.lt.0) return

	phot%i=inext
	phot%j=jnext
	frac=frac*exptau
	if(frac.lt.1d-50) return

	ntrace=ntrace+1
	if(ntrace.gt.D%nR*D%nTheta) then
		print*,'raar!',phot%i,phot%j,phot%x,phot%y,phot%z
		return
	endif

	goto 1
	
	return
	end
	


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine tau2exit(phot,tau)
	use Parameters
	IMPLICIT NONE
	type(Photon) phot
	real*8 tau,v,Kext
	integer inext,jnext,ntrace,i,j,ii,iopac

	tau=0d0

1	continue

	call Trace2edge(phot,v,inext,jnext)
	phot%onEdge=.true.

	phot%x=phot%x+phot%vx*v
	phot%y=phot%y+phot%vy*v
	phot%z=phot%z+phot%vz*v

	Kext=0d0
	do ii=1,ngrains
	do iopac=1,Grain(ii)%nopac
		Kext=Kext+(Grain(ii)%Kext(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kext(iopac,phot%ilam2)*phot%wl2)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
	enddo
	enddo
	tau=tau+v*C(phot%i,phot%j)%dens*AU*Kext

	if(inext.ge.D%nR.or.inext.lt.0) return

	phot%i=inext
	phot%j=jnext

	goto 1

	return
	end

