c		20071027 MM: Fixed the gridding for tracing the observations
c					 when the grid is not finest at the inner edge.
c		20071109 MM: Introduced the weighted emission function to
c					 increase the accuracy of the determination of 
c					 the diffuse field.
c		20071126 MM: Added the (Z)impol output mode which is Q-U

	subroutine TraceFlux(image,lam0,flux,scatflux,fluxQ,Nphot,NphotStar,opening,inc,mask,wmask)
	use Parameters
	use InputOutput
	IMPLICIT NONE
	real*8 lam0,tau,phi,flux,opening,inc,mask,wmask
	integer i,j,k,ii,jj,Nphot,NphotStar
	type(RPhiImage) image
	real*8 tau_e,tau_s,tau_a,w1,w2,nu0,exptau_e
	real*8 wT1,wT2,emis(0:D%nR,D%nTheta),frac,wl1,wl2,Ksca
	integer ilam1,ilam2,iT,ip,jp,kp,jj1,jj2,djj,njj,irg,iopac
	real*8 scat(2,0:D%nR,D%nTheta),fact,scatflux
	real*8 scatQ(2,0:D%nR,D%nTheta),scatU(2,0:D%nR,D%nTheta),ww
	real*8 scatV(2,0:D%nR,D%nTheta),fluxQ,ReadMCScatt,nf,sf,fracirg(D%nTheta,360)
	real*8 x_scat,x_scatQ,x_scatU,x_scatV
	real*8,allocatable :: scatim(:,:),fluxcontr(:,:,:)
	character*500 fluxfile
	integer il10,il100,il1000,il10000,i10
	real*8 i1,il1,Planck,frac_opening
	logical alltrace
	real*8 Tgas,G
	parameter(G=6.67300d-8) ! in cm^3/g/s

	tau_max=1d8

	if(lam0.lt.lam(1).or.lam0.gt.lam(nlam)) then
		write(*,'("Wavelength not in the grid: ",f10.3)') lam0
		return
	endif

	allocate(scatim(image%nr,image%nphi))
	
	if(outfluxcontr) then
		allocate(fluxcontr(0:D%nR,0:D%nTheta,2))
		fluxcontr=0d0
	endif

	alltrace=.true.
	do ii=1,ngrains
		if(.not.Grain(ii)%trace) alltrace=.false.
	enddo

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
	nu0=2.9979d14/lam0

	do i=1,nplanets
		call ReadPlanet(Planets(i),lam0)
	enddo

	if(fastobs) then
		do i=0,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%Kabs=0d0
			C(i,j)%Kext=0d0
			C(i,j)%Ksca=0d0
			do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				C(i,j)%Kabs=C(i,j)%Kabs+(wl1*Grain(ii)%Kabs(iopac,ilam1)
     &			+wl2*Grain(ii)%Kabs(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
				C(i,j)%Kext=C(i,j)%Kext+(wl1*Grain(ii)%Kext(iopac,ilam1)
     &			+wl2*Grain(ii)%Kext(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
				C(i,j)%Ksca=C(i,j)%Ksca+(wl1*Grain(ii)%Ksca(iopac,ilam1)
     &			+wl2*Grain(ii)%Ksca(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
			enddo
			enddo
			C(i,j)%Albedo=C(i,j)%Ksca/C(i,j)%Kext
			do irg=1,C(i,j)%nrg
				C(i,j)%scattfield(irg,0:NPHISCATT/2,1:2)=(wl1*C(i,j)%LRF(ilam1)+wl2*C(i,j)%LRF(ilam2))
     &					*C(i,j)%V/2d0
			enddo
		enddo
		enddo
	else

	if(storescatt) then
		if(.not.scattcomputed(ilam1)) then
			do i=0,D%nR
			do j=1,D%nTheta-1
				C(i,j)%scattfield(1,0,1:2)=0d0
			enddo
			enddo
			call TraceMono(lam(ilam1),Nphot,image%angle,NphotStar)
			do i=0,D%nR
			do j=1,D%nTheta-1
				C(i,j)%scattfield(1,ilam1,1)=C(i,j)%scattfield(1,0,1)
				C(i,j)%scattfield(1,ilam1,2)=C(i,j)%scattfield(1,0,2)
			enddo
			enddo
			scattcomputed(ilam1)=.true.
			nscattcomputed(ilam1)=Nphot+NphotStar
		else if(nscattcomputed(ilam1).lt.Nphot) then
			do i=0,D%nR
			do j=1,D%nTheta-1
				C(i,j)%scattfield(1,0,1:2)=0d0
			enddo
			enddo
			call TraceMono(lam(ilam1),Nphot-nscattcomputed(ilam1),image%angle,NphotStar)
			do i=0,D%nR
			do j=1,D%nTheta-1
				C(i,j)%scattfield(1,ilam1,1:2)=
     &	(C(i,j)%scattfield(1,ilam1,1:2)*real(nscattcomputed(ilam1))+
     &	 C(i,j)%scattfield(1,0,1:2)*real(Nphot-nscattcomputed(ilam1)))/real(Nphot)
			enddo
			enddo
			scattcomputed(ilam1)=.true.
			nscattcomputed(ilam1)=Nphot+NphotStar
		endif
		if(.not.scattcomputed(ilam2)) then
			do i=0,D%nR
			do j=1,D%nTheta-1
				C(i,j)%scattfield(1,0,1:2)=0d0
			enddo
			enddo
			call TraceMono(lam(ilam2),Nphot,image%angle,NphotStar)
			do i=0,D%nR
			do j=0,D%nTheta
				C(i,j)%scattfield(1,ilam2,1:2)=C(i,j)%scattfield(1,0,1:2)
			enddo
			enddo
			scattcomputed(ilam2)=.true.
			nscattcomputed(ilam2)=Nphot+NphotStar
		else if(nscattcomputed(ilam2).lt.Nphot) then
			do i=0,D%nR
			do j=1,D%nTheta-1
				C(i,j)%scattfield(1,0,1:2)=0d0
			enddo
			enddo
			call TraceMono(lam(ilam2),Nphot-nscattcomputed(ilam2),image%angle,NphotStar)
			do i=0,D%nR
			do j=0,D%nTheta
				C(i,j)%scattfield(1,ilam2,1:2)=
     &	(C(i,j)%scattfield(1,ilam2,1:2)*real(nscattcomputed(ilam2))+
     &	 C(i,j)%scattfield(1,0,1:2)*real(Nphot-nscattcomputed(ilam2)))/real(Nphot)
			enddo
			enddo
			scattcomputed(ilam2)=.true.
			nscattcomputed(ilam2)=Nphot+NphotStar
		endif
		do i=0,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%Kabs=0d0
			C(i,j)%Kext=0d0
			C(i,j)%Ksca=0d0
			do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				C(i,j)%Kabs=C(i,j)%Kabs+(wl1*Grain(ii)%Kabs(iopac,ilam1)
     &			+wl2*Grain(ii)%Kabs(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
				C(i,j)%Kext=C(i,j)%Kext+(wl1*Grain(ii)%Kext(iopac,ilam1)
     &			+wl2*Grain(ii)%Kext(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
				C(i,j)%Ksca=C(i,j)%Ksca+(wl1*Grain(ii)%Ksca(iopac,ilam1)
     &			+wl2*Grain(ii)%Ksca(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
			enddo
			enddo
			C(i,j)%Albedo=C(i,j)%Ksca/C(i,j)%Kext
		enddo
		enddo
	else if(scattering.and.(Nphot.ne.0.or.NphotStar.ne.0)) then
		do i=0,D%nR
		do j=1,D%nTheta-1
			C(i,j)%scattfield=0d0
			if(scat_how.eq.2) then
				C(i,j)%scattQ=0d0
				C(i,j)%scattU=0d0
				C(i,j)%scattV=0d0
			endif
		enddo
		enddo
		call TraceMono(lam0,Nphot,image%angle,NphotStar)
	else
		do i=0,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%Kabs=0d0
			C(i,j)%Kext=0d0
			C(i,j)%Ksca=0d0
			do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				C(i,j)%Kabs=C(i,j)%Kabs+(wl1*Grain(ii)%Kabs(iopac,ilam1)
     &			+wl2*Grain(ii)%Kabs(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
				C(i,j)%Kext=C(i,j)%Kext+(wl1*Grain(ii)%Kext(iopac,ilam1)
     &			+wl2*Grain(ii)%Kext(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
				C(i,j)%Ksca=C(i,j)%Ksca+(wl1*Grain(ii)%Ksca(iopac,ilam1)
     &			+wl2*Grain(ii)%Ksca(iopac,ilam2))*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
			enddo
			enddo
			C(i,j)%Albedo=C(i,j)%Ksca/C(i,j)%Kext
		enddo
		enddo
	endif
	endif

!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(j,ii,iT,wT1,wT2,iopac)
!$OMP& SHARED(C,D,BB,Grain,scat,emis,storescatt,wl1,ilam1,wl2,ilam2,traceemis,tracegas,
!$OMP&   gas2dust,usetgas,lam0,tcontact,ngrains)
!$OMP DO
	do i=0,D%nR-1
	do j=1,D%nTheta-1
		scat(1:2,i,j)=0d0
		if(storescatt) then
		scat(1:2,i,j)=(wl1*C(i,j)%scattfield(1,ilam1,1:2)+wl2*C(i,j)%scattfield(1,ilam2,1:2))
     &		/C(i,j)%V
		scat(1:2,i,j)=scat(1:2,i,j)*C(i,j)%Albedo
		endif

		if(traceemis) then
		if(.not.tcontact) then
			emis(i,j)=0d0
			do ii=1,ngrains
				if(Grain(ii)%trace) then
				if(.not.Grain(ii)%qhp) then
				iT=C(i,j)%TP(ii)/dT
				wT1=real(iT+1)-C(i,j)%TP(ii)/dT
				wT2=C(i,j)%TP(ii)/dT-real(iT)
				if(iT.lt.1) iT=1
				if(iT.gt.TMAX-1) iT=TMAX-1
				do iopac=1,Grain(ii)%nopac
					emis(i,j)=emis(i,j)
     &	+(wl1*(wT1*BB(ilam1,iT)+wT2*BB(ilam1,iT+1))*Grain(ii)%Kabs(iopac,ilam1)
     &	+ wl2*(wT1*BB(ilam2,iT)+wT2*BB(ilam2,iT+1))*Grain(ii)%Kabs(iopac,ilam2))
     &	*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)/C(i,j)%Kabs
				enddo
				else
				emis(i,j)=emis(i,j)+C(i,j)%EJvQHP(Grain(ii)%qhpnr)*(wl1*C(i,j)%QHP(Grain(ii)%qhpnr,ilam1)
     &					+wl2*C(i,j)%QHP(Grain(ii)%qhpnr,ilam2))/C(i,j)%Kabs
     			endif
				endif
			enddo
			emis(i,j)=emis(i,j)*(1d0-C(i,j)%Albedo)
		else
			iT=C(i,j)%T/dT
			wT1=real(iT+1)-C(i,j)%T/dT
			wT2=C(i,j)%T/dT-real(iT)
			if(iT.lt.1) iT=1
			if(iT.gt.TMAX-1) iT=TMAX-1
			emis(i,j)=0d0
			do ii=1,ngrains
				if(Grain(ii)%trace) then
				if(.not.Grain(ii)%qhp) then
				do iopac=1,Grain(ii)%nopac
					emis(i,j)=emis(i,j)
     &	+(wl1*(wT1*BB(ilam1,iT)+wT2*BB(ilam1,iT+1))*Grain(ii)%Kabs(iopac,ilam1)
     &	+ wl2*(wT1*BB(ilam2,iT)+wT2*BB(ilam2,iT+1))*Grain(ii)%Kabs(iopac,ilam2))
     &	*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)/C(i,j)%Kabs
				enddo
				else
				emis(i,j)=emis(i,j)+C(i,j)%EJvQHP(Grain(ii)%qhpnr)*(wl1*C(i,j)%QHP(Grain(ii)%qhpnr,ilam1)
     &					+wl2*C(i,j)%QHP(Grain(ii)%qhpnr,ilam2))/C(i,j)%Kabs
				endif
				endif
			enddo
			emis(i,j)=emis(i,j)*(1d0-C(i,j)%Albedo)
		endif
		else
			emis(i,j)=0d0
		endif
		if(useTgas.and.tracegas) then
			iT=C(i,j)%Tgas/dT
			wT1=real(iT+1)-C(i,j)%Tgas/dT
			wT2=C(i,j)%Tgas/dT-real(iT)
			if(iT.lt.1) iT=1
			if(iT.lt.TMAX) then
				emis(i,j)=emis(i,j)
     &	+(wl1*(wT1*BB(ilam1,iT)+wT2*BB(ilam1,iT+1))
     &	+ wl2*(wT1*BB(ilam2,iT)+wT2*BB(ilam2,iT+1)))
     &	*C(i,j)%KappaGas*C(i,j)%gasdens*gas2dust*(1d0-C(i,j)%Albedo)/(C(i,j)%Kabs*C(i,j)%dens)
     		else
				emis(i,j)=emis(i,j)+Planck(C(i,j)%Tgas,lam0)
     &	*C(i,j)%KappaGas*C(i,j)%gasdens*gas2dust*(1d0-C(i,j)%Albedo)/(C(i,j)%Kabs*C(i,j)%dens)
			endif     		
		endif
	enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	write(*,'("Integrating formal solution")')
	write(*,'("Wavelength:",f10.3)') lam0

	write(9,'("Integrating formal solution")')
	write(9,'("Wavelength:",f10.3)') lam0

	do j=1,D%nTheta-1
		do i=1,C(1,j)%nrg
			fracirg(j,i)=(cos(D%thet(j)+real(i-1)*(D%thet(j+1)-D%thet(j))/real(C(1,j)%nrg))-
     &		  cos(D%thet(j)+real(i)*(D%thet(j+1)-D%thet(j))/real(C(1,j)%nrg)))
     &			/(D%Theta(j)-D%Theta(j+1))
		enddo
	enddo
	call tellertje(1,100)
!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,j,k,tau,fact,ip,jp,kp,irg,jj1,jj2,djj,njj,jj,ww,tau_e,Ksca,frac_opening,
!$OMP&    w1,w2,exptau_e,x_scat,x_scatQ,x_scatU,x_scatV,frac,Tgas,iT)
!$OMP& SHARED(image,Nphot,NphotStar,C,scat,scatQ,scatU,scatV,scat_how,scatim,storescatt,
!$OMP&    scattering,fracirg,alltrace,ngrains,Grain,wl1,ilam1,wl2,ilam2,opening,Rinner_gas,
!$OMP&    outfluxcontr,fluxcontr,emis,tau_max,tracestar,D,dimstar,lam0,inner_gas,BB,inc,T_BG,bg_correct)
!$OMP DO
	do i=1,image%nr
!$OMP CRITICAL
	call tellertje(i+1,image%nr+2)
!$OMP END CRITICAL
	do j=1,image%nphi
		image%image(i,j)=0d0
		if(scat_how.eq.2) then
			image%imageQ(i,j)=0d0
			image%imageU(i,j)=0d0
			image%imageV(i,j)=0d0
		endif
		scatim(i,j)=0d0
		tau=0d0
		fact=1d0

		do k=1,image%p(i,j)%n

		ip=image%p(i,j)%i(k)
		jp=image%p(i,j)%j(k)
		kp=image%p(i,j)%k(k)
		irg=image%p(i,j)%irg(k)
		
		x_scat=scat(kp,ip,jp)
		if(scat_how.eq.2) then
			x_scatQ=scatQ(kp,ip,jp)
			x_scatU=scatU(kp,ip,jp)
			x_scatV=scatV(kp,ip,jp)
		endif
		if(.not.storescatt.and.scattering.and.(Nphot+NphotStar).ne.0) then
		jj1=image%p(i,j)%jphi1(k)
		jj2=image%p(i,j)%jphi2(k)
		x_scat=0d0
		if(scat_how.eq.2) then
			x_scatQ=0d0
			x_scatU=0d0
			x_scatV=0d0
		endif
		djj=-1
		if(jj1.le.jj2) djj=1
		njj=0
		jj=jj1
1		ww=C(ip,jp)%Albedo/(C(ip,jp)%V*fracirg(jp,irg))
		x_scat=x_scat+ww*C(ip,jp)%scattfield(irg,jj,kp)
		if(scat_how.eq.2) then
			x_scatQ=x_scatQ+ww*C(ip,jp)%scattQ(irg,jj,kp)
			x_scatU=x_scatU+ww*C(ip,jp)%scattU(irg,jj,kp)
			x_scatV=x_scatV+ww*C(ip,jp)%scattV(irg,jj,kp)
		endif
		njj=njj+1
		if(jj.eq.jj2) goto 2
		jj=jj+djj
		goto 1
2		x_scat=x_scat/real(njj)
		if(scat_how.eq.2) then
			x_scatQ=x_scatQ/real(njj)
			x_scatU=x_scatU/real(njj)
			x_scatV=x_scatV/real(njj)
		endif
		endif

		tau_e=image%p(i,j)%v(k)*C(ip,jp)%dens*C(ip,jp)%Kext*AU

		if(alltrace) then
			Ksca=C(ip,jp)%Ksca
		else
			Ksca=0d0
			do ii=1,ngrains
				if(Grain(ii)%trace) then
					do iopac=1,Grain(ii)%nopac
						Ksca=Ksca+(wl1*Grain(ii)%Ksca(iopac,ilam1)
     &			    +wl2*Grain(ii)%Ksca(iopac,ilam2))*C(ip,jp)%w(ii)*C(ip,jp)%wopac(ii,iopac)
					enddo
				endif
			enddo
			if(scattering.and.(Nphot+NphotStar).ne.0) then
				x_scat=x_scat*Ksca/C(ip,jp)%Ksca
				if(scat_how.eq.2) then
					x_scatQ=x_scatQ*Ksca/C(ip,jp)%Ksca
					x_scatU=x_scatU*Ksca/C(ip,jp)%Ksca
					x_scatV=x_scatV*Ksca/C(ip,jp)%Ksca
				endif
			endif
		endif

		if(real(image%p(i,j)%jphi2(k)).lt.opening.and.real(image%p(i,j)%jphi1(k)).gt.opening) then
			frac_opening=(real(image%p(i,j)%jphi1(k))-opening)
     &				/abs(real(image%p(i,j)%jphi1(k)-image%p(i,j)%jphi2(k)))
			tau_e=tau_e*frac_opening
		else if(real(image%p(i,j)%jphi1(k)).lt.opening.and.real(image%p(i,j)%jphi2(k)).gt.opening) then
			frac_opening=(real(image%p(i,j)%jphi2(k))-opening)
     &				/abs(real(image%p(i,j)%jphi1(k)-image%p(i,j)%jphi2(k)))
			tau_e=tau_e*frac_opening
		else if(real(image%p(i,j)%jphi1(k)).lt.opening.and.real(image%p(i,j)%jphi2(k)).lt.opening) then
			tau_e=0d0
		endif

		if(ip.ne.0) then
		if(outfluxcontr) then
			if(i.lt.image%nr) then
				w1=(image%R(i+1)-image%R(i))*pi*abs(image%R(i))*AU**2/real(image%nPhi)
				if(tau_e.lt.1d-6) then
					fluxcontr(ip,jp,kp)=fluxcontr(ip,jp,kp)+w1*(2d0*x_scat+emis(ip,jp))*tau_e*fact
				else
					exptau_e=exp(-tau_e)
					frac=(1d0-exptau_e)
					fluxcontr(ip,jp,kp)=fluxcontr(ip,jp,kp)+w1*(2d0*x_scat+emis(ip,jp))*frac*fact
				endif
			endif
			if(i.gt.1) then
				w2=(image%R(i)-image%R(i-1))*pi*abs(image%R(i))*AU**2/real(image%nPhi)
				if(tau_e.lt.1d-6) then
					fluxcontr(ip,jp,kp)=fluxcontr(ip,jp,kp)+w2*(2d0*x_scat+emis(ip,jp))*tau_e*fact
				else
					exptau_e=exp(-tau_e)
					frac=(1d0-exptau_e)
					fluxcontr(ip,jp,kp)=fluxcontr(ip,jp,kp)+w2*(2d0*x_scat+emis(ip,jp))*frac*fact
				endif
			endif
		endif

		if(tau_e.lt.1d-6) then
			image%image(i,j)=image%image(i,j)+(2d0*x_scat+emis(ip,jp))*tau_e*fact
			if(scat_how.eq.2) then
				image%imageQ(i,j)=image%imageQ(i,j)+2d0*x_scatQ*tau_e*fact
				image%imageU(i,j)=image%imageU(i,j)+2d0*x_scatU*tau_e*fact
				image%imageV(i,j)=image%imageV(i,j)+2d0*x_scatV*tau_e*fact
			endif
			scatim(i,j)=scatim(i,j)+2d0*x_scat*tau_e*fact

			fact=fact*(1d0-tau_e)
		else
			exptau_e=exp(-tau_e)
			frac=(1d0-exptau_e)
			image%image(i,j)=image%image(i,j)+(2d0*x_scat+emis(ip,jp))*frac*fact
			if(scat_how.eq.2) then
				image%imageQ(i,j)=image%imageQ(i,j)+2d0*x_scatQ*frac*fact
				image%imageU(i,j)=image%imageU(i,j)+2d0*x_scatU*frac*fact
				image%imageV(i,j)=image%imageV(i,j)+2d0*x_scatV*frac*fact
			endif
			scatim(i,j)=scatim(i,j)+2d0*x_scat*frac*fact

			fact=fact*exptau_e
		endif
		tau=tau+tau_e
		else
			if(inner_gas.and.jp.eq.D%nTheta-1.and.irg.eq.1) then
				if(image%p(i,j)%rad(k).gt.(Rinner_gas*D%Rstar/AU)) then
					Tgas=(3d0*G*D%Mstar*D%Mdot*(1d0-sqrt(D%Rstar/(image%p(i,j)%rad(k)*AU)))/(8d0*pi*(image%p(i,j)%rad(k)*AU)**3*sigma))**0.25
					image%image(i,j)=image%image(i,j)+Planck(Tgas,lam0)*fact/(8d0*cos(inc))
				endif
			endif
		endif
		if(tau.gt.tau_max) goto 10
		enddo
		if(image%p(i,j)%hitstar.and.tracestar) then
			image%image(i,j)=image%image(i,j)+(D%Fstar(ilam1)*wl1+D%Fstar(ilam2)*wl2)
     &								*fact*dimstar/(pi*(D%R(0)*AU)**2)
		endif
10		continue
		if(bg_correct) image%image(i,j)=image%image(i,j)+Planck(T_BG,lam0)*(fact-1d0)
	enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(100,100)

	flux=0d0
	scatflux=0d0
	fluxQ=0d0

	do i=1,image%nr-1
	do k=1,image%nPhi
		w1=2d0*pi*abs(image%R(i))*AU**2/real(image%nPhi)
		w2=2d0*pi*abs(image%R(i+1))*AU**2/real(image%nPhi)
		if(mask.lt.1d0.and.(image%R(i+1)/(D%distance/parsec)).lt.wmask) then
			w1=w1*mask
			w2=w2*mask
		endif
		flux=flux+(image%R(i+1)-image%R(i))*
     &		(w1*image%image(i,k)+w2*image%image(i+1,k))/2d0
		scatflux=scatflux+(image%R(i+1)-image%R(i))*
     &		(w1*scatim(i,k)+w2*scatim(i+1,k))/2d0
		if(scat_how.eq.2) then
			fluxQ=fluxQ+(image%R(i+1)-image%R(i))*
     &		(w1*image%imageQ(i,k)+w2*image%imageQ(i+1,k))/2d0
		endif
	enddo
	enddo

	if(readmcscat) then
		sf=ReadMCScatt(image%angle,lam0,nf)
		flux=flux-scatflux
		if((real(Nphot)**0.25).lt.nf) scatflux=sf
		flux=flux+scatflux
	endif

	image%flux=flux

	write(*,'("Total flux:     ",e15.3," Jy")') flux*1e23/D%distance**2
	write(9,'("Total flux:     ",e15.3," Jy")') flux*1e23/D%distance**2


	if(outfluxcontr) then
		i10=(180d0*image%angle/pi)/10d0
		i1=(180d0*image%angle/pi)-10d0*i10
		if(i1.ge.9.95d0) then
			i1=i1-9.95d0
			i10=i10+1
		endif

		il10000=image%lam/10000d0
		il1000=(image%lam-10000d0*il10000)/1000d0
		il100=(image%lam-1000d0*il1000-10000d0*il10000)/100d0
		il10=(image%lam-100d0*il100-1000d0*il1000-10000d0*il10000)/10d0
		il1=image%lam-10d0*il10-100d0*il100-1000d0*il1000-10000d0*il10000

		if(outputfits) then
			write(fluxfile,'(a,"fluxfile","_i",i1,f3.1,"_l",i1,i1,i1,i1,f4.2,".fits")')
     & outdir(1:len_trim(outdir)),i10,i1,il10000,il1000,il100,il10,il1
		else
			write(fluxfile,'(a,"fluxfile","_i",i1,f3.1,"_l",i1,i1,i1,i1,f4.2,".dat")')
     & outdir(1:len_trim(outdir)),i10,i1,il10000,il1000,il100,il10,il1
		endif
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		write(*,'("Writing flux contributions to: ",a)') fluxfile(1:len_trim(fluxfile))
		write(9,'("Writing flux contributions to: ",a)') fluxfile(1:len_trim(fluxfile))

c		call MakeFluxImage(fluxcontr,image%R(image%nr),512,2000,image%angle,image%lam,'fluximage')

		do i=1,D%nR-1
		do j=1,D%nTheta-1
			fluxcontr(i,j,1)=fluxcontr(i,j,1)/C(i,j)%V
			fluxcontr(i,j,2)=fluxcontr(i,j,2)/C(i,j)%V
		enddo
		enddo

		call outputstruct(fluxfile,(/'ARRAY  '/),1,0,fluxcontr(1:D%nR-1,1:D%nTheta-1,1:2),2)

		deallocate(fluxcontr)
	endif

	deallocate(scatim)

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine TraceMono(lam0,Nphot,angle,NphotStar)
	use Parameters
	IMPLICIT NONE
	character*500 input,tmp,specfile
	integer i,j,k,l,t1,t2,iangle,nangle,Nphot,iphot,ii,NphotStar
	real*8 ran2,tau,R0,R1,angle,r,distance,tottime,ct,lam0,E
	real*8 th,ph,inp,determineT,fstop,fact,s1,s2
	integer starttime,stoptime,starttrace,cr,ia
	logical escape,hitstar,hitmid,ignore
	type(Photon) phot,phot2,photinit
	integer ilam,nabs,iopac,idumstart
	real*8 x,y,z,phi,theta,Emin,rho,dangle,EnergyTot2
	real*8 EnergyTot,Estar,Rad,VETot,tot,tot2,thet,Eirf
	real*8,allocatable :: VisEmisDis(:,:),EmisDis(:,:),vismass(:,:)
	integer,allocatable :: NEmisDis(:,:)
	real*8 Einner,fact_IRF
	integer Nstar,ip,np,jj,Nmin,iscat,omp_get_thread_num
	type(Mueller) M,ML(ngrains,ngrains2)
	real*8 KabsL(ngrains,ngrains2),KscaL(ngrains,ngrains2),KextL(ngrains,ngrains2)
	type(Cell), pointer :: CC
	type(Particle), pointer :: GG
	
	allocate(VisEmisDis(0:D%nR+1,0:D%nTheta+1))
	allocate(EmisDis(0:D%nR+1,0:D%nTheta+1))
	allocate(vismass(0:D%nR+1,0:D%nTheta+1))
	allocate(NEmisDis(0:D%nR+1,0:D%nTheta+1))
	
	nemit=0

	tautot=0d0

	if(Nphot.eq.0.and.NphotStar.eq.0) return

	write(*,'("Creating local field")')
	write(*,'("Wavelength:",f10.3)') lam0

	write(9,'("Creating local field")')
	write(9,'("Wavelength:",f10.3)') lam0

	if(lam0.le.lam(1)) then
		phot%ilam1=1
		phot%ilam2=2
		phot%wl1=1d0
		phot%wl2=0d0
		phot%lam=lam0
		phot%nu=1d0/lam0
	endif
	do ilam=1,nlam-1
		if(lam0.ge.lam(ilam).and.lam0.le.lam(ilam+1)) then
			phot%ilam1=ilam
			phot%ilam2=ilam+1
			phot%wl1=(lam(ilam+1)-lam0)/(lam(ilam+1)-lam(ilam))
			phot%wl2=(lam0-lam(ilam))/(lam(ilam+1)-lam(ilam))
			phot%lam=lam0
			phot%nu=1d0/lam0
		endif
	enddo
	if(lam0.ge.lam(nlam)) then
		phot%ilam1=nlam-1
		phot%ilam2=nlam
		phot%wl1=0d0
		phot%wl2=1d0
		phot%lam=lam0
		phot%nu=1d0/lam0
	endif

!$OMP PARALLEL IF(ngrains.gt.(omp_get_thread_num()*3))
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ii,iopac,ia,GG)
!$OMP& SHARED(ngrains,Grain,KabsL,KextL,KscaL,phot,ML)
!$OMP DO
	do ii=1,ngrains
		GG => Grain(ii)
		do iopac=1,GG%nopac
			KabsL(ii,iopac)=phot%wl1*GG%Kabs(iopac,phot%ilam1)+phot%wl2*GG%Kabs(iopac,phot%ilam2)
			KextL(ii,iopac)=phot%wl1*GG%Kext(iopac,phot%ilam1)+phot%wl2*GG%Kext(iopac,phot%ilam2)
			KscaL(ii,iopac)=phot%wl1*GG%Ksca(iopac,phot%ilam1)+phot%wl2*GG%Ksca(iopac,phot%ilam2)

			ML(ii,iopac)%IF11=(GG%F(iopac,phot%ilam1)%IF11*GG%Ksca(iopac,phot%ilam1))*phot%wl1
     &			+(GG%F(iopac,phot%ilam2)%IF11*GG%Ksca(iopac,phot%ilam2))*phot%wl2
			ML(ii,iopac)%IF12=(GG%F(iopac,phot%ilam1)%IF12*GG%Ksca(iopac,phot%ilam1))*phot%wl1
     &			+(GG%F(iopac,phot%ilam2)%IF12*GG%Ksca(iopac,phot%ilam2))*phot%wl2
			do ia=1,180
				ML(ii,iopac)%F11(ia)=(GG%F(iopac,phot%ilam1)%F11(ia)*GG%Ksca(iopac,phot%ilam1))*phot%wl1
     &			+(GG%F(iopac,phot%ilam2)%F11(ia)*GG%Ksca(iopac,phot%ilam2))*phot%wl2
				ML(ii,iopac)%F12(ia)=(GG%F(iopac,phot%ilam1)%F12(ia)*GG%Ksca(iopac,phot%ilam1))*phot%wl1
     &			+(GG%F(iopac,phot%ilam2)%F12(ia)*GG%Ksca(iopac,phot%ilam2))*phot%wl2
				ML(ii,iopac)%F22(ia)=(GG%F(iopac,phot%ilam1)%F22(ia)*GG%Ksca(iopac,phot%ilam1))*phot%wl1
     &			+(GG%F(iopac,phot%ilam2)%F22(ia)*GG%Ksca(iopac,phot%ilam2))*phot%wl2
				ML(ii,iopac)%F33(ia)=(GG%F(iopac,phot%ilam1)%F33(ia)*GG%Ksca(iopac,phot%ilam1))*phot%wl1
     &			+(GG%F(iopac,phot%ilam2)%F33(ia)*GG%Ksca(iopac,phot%ilam2))*phot%wl2
				ML(ii,iopac)%F34(ia)=(GG%F(iopac,phot%ilam1)%F34(ia)*GG%Ksca(iopac,phot%ilam1))*phot%wl1
     &			+(GG%F(iopac,phot%ilam2)%F34(ia)*GG%Ksca(iopac,phot%ilam2))*phot%wl2
				ML(ii,iopac)%F44(ia)=(GG%F(iopac,phot%ilam1)%F44(ia)*GG%Ksca(iopac,phot%ilam1))*phot%wl1
     &			+(GG%F(iopac,phot%ilam2)%F44(ia)*GG%Ksca(iopac,phot%ilam2))*phot%wl2
			enddo
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	

!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(j,ii,iopac,ia,CC)
!$OMP& SHARED(D,C,ngrains,Grain,phot,scattering,scat_how,useobspol,KabsL,KextL,KscaL,ML)
!$OMP DO
	do i=0,D%nR-1
	do j=1,D%nTheta-1
		CC => C(i,j)
		CC%Kabs=0d0
		CC%Kext=0d0
		CC%Ksca=0d0
		do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			CC%Kabs=CC%Kabs+KabsL(ii,iopac)*CC%w(ii)*CC%wopac(ii,iopac)
			CC%Kext=CC%Kext+KextL(ii,iopac)*CC%w(ii)*CC%wopac(ii,iopac)
			CC%Ksca=CC%Ksca+KscaL(ii,iopac)*CC%w(ii)*CC%wopac(ii,iopac)
		enddo
		enddo
		CC%Albedo=CC%Ksca/CC%Kext
		if(scattering.and.scat_how.eq.2) then
		if(CC%Ksca.ne.0d0) then
		CC%F%IF11=0d0
		CC%F%IF12=0d0
		do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			CC%F%IF11=CC%F%IF11+ML(ii,iopac)%IF11*CC%w(ii)*CC%wopac(ii,iopac)
			CC%F%IF12=CC%F%IF12+ML(ii,iopac)%IF12*CC%w(ii)*CC%wopac(ii,iopac)
    	enddo
		enddo
		CC%F%IF11=CC%F%IF11/CC%Ksca
		CC%F%IF12=CC%F%IF12/CC%Ksca
		do ia=1,180
			CC%F%F11(ia)=0d0
			CC%F%F12(ia)=0d0
			CC%F%F22(ia)=0d0
			CC%F%F33(ia)=0d0
			CC%F%F34(ia)=0d0
			CC%F%F44(ia)=0d0
			do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				CC%F%F11(ia)=CC%F%F11(ia)+ML(ii,iopac)%F11(ia)*CC%w(ii)*CC%wopac(ii,iopac)
				CC%F%F12(ia)=CC%F%F12(ia)+ML(ii,iopac)%F12(ia)*CC%w(ii)*CC%wopac(ii,iopac)
     			if(useobspol) then
					CC%F%F22(ia)=CC%F%F22(ia)+ML(ii,iopac)%F22(ia)*CC%w(ii)*CC%wopac(ii,iopac)
					CC%F%F33(ia)=CC%F%F33(ia)+ML(ii,iopac)%F33(ia)*CC%w(ii)*CC%wopac(ii,iopac)
					CC%F%F34(ia)=CC%F%F34(ia)+ML(ii,iopac)%F34(ia)*CC%w(ii)*CC%wopac(ii,iopac)
					CC%F%F44(ia)=CC%F%F44(ia)+ML(ii,iopac)%F44(ia)*CC%w(ii)*CC%wopac(ii,iopac)
				endif
			enddo
			enddo
			CC%F%F11(ia)=CC%F%F11(ia)/CC%Ksca
			CC%F%F12(ia)=CC%F%F12(ia)/CC%Ksca
			if(useobspol) then
				CC%F%F22(ia)=CC%F%F22(ia)/CC%Ksca
				CC%F%F33(ia)=CC%F%F33(ia)/CC%Ksca
				CC%F%F34(ia)=CC%F%F34(ia)/CC%Ksca
				CC%F%F44(ia)=CC%F%F44(ia)/CC%Ksca
			endif
		enddo
		else
			CC%F%F11(1:180)=1d0
		endif
		endif
	enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	do i=1,nplanets
		Planets(i)%E=0d0
		Planets(i)%Q=0d0
		Planets(i)%U=0d0
		Planets(i)%V=0d0

		do k=1,D%nR-1
			if(Planets(i)%d.gt.D%R(k).and.Planets(i)%d.le.D%R(k+1)) goto 500
		enddo
500		continue
		Planets(i)%i=k
		do k=1,D%nTheta-1
			if(Planets(i)%theta.gt.D%thet(k).and.Planets(i)%theta.le.D%thet(k+1)) then
				Planets(i)%j=k
				Planets(i)%k=1
			endif
			if((pi-Planets(i)%theta).gt.D%thet(k).and.(pi-Planets(i)%theta).le.D%thet(k+1)) then
				Planets(i)%j=k
				Planets(i)%k=1
			endif
		enddo
		print*,Planets(i)%j
	enddo

	call EmissionDistribution(phot,EmisDis,EnergyTot,EnergyTot2,Estar,Eirf,Einner,vismass)

	call MakeStarScatter(Estar,NphotStar)

	call determine_fact_IRF(fact_IRF,Estar,Eirf,angle)

	write(*,'("Emitting   ",i10," photon packages")') Nphot
	write(9,'("Emitting   ",i10," photon packages")') Nphot

c Start tracing the photons

	photinit=phot
	call tellertje(1,100)
	idumstart=idum
!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(phot,x,y,z,r,ignore,tautot,tau,hitstar,escape,fstop,fact,xsn,ysn,zsn,
!$OMP&   s1,s2,phot2,ninteract,iscat)
!$OMP& SHARED(scat_how,C,EmisDis,EnergyTot,EnergyTot2,Estar,Eirf,Einner,vismass,
!$OMP&   xsf,ysf,zsf,Nphot,forcefirst,photinit,fact_IRF,idumstart)
!$OMP DO
	do iphot=1,Nphot
!$OMP CRITICAL
	call tellertje(iphot+1,Nphot+2)
!$OMP END CRITICAL
	phot=photinit

	idum=idumstart+iphot

	phot%nr=iphot
	call EmitPosition(phot,EmisDis,EnergyTot,EnergyTot2,Estar,Eirf,Einner,vismass,fact_IRF)
	phot%E=phot%E/real(Nphot)
	if(scat_how.eq.2) then
		phot%pol=.false.
		phot%Q=0d0
		phot%U=0d0
		phot%V=0d0
		if(phot%vz.ne.1d0) then
			x=-phot%vy
			y=phot%vx
			z=0d0
		else
			x=1d0
			y=0d0
			z=0d0
		endif
		r=sqrt(x**2+y**2+z**2)
		phot%Sx=x/r
		phot%Sy=y/r
		phot%Sz=z/r
	endif

	tautot=0d0
	ninteract=0

1	continue

	if(.not.forcefirst.or.phot%scatt) then
		tau=-log(ran2(idum))
	else
		phot2=phot
		call trace2exit(phot2,tau,.true.)
		if(tau.lt.1d-9.or.tau.gt.15d0) then
			s1=0d0
			s2=1d0
		else 
			if(tau.gt.1d-6) then
				s1=exp(-tau)
				s2=1d0-s1
			else
				s1=1d0-tau
				s2=tau
			endif
		endif
		tau=abs(-log(ran2(idum)*s2+s1))
		phot%E=phot%E*s2
	endif
		
	call trace2dmono(phot,tau,escape,hitstar)
	if(hitstar) goto 3
	if(escape) goto 3

c	fstop=((C(phot%i,phot%j)%Albedo+1d0)/2d0)
	fstop=C(phot%i,phot%j)%Albedo**0.25d0	! fstop is the chance the photon goes trhough
c	fstop=C(phot%i,phot%j)%Albedo
	fact=C(phot%i,phot%j)%Albedo/fstop

	phot%E=phot%E*fact
	if(scat_how.eq.2) then
		phot%Q=phot%Q*fact
		phot%U=phot%U*fact
		phot%V=phot%V*fact
	endif
	if(ran2(idum).gt.fstop) goto 3

	ninteract=ninteract+1
	phot%scatt=.true.
	if(scat_how.eq.1) then
		call randomdirection(phot%vx,phot%vy,phot%vz)
	else
		call scatangle(phot,C(phot%i,phot%j)%F,iscat)
	endif
	
	goto 1

3	continue
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(100,100)


	if(.not.storescatt.and.scat_how.eq.1) then
	do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%scattfield(1,1:NPHISCATT/2,1)=C(i,j)%scattfield(1,0,1)
			C(i,j)%scattfield(1,1:NPHISCATT/2,2)=C(i,j)%scattfield(1,0,2)
		enddo
	enddo
	endif

	deallocate(VisEmisDis)
	deallocate(EmisDis)
	deallocate(vismass)
	deallocate(NEmisDis)

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine determine_fact_IRF(fact_IRF,Estar,Eirf,angle)
	use Parameters
	IMPLICIT NONE
	real*8 fact_IRF,Estar,Eirf,angle,tau,f
	integer i,j

	if(Eirf.eq.0d0.or..not.use_IRF) then
		fact_IRF=1d0
		return
	endif

	do j=1,D%nTheta-1
		if(cos(angle).lt.D%Theta(j).and.cos(angle).ge.D%Theta(j+1)) exit
	enddo
	tau=0d0
	do i=1,D%nR-1
		tau=tau+(D%R(i+1)-D%R(i))*AU*C(i,j)%dens*C(i,j)%Kext
	enddo
	f=exp(-tau)

	fact_IRF=1d0+(1d0-f)*(Estar/Eirf-1d0)

	if(fact_IRF.lt.(0.1d0*(Estar+Eirf)/Eirf)) fact_IRF=0.1d0*(Estar+Eirf)/Eirf

	if(fact_IRF.lt.10d0) fact_IRF=10d0
	if(fact_IRF.gt.1d4) fact_IRF=1d4

	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine EmissionDistribution(phot,EmisDis,EnergyTot,EnergyTot2,Estar,Eirf,Einner,vismass)
	use Parameters
	IMPLICIT NONE
	real*8 lam0
	real*8 EmisDis(0:D%nR+1,0:D%nTheta+1),EnergyTot,Estar,EnergyTot2
	real*8 Planck,vismass(0:D%nR+1,0:D%nTheta+1),Eirf,Einner
	real*8 tau,taumin,ran2,Rad,Theta,phi
	integer i,j,ii,k,iopac
	type(Photon) phot,phot2

	integer iT
	real*8 wT1,wT2
	
c	if(.not.C(1,1)%exittau) then
c		write(*,'("Finding shortest route to observer")')
c		write(9,'("Finding shortest route to observer")')
c	endif
c	Estar=pi*Planck(D%Tstar,phot%lam)*D%Rstar**2
	Estar=D%Fstar(phot%ilam1)*phot%wl1+D%Fstar(phot%ilam2)*phot%wl2
	if(inner_gas) then
		call InitInnerGasDisk(phot%lam,Einner)
	else
		Einner=0d0
	endif
	if(use_IRF) then
		Eirf=IRF(phot%ilam1)*phot%wl1+IRF(phot%ilam2)*phot%wl2
	else
		Eirf=0d0
	endif
	EnergyTot=0d0
	EnergyTot2=0d0
	if(nexits.ne.0) call tellertje(1,100)
!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(j,ii,iopac,iT,wT1,wT2,phot2,Rad,Theta,phi,tau,k,taumin)
!$OMP& SHARED(nexits,D,C,Grain,EmisDis,EnergyTot,EnergyTot2,Tcontact,ngrains,
!$OMP&    phot,useTgas,tracegas,gas2dust,vismass,BB)

!$OMP DO
	do i=1,D%nR-1
!$OMP CRITICAL
		if(nexits.ne.0) call tellertje(i+1,D%nR+1)
!$OMP END CRITICAL
		do j=1,D%nTheta-1
			if(C(i,j)%T.ne.0d0) then
			if(.not.tcontact) then
				EmisDis(i,j)=0d0
				do ii=1,ngrains
					if(.not.Grain(ii)%qhp) then
						iT=C(i,j)%TP(ii)/dT
						wT1=real(iT+1)-C(i,j)%TP(ii)/dT
						wT2=C(i,j)%TP(ii)/dT-real(iT)
						if(iT.lt.1) iT=1
						if(iT.gt.TMAX-1) iT=TMAX-1
						do iopac=1,Grain(ii)%nopac
							EmisDis(i,j)=EmisDis(i,j)
     &	+(phot%wl1*(wT1*BB(phot%ilam1,iT)+wT2*BB(phot%ilam1,iT+1))*Grain(ii)%Kabs(iopac,phot%ilam1)
     &	+ phot%wl2*(wT1*BB(phot%ilam2,iT)+wT2*BB(phot%ilam2,iT+1))*Grain(ii)%Kabs(iopac,phot%ilam2))
     &	*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
						enddo
					else
						EmisDis(i,j)=EmisDis(i,j)+C(i,j)%EJvQHP(Grain(ii)%qhpnr)*(phot%wl1*C(i,j)%QHP(Grain(ii)%qhpnr,phot%ilam1)
     &					+phot%wl2*C(i,j)%QHP(Grain(ii)%qhpnr,phot%ilam2))
     				endif
				enddo
			else
				iT=C(i,j)%T/dT
				wT1=real(iT+1)-C(i,j)%T/dT
				wT2=C(i,j)%T/dT-real(iT)
				if(iT.lt.1) iT=1
				if(iT.gt.TMAX-1) iT=TMAX-1
				EmisDis(i,j)=0d0
				do ii=1,ngrains
					if(.not.Grain(ii)%qhp) then
						do iopac=1,Grain(ii)%nopac
							EmisDis(i,j)=EmisDis(i,j)
     &	+(phot%wl1*(wT1*BB(phot%ilam1,iT)+wT2*BB(phot%ilam1,iT+1))*Grain(ii)%Kabs(iopac,phot%ilam1)
     &	+ phot%wl2*(wT1*BB(phot%ilam2,iT)+wT2*BB(phot%ilam2,iT+1))*Grain(ii)%Kabs(iopac,phot%ilam2))
     &	*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
						enddo
					else
						EmisDis(i,j)=EmisDis(i,j)+C(i,j)%EJvQHP(Grain(ii)%qhpnr)*(phot%wl1*C(i,j)%QHP(Grain(ii)%qhpnr,phot%ilam1)
     &					+phot%wl2*C(i,j)%QHP(Grain(ii)%qhpnr,phot%ilam2))
					endif
				enddo
			endif

			else
				EmisDis(i,j)=0d0
			endif
			EmisDis(i,j)=EmisDis(i,j)*C(i,j)%dens*C(i,j)%V
			if(useTgas.and.tracegas) then
				iT=C(i,j)%Tgas/dT
				wT1=real(iT+1)-C(i,j)%Tgas/dT
				wT2=C(i,j)%Tgas/dT-real(iT)
				if(iT.lt.1) iT=1
				if(iT.lt.TMAX) then
					EmisDis(i,j)=EmisDis(i,j)
     &	+(phot%wl1*(wT1*BB(phot%ilam1,iT)+wT2*BB(phot%ilam1,iT+1))
     &	+ phot%wl2*(wT1*BB(phot%ilam2,iT)+wT2*BB(phot%ilam2,iT+1)))*C(i,j)%KappaGas*gas2dust
     &	*C(i,j)%gasdens*C(i,j)%V
				else
					EmisDis(i,j)=EmisDis(i,j)
     &	+Planck(C(i,j)%Tgas,phot%lam)*C(i,j)%KappaGas*gas2dust*C(i,j)%gasdens*C(i,j)%V
				endif
			endif

c eliminating 'dark-zone'
			if(nexits.eq.0) then
				C(i,j)%tauexit=0d0
			else
				phot2=phot
				do k=1,abs(nexits)
					Rad=0.5d0
					Rad=sqrt(D%R(i)**2*Rad+D%R(i+1)**2*(1d0-Rad))
					Theta=0.5d0
					Theta=D%Theta(j)*Theta+D%Theta(j+1)*(1d0-Theta)
					phot2%z=Rad*Theta
					phi=1d-4
					phot2%x=Rad*sqrt(1d0-Theta**2)*sin(phi)
					phot2%y=Rad*sqrt(1d0-Theta**2)*cos(phi)
					phi=pi*real(k-1)/real(nexits-1)
					phot2%vx=0d0
					phot2%vy=cos(phi)
					phot2%vz=sin(phi)
					phot2%i=i
					phot2%j=j
					phot2%onEdge=.false.
					call trace2exit(phot2,tau,.false.)
					if(tau.lt.taumin.or.k.eq.1) taumin=tau
				enddo
				C(i,j)%tauexit=taumin
			endif

			if(nexits.ge.0) then
				vismass(i,j)=dexp(-C(i,j)%tauexit)
			else
				vismass(i,j)=(1d0-dexp(-C(i,j)%tauexit))
			endif

!$OMP FLUSH(EnergyTot)
			EnergyTot=EnergyTot+EmisDis(i,j)
			EmisDis(i,j)=EmisDis(i,j)*vismass(i,j)
!$OMP FLUSH(EnergyTot2)
			EnergyTot2=EnergyTot2+EmisDis(i,j)
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	if(nexits.ne.0) call tellertje(100,100)

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine emitposition(phot,EmisDis,EnergyTot,EnergyTot2,Estar,Eirf,Einner,vismass,fact_IRF)
	use Parameters
	IMPLICIT NONE
	real*8 EmisDis(0:D%nR+1,0:D%nTheta+1),EnergyTot,Estar,Er,Et,thet,Eirf,fact_IRF
	real*8 inp,ct,r,Rad,Theta,ran2,phi,EnergyTot2,vismass(0:D%nR+1,0:D%nTheta+1)
	real*8 Einner,Etot,Er2,tot,RadInnerGasDisk,R1,R2
	integer i,j,iter
	type(photon) phot
	logical ignore
	
	phot%pol=.false.

	Etot=(EnergyTot2+Estar+Eirf*fact_IRF+Einner)
	
	Er=Etot*ran2(idum)
	if(Er.lt.Estar) then
		phot%E=Etot
		phot%scatt=.false.
		call randomdirection(phot%x,phot%y,phot%z)
		phot%x=D%R(0)*phot%x
		phot%y=D%R(0)*phot%y
		phot%z=D%R(0)*phot%z
		call randomdirection(phot%vx,phot%vy,phot%vz)
		inp=(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz)
		if(inp.lt.0d0) then
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
				if(C(phot%i,phot%j)%nrg.gt.1) then
					thet=acos(abs(phot%z)/sqrt(phot%x**2+phot%y**2+phot%z**2))
					phot%irg=int(real(C(phot%i,phot%j)%nrg)*(thet-D%thet(phot%j))/(D%thet(phot%j+1)-D%thet(phot%j)))+1
					if(phot%irg.gt.C(phot%i,phot%j)%nrg) phot%irg=C(phot%i,phot%j)%nrg
					if(phot%irg.lt.1) phot%irg=1
				else
					phot%irg=1
				endif
				return
			endif
		enddo
	else if(Er.lt.(Estar+Eirf*fact_IRF)) then
		phot%E=Etot/fact_IRF
		phot%scatt=.true.
		call randomdirection(phot%x,phot%y,phot%z)
		phot%x=D%R(D%nR)*phot%x
		phot%y=D%R(D%nR)*phot%y
		phot%z=D%R(D%nR)*phot%z
7		call randomdirection(phot%vx,phot%vy,phot%vz)
		inp=(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz)
		if(inp.gt.0d0) then
			phot%vx=-phot%vx
			phot%vy=-phot%vy
			phot%vz=-phot%vz
		endif
		if(abs(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz)/D%R(D%nR).lt.ran2(idum)) goto 7

		phot%edgeNr=2
		phot%onEdge=.true.
		ct=abs(phot%z)/D%R(D%nR)
		do j=1,D%nTheta-1
			if(ct.lt.D%Theta(j).and.ct.ge.D%Theta(j+1)) then
				phot%i=D%nR-1
				phot%j=j
				if(C(phot%i,phot%j)%nrg.gt.1) then
					thet=acos(abs(phot%z)/sqrt(phot%x**2+phot%y**2+phot%z**2))
					phot%irg=int(real(C(phot%i,phot%j)%nrg)*(thet-D%thet(phot%j))/(D%thet(phot%j+1)-D%thet(phot%j)))+1
					if(phot%irg.gt.C(phot%i,phot%j)%nrg) phot%irg=C(phot%i,phot%j)%nrg
					if(phot%irg.lt.1) phot%irg=1
				else
					phot%irg=1
				endif
				return
			endif
		enddo
	else if(Er.lt.(Estar+Eirf*fact_IRF+Einner)) then
		phot%E=Etot
		phot%scatt=.true.

		Er2=ran2(idum)*Einner
		Rad=RadInnerGasDisk(Er2)

		j=D%nTheta-1
		Theta=0.5
		Theta=D%Theta(j)*Theta+D%Theta(j+1)*(1d0-Theta)
		phot%z=Rad*Theta
		phi=ran2(idum)*pi*2d0
		phot%x=Rad*sqrt(1d0-Theta**2)*sin(phi)
		phot%y=Rad*sqrt(1d0-Theta**2)*cos(phi)
		call randomdirection(phot%vx,phot%vy,phot%vz)
		phot%i=0
		phot%j=j
		phot%onEdge=.false.

		call randomdirection(phot%vx,phot%vy,phot%vz)

		phot%z=abs(phot%z)
		if(phot%vz.lt.0d0) phot%z=-phot%z

		ct=abs(phot%z)/Rad
		do j=1,D%nTheta-1
			if(ct.lt.D%Theta(j).and.ct.ge.D%Theta(j+1)) then
				phot%j=j
				if(C(phot%i,phot%j)%nrg.gt.1) then
					thet=acos(abs(phot%z)/sqrt(phot%x**2+phot%y**2+phot%z**2))
					phot%irg=int(real(C(phot%i,phot%j)%nrg)*(thet-D%thet(phot%j))/(D%thet(phot%j+1)-D%thet(phot%j)))+1
					if(phot%irg.gt.C(phot%i,phot%j)%nrg) phot%irg=C(phot%i,phot%j)%nrg
					if(phot%irg.lt.1) phot%irg=1
				else
					phot%irg=1
				endif
				return
			endif
		enddo
	else
		phot%scatt=.true.
		Er=Er-Estar-Eirf*fact_IRF-Einner
		Et=0d0
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			Et=Et+EmisDis(i,j)
			if(Er.lt.Et) then
				Rad=ran2(idum)
				Rad=sqrt(D%R(i)**2*Rad+D%R(i+1)**2*(1d0-Rad))
				Theta=ran2(idum)
				Theta=D%Theta(j)*Theta+D%Theta(j+1)*(1d0-Theta)
				phot%z=Rad*Theta
				if(ran2(idum).lt.0.5) phot%z=-phot%z
				phi=ran2(idum)*pi*2d0
				phot%x=Rad*sqrt(1d0-Theta**2)*sin(phi)
				phot%y=Rad*sqrt(1d0-Theta**2)*cos(phi)
				call randomdirection(phot%vx,phot%vy,phot%vz)
				phot%i=i
				phot%j=j
				if(C(phot%i,phot%j)%nrg.gt.1) then
					thet=acos(abs(phot%z)/sqrt(phot%x**2+phot%y**2+phot%z**2))
					phot%irg=int(real(C(phot%i,phot%j)%nrg)*(thet-D%thet(phot%j))/(D%thet(phot%j+1)-D%thet(phot%j)))+1
					if(phot%irg.gt.C(phot%i,phot%j)%nrg) phot%irg=C(phot%i,phot%j)%nrg
					if(phot%irg.lt.1) phot%irg=1
				else
					phot%irg=1
				endif
				phot%onEdge=.false.
				phot%E=Etot/vismass(i,j)
				return
			endif
		enddo
		enddo
	endif
	
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine trace2dmono(phot,tau0,escape,hitstar)
	use Parameters
	IMPLICIT NONE
	type(Photon) phot,phot1(NPHISCATT)
	real*8 tau0,tau,a,b,cc,tau1,det,v
	real*8 at,bt,ct,EJv,r,sin2t(NPHISCATT),cos2t(NPHISCATT)
	logical escape
	logical hitstar,hitmid
	integer inext,jnext,ntrace,i,iangle(NPHISCATT),irgnext
	integer ncells_visit,maxcells_visit
	
	maxcells_visit=(D%nR+D%nTheta+90)*2

	escape=.false.
	hitstar=.false.
	hitmid=.false.
	tau=0d0

	if(scat_how.ne.1.and..not.makeangledependence.and.phot%scatt) then
	do i=1,NPHISCATT
		iangle(i)=180d0*acos(phot%vx*xsf(i)+phot%vy*ysf(i)
     &						+phot%vz*zsf(i))/pi+1
		xsn(i)=phot%vy*zsf(i)-phot%vz*ysf(i)
		ysn(i)=phot%vz*xsf(i)-phot%vx*zsf(i)
		zsn(i)=phot%vx*ysf(i)-phot%vy*xsf(i)
		r=sqrt(xsn(i)**2+ysn(i)**2+zsn(i)**2)
		xsn(i)=xsn(i)/r
		ysn(i)=ysn(i)/r
		zsn(i)=zsn(i)/r
		phot1(i)=phot
		if(phot%pol) then
			call RotateStokes(phot1(i),xsn(i),ysn(i),zsn(i))
		endif
		phot1(i)%Sx=xsn(i)
		phot1(i)%Sy=ysn(i)
		phot1(i)%Sz=zsn(i)
		phot1(i)%vx=xsf(i)
		phot1(i)%vy=ysf(i)
		phot1(i)%vz=zsf(i)
		call MakeRotateStokes(phot1(i),xin(i),yin(i),zin(i),sin2t(i),cos2t(i))
		phot1(i)%vx=phot%vx
		phot1(i)%vy=phot%vy
		phot1(i)%vz=phot%vz
   	enddo
	endif

	ncells_visit=0
1	continue
	ncells_visit=ncells_visit+1
	if(ncells_visit.gt.maxcells_visit) then
		print*,'something went wrong here...'
		escape=.true.
		return
	endif
	C(phot%i,phot%j)%Ni=C(phot%i,phot%j)%Ni+1
	call Trace2edgeRG(phot,v,inext,jnext,irgnext)
	tau1=v*C(phot%i,phot%j)%dens*C(phot%i,phot%j)%Kext*AU
	if((tau+tau1).gt.tau0) then
		v=v*(tau0-tau)/tau1
		tau1=tau0-tau
		call makescattfield(phot,v,iangle,phot1,sin2t,cos2t)
		phot%x=phot%x+phot%vx*v
		phot%y=phot%y+phot%vy*v
		phot%z=phot%z+phot%vz*v
		tautot=tautot+(tau0-tau)
		phot%onEdge=.false.
		return
	endif
	call makescattfield(phot,v,iangle,phot1,sin2t,cos2t)
	phot%x=phot%x+phot%vx*v
	phot%y=phot%y+phot%vy*v
	phot%z=phot%z+phot%vz*v

	phot%onEdge=.true.
	tau=tau+tau1
	tautot=tautot+tau1
	if(inext.ge.D%nR) then
		escape=.true.
		return
	endif
	if(inext.lt.0) then
		escape=.true.
		hitstar=.true.
		return
	endif
	
	if(emptylower) then
		if(jnext.eq.D%nTheta-1..and.phot%j.eq.D%nTheta-1.and.inext.eq.phot%i.and.phot%irg.eq.phot%irg) then
			escape=.true.
			return
		endif
	endif
	
	phot%i=inext
	phot%j=jnext
	phot%irg=irgnext
	goto 1

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine makescattfield(phot,v,iangle,phot1,sin2t,cos2t)
	use Parameters
	IMPLICIT NONE
	real*8 v,x,y,z,phi,phi0,sin2t(NPHISCATT),cos2t(NPHISCATT),tauplanet
	real*8 F11,F12,F22,F33,F34,F44,sI,sQ,sU,sV,Qt,Ut,vAU,P1,P2,theta
	type(photon) phot,phot1(NPHISCATT)
	integer i,j,iangle(NPHISCATT),ia,side,irg
	
	if(.not.phot%scatt) return
	
	if(scat_how.eq.1.or.makeangledependence) then
!$OMP FLUSH(C)
		C(phot%i,phot%j)%scattfield(1,0,1)=C(phot%i,phot%j)%scattfield(1,0,1)+phot%E*v*AU/2d0
!$OMP FLUSH(C)
		C(phot%i,phot%j)%scattfield(1,0,2)=C(phot%i,phot%j)%scattfield(1,0,2)+phot%E*v*AU/2d0
	else
		x=phot%x+phot%vx*v/2d0
		y=phot%y+phot%vy*v/2d0
		z=phot%z+phot%vz*v/2d0
		irg=phot%irg
		if(phot%z.gt.0d0) then
			C(phot%i,phot%j)%scattfield(irg,0,1)=C(phot%i,phot%j)%scattfield(irg,0,1)+phot%E*v*AU
		else
			C(phot%i,phot%j)%scattfield(irg,0,2)=C(phot%i,phot%j)%scattfield(irg,0,2)+phot%E*v*AU
		endif
		phi0=180d0*acos(x/sqrt(x**2+y**2))/pi
		if(y.lt.0d0) phi0=360d0-phi0

		phi=180d0/real(NPHISCATT)
		j=(real(NPHISCATT)*(phi+phi0)/360d0)
		if(phot%z.gt.0d0) then
			side=1
		else
			side=2
		endif

		vAU=v*AU/2d0
		do i=1,NPHISCATT/2
			if(j.lt.1) j=j+NPHISCATT
			if(j.gt.NPHISCATT) j=j-NPHISCATT
			ia=iangle(j)

			F11=C(phot%i,phot%j)%F%F11(ia)
			F12=C(phot%i,phot%j)%F%F12(ia)
			F22=C(phot%i,phot%j)%F%F22(ia)
			F33=C(phot%i,phot%j)%F%F33(ia)
			F34=C(phot%i,phot%j)%F%F34(ia)
			F44=C(phot%i,phot%j)%F%F44(ia)

			sI=F11*phot1(j)%E+F12*phot1(j)%Q
			if(useobspol) then
				Qt=F12*phot1(j)%E+F22*phot1(j)%Q
				Ut=F34*phot1(j)%V+F33*phot1(j)%U
				sQ=Qt*cos2t(j)+Ut*sin2t(j)
				sU=-Qt*sin2t(j)+Ut*cos2t(j)
				sV=-F34*phot1(j)%U+F44*phot1(j)%V
			endif

!$OMP FLUSH(C)
			C(phot%i,phot%j)%scattfield(irg,i,side)=C(phot%i,phot%j)%scattfield(irg,i,side)+sI*vAU
			if(useobspol) then
!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattQ(irg,i,side)=C(phot%i,phot%j)%scattQ(irg,i,side)+sQ*vAU
!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattU(irg,i,side)=C(phot%i,phot%j)%scattU(irg,i,side)+sU*vAU
!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattV(irg,i,side)=C(phot%i,phot%j)%scattV(irg,i,side)+sV*vAU
			endif
			j=j+1
		enddo
		do i=NPHISCATT/2+1,NPHISCATT
			if(j.lt.1) j=j+NPHISCATT
			if(j.gt.NPHISCATT) j=j-NPHISCATT
			ia=iangle(j)

			F11=C(phot%i,phot%j)%F%F11(ia)
			F12=C(phot%i,phot%j)%F%F12(ia)
			F22=C(phot%i,phot%j)%F%F22(ia)
			F33=C(phot%i,phot%j)%F%F33(ia)
			F34=C(phot%i,phot%j)%F%F34(ia)
			F44=C(phot%i,phot%j)%F%F44(ia)

			sI=F11*phot1(j)%E+F12*phot1(j)%Q
			if(useobspol) then
				Qt=F12*phot1(j)%E+F22*phot1(j)%Q
				Ut=F34*phot1(j)%V+F33*phot1(j)%U
				sQ=Qt*cos2t(j)+Ut*sin2t(j)
				sU=-Qt*sin2t(j)+Ut*cos2t(j)
				sV=-F34*phot1(j)%U+F44*phot1(j)%V
			endif

!$OMP FLUSH(C)
			C(phot%i,phot%j)%scattfield(irg,NPHISCATT+1-i,side)=C(phot%i,phot%j)%scattfield(irg,NPHISCATT+1-i,side)+sI*vAU
			if(useobspol) then
!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattQ(irg,NPHISCATT+1-i,side)=C(phot%i,phot%j)%scattQ(irg,NPHISCATT+1-i,side)+sQ*vAU
!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattU(irg,NPHISCATT+1-i,side)=C(phot%i,phot%j)%scattU(irg,NPHISCATT+1-i,side)-sU*vAU
!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattV(irg,NPHISCATT+1-i,side)=C(phot%i,phot%j)%scattV(irg,NPHISCATT+1-i,side)-sV*vAU
			endif
			j=j+1
		enddo
		do i=1,nplanets
			if(Planets(i)%i.eq.phot%i.and.Planets(i)%j.eq.phot%j) then
				j=(real(NPHISCATT)*(phi+phi0-Planets(i)%phi)/360d0)
				if(j.lt.1) j=j+NPHISCATT
				if(j.gt.NPHISCATT) j=j-NPHISCATT
				ia=iangle(j)

				F11=Planets(i)%F%F11(ia)
				F12=Planets(i)%F%F12(ia)
				F22=Planets(i)%F%F22(ia)
				F33=Planets(i)%F%F33(ia)
				F34=Planets(i)%F%F34(ia)
				F44=Planets(i)%F%F44(ia)

				sI=F11*phot1(j)%E+F12*phot1(j)%Q
				if(useobspol) then
					Qt=F12*phot1(j)%E+F22*phot1(j)%Q
					Ut=F34*phot1(j)%V+F33*phot1(j)%U
					sQ=Qt*cos2t(j)+Ut*sin2t(j)
					sU=-Qt*sin2t(j)+Ut*cos2t(j)
					sV=-F34*phot1(j)%U+F44*phot1(j)%V
				endif

				tauplanet=v*AU*pi*Planets(i)%R**2*Planets(i)%A/C(phot%i,phot%j)%V

				if(Planets(i)%phi.lt.180d0) then
					Planets(i)%E=Planets(i)%E+sI*tauplanet
					if(useobspol) then
						Planets(i)%Q=Planets(i)%Q+sQ*tauplanet
						Planets(i)%U=Planets(i)%U-sU*tauplanet
						Planets(i)%V=Planets(i)%V-sV*tauplanet
					endif
				else
					Planets(i)%E=Planets(i)%E+sI*tauplanet
					if(useobspol) then
						Planets(i)%Q=Planets(i)%Q+sQ*tauplanet
						Planets(i)%U=Planets(i)%U+sU*tauplanet
						Planets(i)%V=Planets(i)%V+sV*tauplanet
					endif
				endif
			endif
		enddo
	endif
	return
	end		

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine MakeRotateStokes(phot,x,y,z,sin2t,cos2t)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 x,y,z,cost,sint,Q,U,cos2t,sin2t,een
	
	cost=phot%Sx*x+phot%Sy*y+phot%Sz*z
	sint=(x-
     &	(phot%Sx*(phot%vy**2+phot%vz**2)-phot%vx*(phot%vy*phot%Sy+phot%vz*phot%Sz))*cost)/
     &	(phot%vy*phot%Sz-phot%vz*phot%Sy)
	
	cos2t=cost**2-sint**2
	sin2t=2d0*sint*cost

	een=sqrt(sin2t**2+cos2t**2)
	sin2t=sin2t/een
	cos2t=cos2t/een

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine MakeStarScatter(Estar,Nphot)
	use Parameters
	IMPLICIT NONE
	real*8 Estar,r,ct,inp,sin2t(NPHISCATT),cos2t(NPHISCATT)
	integer Nphot,i,j,iphot,iangle(NPHISCATT),ia,inext,jnext,side,jj,irg,irgnext
	real*8 v,x,y,z,phi,phi0,Qt,sQ,sU,sI,w,tau,F11,F12,vAU,P
	type(photon) phot,phot1(NPHISCATT),phot2
	real*8 theta,tauplanet,Emin

	write(*,'("Single scattered starlight:",i8," photon packages")') Nphot
	write(9,'("Single scattered starlight:",i8," photon packages")') Nphot

	Emin=1d-10*Estar/real(Nphot)

	call tellertje(1,100)
!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(phot,x,y,z,r,tautot,tau,xsn,ysn,zsn,
!$OMP&   inp,ct,theta,iangle,v,w,side,irg,phi0,phi,j,vAU,F11,F12,sI,Qt,sQ,sU,ia,
!$OMP&   phot1,inext,jnext,irgnext,cos2t,sin2t)
!$OMP& SHARED(scat_how,C,D,Estar,makeangledependence,xin,yin,zin,
!$OMP&   ninteract,useobspol,xsf,ysf,zsf,Nphot)
!$OMP DO
	do iphot=1,Nphot
!$OMP CRITICAL
		call tellertje(iphot+1,Nphot+2)
!$OMP END CRITICAL
		phot%E=Estar/real(Nphot)
		phot%scatt=.false.

		call randomdirection(phot%x,phot%y,phot%z)
		phot%x=D%R(0)*phot%x
		phot%y=D%R(0)*phot%y
		phot%z=D%R(0)*phot%z
		call randomdirection(phot%vx,phot%vy,phot%vz)
		inp=(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz)
		if(inp.lt.0d0) then
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
		if(C(phot%i,phot%j)%nrg.gt.1) then
			theta=acos(ct)
			phot%irg=int(real(C(phot%i,phot%j)%nrg)*(theta-D%thet(phot%j))/(D%thet(phot%j+1)-D%thet(phot%j)))+1
			if(phot%irg.gt.C(phot%i,phot%j)%nrg) phot%irg=C(phot%i,phot%j)%nrg
			if(phot%irg.lt.1) phot%irg=1
		else
			phot%irg=1
		endif

		if(scat_how.ne.1.and..not.makeangledependence) then
		do i=1,NPHISCATT
			iangle(i)=180d0*acos(phot%vx*xsf(i)+phot%vy*ysf(i)
     &						+phot%vz*zsf(i))/pi+1
			xsn(i)=phot%vy*zsf(i)-phot%vz*ysf(i)
			ysn(i)=phot%vz*xsf(i)-phot%vx*zsf(i)
			zsn(i)=phot%vx*ysf(i)-phot%vy*xsf(i)
			r=sqrt(xsn(i)**2+ysn(i)**2+zsn(i)**2)
			xsn(i)=xsn(i)/r
			ysn(i)=ysn(i)/r
			zsn(i)=zsn(i)/r
			phot1(i)=phot
			phot1(i)%Sx=xsn(i)
			phot1(i)%Sy=ysn(i)
			phot1(i)%Sz=zsn(i)
			phot1(i)%vx=xsf(i)
			phot1(i)%vy=ysf(i)
			phot1(i)%vz=zsf(i)
			call MakeRotateStokes(phot1(i),xin(i),yin(i),zin(i),sin2t(i),cos2t(i))
			phot1(i)%vx=phot%vx
			phot1(i)%vy=phot%vy
			phot1(i)%vz=phot%vz
   		enddo
		endif

1		continue
		call Trace2edgeRG(phot,v,inext,jnext,irgnext)
		tau=C(phot%i,phot%j)%dens*C(phot%i,phot%j)%Kext*AU

!$OMP FLUSH(C)
		if(phot%E.gt.(1d-3*Estar/real(Nphot))) C(phot%i,phot%j)%Ni=C(phot%i,phot%j)%Ni+1

		if((tau*v).gt.1d-5) then
			w=(1d0-exp(-tau*v))/tau
		else
			w=v
		endif

		if(scat_how.eq.1.or.makeangledependence) then
!$OMP FLUSH(C)
			C(phot%i,phot%j)%scattfield(1,0,1)=C(phot%i,phot%j)%scattfield(1,0,1)+phot%E*w*AU/2d0
!$OMP FLUSH(C)
			C(phot%i,phot%j)%scattfield(1,0,2)=C(phot%i,phot%j)%scattfield(1,0,2)+phot%E*w*AU/2d0
		else
			x=phot%x+phot%vx*v/2d0
			y=phot%y+phot%vy*v/2d0
			z=phot%z+phot%vz*v/2d0
			irg=phot%irg
			if(phot%z.gt.0d0) then
!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattfield(irg,0,1)=C(phot%i,phot%j)%scattfield(irg,0,1)+phot%E*w*AU
				side=1
			else
!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattfield(irg,0,2)=C(phot%i,phot%j)%scattfield(irg,0,2)+phot%E*w*AU
				side=2
			endif
			
			phi0=180d0*acos(x/sqrt(x**2+y**2))/pi
			if(y.lt.0d0) phi0=360d0-phi0

			phi=180d0/real(NPHISCATT)
			j=(real(NPHISCATT)*(phi+phi0)/360d0)

			vAU=v*AU/2d0
			do i=1,NPHISCATT/2
				if(j.lt.1) j=j+NPHISCATT
				if(j.gt.NPHISCATT) j=j-NPHISCATT
				ia=iangle(j)

				F11=C(phot%i,phot%j)%F%F11(ia)
				F12=C(phot%i,phot%j)%F%F12(ia)

				sI=F11*phot1(j)%E
				if(useobspol) then
					Qt=F12*phot1(j)%E
					sQ=Qt*cos2t(j)
					sU=-Qt*sin2t(j)
				endif

!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattfield(irg,i,side)=C(phot%i,phot%j)%scattfield(irg,i,side)+sI*vAU
				if(useobspol) then
!$OMP FLUSH(C)
					C(phot%i,phot%j)%scattQ(irg,i,side)=C(phot%i,phot%j)%scattQ(irg,i,side)+sQ*vAU
!$OMP FLUSH(C)
					C(phot%i,phot%j)%scattU(irg,i,side)=C(phot%i,phot%j)%scattU(irg,i,side)+sU*vAU
				endif
				j=j+1
			enddo
			do i=NPHISCATT/2+1,NPHISCATT
				if(j.lt.1) j=j+NPHISCATT
				if(j.gt.NPHISCATT) j=j-NPHISCATT
				ia=iangle(j)

				F11=C(phot%i,phot%j)%F%F11(ia)
				F12=C(phot%i,phot%j)%F%F12(ia)

				sI=F11*phot1(j)%E
				if(useobspol) then
					Qt=F12*phot1(j)%E
					sQ=Qt*cos2t(j)
					sU=-Qt*sin2t(j)
				endif

!$OMP FLUSH(C)
				C(phot%i,phot%j)%scattfield(irg,NPHISCATT+1-i,side)=C(phot%i,phot%j)%scattfield(irg,NPHISCATT+1-i,side)+sI*vAU
				if(useobspol) then
!$OMP FLUSH(C)
					C(phot%i,phot%j)%scattQ(irg,NPHISCATT+1-i,side)=C(phot%i,phot%j)%scattQ(irg,NPHISCATT+1-i,side)+sQ*vAU
!$OMP FLUSH(C)
					C(phot%i,phot%j)%scattU(irg,NPHISCATT+1-i,side)=C(phot%i,phot%j)%scattU(irg,NPHISCATT+1-i,side)-sU*vAU
				endif
				j=j+1
			enddo
		endif

		w=exp(-tau*v)
		phot%E=phot%E*w
		phot1(1:NPHISCATT)%E=phot1(1:NPHISCATT)%E*w
		phot%x=phot%x+phot%vx*v
		phot%y=phot%y+phot%vy*v
		phot%z=phot%z+phot%vz*v

		phot%onEdge=.true.

		if(inext.ge.D%nR.or.inext.lt.0) goto 2
c		if(phot%E.lt.Emin) goto 2
		phot%i=inext
		phot%j=jnext
		phot%irg=irgnext
		goto 1

2		continue
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(100,100)

	if(nplanets.gt.0) then
	write(*,'("Planet scattered starlight")')
	write(9,'("Planet scattered starlight")')

	do iphot=1,D%nTheta-1
	call tellertje(iphot,D%nTheta-1)
	do side=1,2

c		phot%E=Estar*(D%thet(iphot+1)-D%thet(iphot))/(pi/2d0)
		phot%E=Estar*(D%Theta(iphot)-D%Theta(iphot+1))

		phot%scatt=.false.

		phot%z=D%R(0)*cos(D%theta_av(iphot))
		if(side.eq.2) phot%z=-phot%z
		phot%x=sqrt((D%R(0)**2-phot%z**2)/2d0)
		phot%y=phot%x
		phot%vx=phot%x/D%R(0)
		phot%vy=phot%y/D%R(0)
		phot%vz=phot%z/D%R(0)

		phot%edgeNr=1
		phot%onEdge=.true.
		ct=abs(phot%z)/D%R(0)
		do j=1,D%nTheta-1
			if(ct.lt.D%Theta(j).and.ct.gt.D%Theta(j+1)) then
				phot%i=0
				phot%j=j
			endif
		enddo
		if(C(phot%i,phot%j)%nrg.gt.1) then
			theta=acos(ct)
			phot%irg=int(real(C(phot%i,phot%j)%nrg)*(theta-D%thet(phot%j))/(D%thet(phot%j+1)-D%thet(phot%j)))+1
			if(phot%irg.gt.C(phot%i,phot%j)%nrg) phot%irg=C(phot%i,phot%j)%nrg
			if(phot%irg.lt.1) phot%irg=1
		else
			phot%irg=1
		endif

		if(scat_how.ne.1.and..not.makeangledependence) then
		do i=1,NPHISCATT
			iangle(i)=180d0*acos(phot%vx*xsf(i)+phot%vy*ysf(i)
     &						+phot%vz*zsf(i))/pi+1
			xsn(i)=phot%vy*zsf(i)-phot%vz*ysf(i)
			ysn(i)=phot%vz*xsf(i)-phot%vx*zsf(i)
			zsn(i)=phot%vx*ysf(i)-phot%vy*xsf(i)
			r=sqrt(xsn(i)**2+ysn(i)**2+zsn(i)**2)
			xsn(i)=xsn(i)/r
			ysn(i)=ysn(i)/r
			zsn(i)=zsn(i)/r
			phot1(i)=phot
			phot1(i)%Sx=xsn(i)
			phot1(i)%Sy=ysn(i)
			phot1(i)%Sz=zsn(i)
			phot1(i)%vx=xsf(i)
			phot1(i)%vy=ysf(i)
			phot1(i)%vz=zsf(i)
			call MakeRotateStokes(phot1(i),xin(i),yin(i),zin(i),sin2t(i),cos2t(i))
			phot1(i)%vx=phot%vx
			phot1(i)%vy=phot%vy
			phot1(i)%vz=phot%vz
   		enddo
		endif

3		continue
		call Trace2edgeRG(phot,v,inext,jnext,irgnext)
		tau=C(phot%i,phot%j)%dens*C(phot%i,phot%j)%Kext*AU

		if((tau*v).gt.1d-5) then
			w=(1d0-exp(-tau*v))/tau
		else
			w=v
		endif

		x=phot%x+phot%vx*v/2d0
		y=phot%y+phot%vy*v/2d0
		z=phot%z+phot%vz*v/2d0
		irg=phot%irg

		phi0=180d0*acos(x/sqrt(x**2+y**2))/pi
		if(y.lt.0d0) phi0=360d0-phi0

		phi=180d0/real(NPHISCATT)
		j=(real(NPHISCATT)*(phi+phi0)/360d0)

		vAU=v*AU/2d0
		do i=1,nplanets
			if(Planets(i)%i.eq.phot%i.and.Planets(i)%j.eq.phot%j.and.Planets(i)%k.eq.side) then
				if(Planets(i)%phi.lt.180d0) then
					j=(real(NPHISCATT)*(phi+phi0+Planets(i)%phi)/360d0)
				else
					j=(real(NPHISCATT)*(phi+phi0+360d0-Planets(i)%phi)/360d0)
				endif
				if(j.lt.1) j=j+NPHISCATT
				if(j.gt.NPHISCATT) j=j-NPHISCATT
				ia=iangle(j)
				
				F11=Planets(i)%F%F11(ia)
				F12=Planets(i)%F%F12(ia)

				sI=F11*phot1(j)%E
				if(useobspol) then
					Qt=F12*phot1(j)%E
					sQ=Qt*cos2t(j)
					sU=-Qt*sin2t(j)
				endif

				tauplanet=v*AU*pi*Planets(i)%R**2*Planets(i)%A/C(phot%i,phot%j)%V

				if(Planets(i)%phi.lt.180d0) then
					Planets(i)%E=Planets(i)%E+sI*tauplanet
					if(useobspol) then
						Planets(i)%Q=Planets(i)%Q+sQ*tauplanet
						Planets(i)%U=Planets(i)%U+sU*tauplanet
					endif
				else
					Planets(i)%E=Planets(i)%E+sI*tauplanet
					if(useobspol) then
						Planets(i)%Q=Planets(i)%Q+sQ*tauplanet
						Planets(i)%U=Planets(i)%U-sU*tauplanet
					endif
				endif
			endif
		enddo

		w=exp(-tau*v)
		phot%E=phot%E*w
		phot1(1:NPHISCATT)%E=phot1(1:NPHISCATT)%E*w
		phot%x=phot%x+phot%vx*v
		phot%y=phot%y+phot%vy*v
		phot%z=phot%z+phot%vz*v

		phot%onEdge=.true.

		if(inext.ge.D%nR.or.inext.lt.0) goto 4
		phot%i=inext
		phot%j=jnext
		phot%irg=irgnext
		goto 3

4		continue
	enddo
	enddo
	endif

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine TracePath(image,angle,nphi,nj0,lam0)
	use Parameters
	IMPLICIT NONE
	character*500 input,tmp,specfile
	integer i,j,k,l,t1,t2,nstar,ninner
	real*8 tau,angle,spec(nlam),distance,tottime,ct,lam0
	real*8 ph,inp,determineT,lam1,lam2,v
	integer starttime,stoptime,starttrace,cr,l1,l2,nj0,ip,jp
	type(Photon) phot,photcount
	integer ilam,MinPhot,nabs,nj,nj2,nphi
	real*8 x,y,z,phi,theta,Albedo,w1,w2,fact(nlam),Rad,r
	logical hitstar
	real*8 extstar(nlam),R01(D%nTheta)
	type(RPhiImage) image

	nstar=30
	ninner=0
	if(inner_gas) ninner=30

	write(*,'("Creating photon paths for image")')
	write(*,'("Inclination angle:",f8.1)') 180d0*angle/pi

	write(9,'("Creating photon paths for image")')
	write(9,'("Inclination angle:",f8.1)') 180d0*angle/pi

	if(scat_how.eq.2) then
	do i=1,NPHISCATT
		phi=pi*(real(i)-0.5d0)/(real(NPHISCATT)/2d0)
		xsf(i)=0d0
		ysf(i)=0d0
		zsf(i)=1d0
		call rotate(xsf(i),ysf(i),zsf(i),0d0,1d0,0d0,angle)
		call rotate(xsf(i),ysf(i),zsf(i),0d0,0d0,1d0,phi)
		r=sqrt(xsf(i)**2+ysf(i)**2+zsf(i)**2)
		xsf(i)=xsf(i)/r
		ysf(i)=ysf(i)/r
		zsf(i)=zsf(i)/r
		xin(i)=0d0
		yin(i)=1d0
		zin(i)=0d0
		call rotate(xin(i),yin(i),zin(i),0d0,0d0,1d0,phi)
		r=sqrt(xin(i)**2+yin(i)**2+zin(i)**2)
		xin(i)=xin(i)/r
		yin(i)=yin(i)/r
		zin(i)=zin(i)/r
	enddo
	endif


	image%nphi=nphi

	image%angle=angle

c	image%nr=0
c	nj=nj0*50
c	do j=1,nj
c		image%nr=image%nr+1
c	enddo
c	nj=nj0
c	do i=1,D%nR-1
c	do j=1,nj
c		image%nr=image%nr+1
c	enddo
c	enddo
c	image%nr=image%nr+1

	if(nj0.gt.0) then
		image%nr=nj0*(D%nR-1)+(D%nTheta-1)*2*(D%nRfix+2)+1+nstar+ninner
	else
		image%nr=(D%nR-2)/(-nj0)+1+(D%nTheta-1)*2*(D%nRfix+2)+1+nstar+ninner
	endif

	if(scat_how.eq.2) then
		image%nr=image%nr+180*(D%nRfix+1)
	endif
	
	allocate(image%image(image%nr,image%nphi))
	if(scat_how.eq.2) then
		allocate(image%imageQ(image%nr,image%nphi))
		allocate(image%imageU(image%nr,image%nphi))
		allocate(image%imageV(image%nr,image%nphi))
	endif
	allocate(image%Phi(image%nphi))
	allocate(image%R(image%nr))
	allocate(image%p(image%nr,image%nphi))
c	do i=1,image%nr
c		do j=1,image%nphi
c			allocate(image%p(i,j)%i((D%nR+2)*2+(D%nTheta+2)*2+360))
c			allocate(image%p(i,j)%j((D%nR+2)*2+(D%nTheta+2)*2+360))
c			allocate(image%p(i,j)%k((D%nR+2)*2+(D%nTheta+2)*2+360))
c			allocate(image%p(i,j)%irg((D%nR+2)*2+(D%nTheta+2)*2+360))
c			allocate(image%p(i,j)%v((D%nR+2)*2+(D%nTheta+2)*2+360))
c			allocate(image%p(i,j)%phi1((D%nR+2)*2+(D%nTheta+2)*2+360))
c			allocate(image%p(i,j)%phi2((D%nR+2)*2+(D%nTheta+2)*2+360))
c			allocate(image%p(i,j)%jphi1((D%nR+2)*2+(D%nTheta+2)*2+360))
c			allocate(image%p(i,j)%jphi2((D%nR+2)*2+(D%nTheta+2)*2+360))
c		enddo
c	enddo
	
	write(*,'("Number of radial points:",i10)') image%nr
	write(*,'("Number of angles:       ",i10)') image%nphi
	write(9,'("Number of radial points:",i10)') image%nr
	write(9,'("Number of angles:       ",i10)') image%nphi

c	image%nr=0
c	nj=nj0*50
c	do j=1,nj
c		image%nr=image%nr+1
c		image%R(image%nr)=D%R(1)*(real(j)-0.999d0)/real(nj)
c	enddo
c	nj=nj0
c	do i=1,D%nR-1
c	do j=1,nj
c		image%nr=image%nr+1
c		image%R(image%nr)=D%R(i)+(D%R(i+1)-D%R(i))*(real(j)-0.5d0)/real(nj)
c	enddo
c	enddo
c	image%nr=image%nr+1
c	image%R(image%nr)=D%R(D%nR)*0.9999

	call tau01R(R01,lam0)
	nj=0
	do j=1,D%nTheta-1
		nj=nj+1
		image%R(nj)=R01(D%nTheta-1)*abs(cos(pi/2d0-D%theta_av(j)+angle))
		nj=nj+1
		image%R(nj)=abs(R01(D%nTheta-1)*cos(D%theta_av(j)-pi/2d0+angle))
	enddo
	do j=1,D%nTheta-1
		nj=nj+1
		image%R(nj)=R01(D%nTheta-1)*abs(cos(pi/2d0-D%theta_av(j)))
		nj=nj+1
		image%R(nj)=abs(R01(D%nTheta-1)*cos(D%theta_av(j)-pi/2d0))
	enddo
	if(scat_how.eq.2) then
		do j=1,180
			theta=real(j)*(pi/2d0)/181d0
			nj=nj+1
			image%R(nj)=D%R(1)*abs(cos(pi/2d0-theta))
			do i=1,D%nRfix
				nj=nj+1
				image%R(nj)=D%Rfix(i)*abs(cos(pi/2d0-theta))
			enddo
		enddo
	endif

	do i=1,D%nRfix
	do j=1,D%nTheta-1
		nj=nj+1
		image%R(nj)=D%Rfix(i)*abs(cos(pi/2d0-D%theta_av(j)+angle))
		nj=nj+1
		image%R(nj)=abs(D%Rfix(i)*cos(D%theta_av(j)-pi/2d0+angle))
	enddo
	enddo
	call sort(image%R,nj)
	if(nj0.gt.0) then
		do i=1,D%nR-1
		do j=1,nj0
			nj=nj+1
			image%R(nj)=D%R(i)+(D%R(i+1)-D%R(i))*(real(j)-0.5d0)/real(nj0)
		enddo
		enddo
	else
		do i=1,D%nR-1,-nj0
			nj=nj+1
			image%R(nj)=10d0**(log10(D%R(i)*D%R(i+1))/2d0)
		enddo
	endif
	do i=1,nstar-1
		nj=nj+1
		image%R(nj)=0.999d0*D%R(0)*real(i-1)/real(nstar-2)
	enddo
	nj=nj+1
	image%R(nj)=D%R(0)*1.001

	if(inner_gas) then
		do i=1,ninner
			nj=nj+1
			image%R(nj)=D%R(0)+(D%R(1)-D%R(0))*real(i)/real(ninner+1)
		enddo
	endif

	image%R(image%nr)=D%R(D%nR)*0.9999
	call sort(image%R,image%nr)

	k=0
1	continue
	k=k+1
	do i=image%nr,2,-1
		if(image%R(i).eq.image%R(i-1)) then
			if(i.ne.image%nr) then
				image%R(i)=(image%R(i+1)+image%R(i))/2d0
			else
				image%R(i-1)=(image%R(i-1)+image%R(i))/2d0
			endif
			if(k.lt.image%nr*10) goto 1
		endif
	enddo
	
	call tellertje(1,100)
!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(Rad,phot,ct,j,k,theta,photcount)
!$OMP& SHARED(C,D,image,angle)
!$OMP DO
	do i=1,image%nr
!$OMP CRITICAL
	call tellertje(i+1,image%nr+2)
!$OMP END CRITICAL
	do k=1,image%nPhi
		image%phi(k)=1d-5+pi*(real(k)-0.5)/real(image%nPhi)
		Rad=image%R(i)
		phot%x=Rad*cos(image%phi(k))
		phot%y=Rad*sin(image%phi(k))
		phot%z=sqrt(D%R(D%nR)**2-Rad**2)
		phot%edgeNr=2
		phot%onEdge=.true.
		call rotate(phot%x,phot%y,phot%z,0d0,1d0,0d0,angle)
		ct=abs(phot%z)/D%R(D%nR)
		phot%i=D%nR-1
		do j=1,D%nTheta-1
			if(ct.lt.D%Theta(j).and.ct.ge.D%Theta(j+1)) then
				phot%j=j
			endif
		enddo

		if(C(phot%i,phot%j)%nrg.gt.1) then
			theta=acos(ct)
			phot%irg=int(real(C(phot%i,phot%j)%nrg)*(theta-D%thet(phot%j))/(D%thet(phot%j+1)-D%thet(phot%j)))+1
			if(phot%irg.gt.C(phot%i,phot%j)%nrg) phot%irg=C(phot%i,phot%j)%nrg
			if(phot%irg.lt.1) phot%irg=1
		else
			phot%irg=1
		endif

		phot%vx=0d0
		phot%vy=0d0
		phot%vz=-1d0
		call rotate(phot%vx,phot%vy,phot%vz,0d0,1d0,0d0,angle)
		photcount=phot
		call count2dpath(photcount,image%p(i,k))
		allocate(image%p(i,k)%i(image%p(i,k)%n))
		allocate(image%p(i,k)%j(image%p(i,k)%n))
		allocate(image%p(i,k)%k(image%p(i,k)%n))
		allocate(image%p(i,k)%irg(image%p(i,k)%n))
		allocate(image%p(i,k)%v(image%p(i,k)%n))
		allocate(image%p(i,k)%phi1(image%p(i,k)%n))
		allocate(image%p(i,k)%phi2(image%p(i,k)%n))
		allocate(image%p(i,k)%jphi1(image%p(i,k)%n))
		allocate(image%p(i,k)%jphi2(image%p(i,k)%n))
		allocate(image%p(i,k)%rad(image%p(i,k)%n))

		call trace2dpath(phot,image%p(i,k))
	enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(100,100)

	return
	end
	
c-----------------------------------------------------------------------
c This subroutine outputs the tau=1 height of the disk at visual
c wavelength (lambda=0.55 micron) when viewed from the central star
c to the file filename
c-----------------------------------------------------------------------
	subroutine tau01R(R01,lam0)
	use Parameters
	IMPLICIT NONE
	real*8 tau1x(D%nTheta,10),tau1z(D%nTheta,10),r(D%nR*10)
	real*8 ct,tau,lr0,lr1,R01(D%nTheta)
	real*8 Mtot,Mtot2,Fsc,p,p0,p1,Kext,lam0
	logical escape,hitstar
	type(photon) phot
	integer i,j,n,ii
	
	phot%nr=1
	tautot=0d0
	do i=1,D%nTheta-1
		r(i)=D%R(1)
		ct=(D%Theta(i)+D%Theta(i+1))/2d0
		phot%x=sqrt(1d0-ct**2)*r(i)
		phot%y=0d0
		phot%z=ct*r(i)

		phot%edgenr=1
		phot%onEdge=.true.
		phot%vx=phot%x/r(i)
		phot%vy=0d0
		phot%vz=phot%z/r(i)
		tau=0.1d0
		phot%lam=lam0
		phot%i=1
		phot%j=i
		do j=1,nlam-1
			if(phot%lam.ge.lam(j).and.phot%lam.le.lam(j+1)) then
				phot%ilam1=j
				phot%ilam2=j+1
				phot%wl1=(lam(j+1)-phot%lam)/(lam(j+1)-lam(j))
				phot%wl2=(phot%lam-lam(j))/(lam(j+1)-lam(j))
				phot%nu=nu(j)*phot%wl1+nu(j+1)*phot%wl2
			endif
		enddo
		phot%E=1d0
		call Trace2D(phot,tau,escape,hitstar,.false.)
		if(escape) phot%z=0d0!-1d0/0d0!sqrt(D%R(D%nR)**2-phot%x**2)
		R01(i)=sqrt(phot%x**2+phot%z**2)
		if(R01(i).gt.D%R(D%nR-1)) R01(i)=D%R(1)
	enddo

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine trace2dpath(phot,p)
	use Parameters
	IMPLICIT NONE
	type(Photon) phot
	real*8 tau(nlam),tau_a(nlam),tau_e(nlam),v,fact(nlam)
	real*8 Kext(nlam),Kabs(nlam),emis(nlam),spec(nlam)
	real*8 absfrac(nlam),wT1,wT2,x,y,z,phi
	integer inext,jnext,ntrace,iT,i,j,l1,l2,irgnext
	logical hitstar
	type(path) p

	p%hitstar=.false.

	p%n=0
	p%x=phot%x
	p%y=phot%y
	p%z=phot%z
	p%vx=phot%vx
	p%vy=phot%vy
	p%vz=phot%vz

1	continue

	call Trace2edgeRG(phot,v,inext,jnext,irgnext)
	phot%onEdge=.true.

	p%n=p%n+1
	p%v(p%n)=v
	
	if(emptylower.and.phot%z.gt.0d0) p%v(p%n)=0d0

	p%i(p%n)=phot%i
	p%j(p%n)=phot%j
	p%irg(p%n)=phot%irg
	
	x=phot%x
	y=phot%y
	z=phot%z

	p%phi1(p%n)=180d0*acos(x/sqrt(x**2+y**2))/pi
	if(y.lt.0d0) p%phi1(p%n)=-p%phi1(p%n)
	p%jphi1(p%n)=real(NPHISCATT)*p%phi1(p%n)/360d0
	p%jphi1(p%n)=p%jphi1(p%n)+1
	if(p%jphi1(p%n).lt.1) p%jphi1(p%n)=1-p%jphi1(p%n)
	if(p%jphi1(p%n).gt.NPHISCATT/2) p%jphi1(p%n)=NPHISCATT+1-p%jphi1(p%n)
	x=phot%x+phot%vx*v
	y=phot%y+phot%vy*v
	z=phot%z+phot%vz*v
	p%rad(p%n)=sqrt(x*x+y*y+z*z)
	if(z.gt.0d0) then
		p%k(p%n)=1
	else
		p%k(p%n)=2
	endif
	p%phi2(p%n)=180d0*acos(x/sqrt(x**2+y**2))/pi
	if(y.lt.0d0) p%phi2(p%n)=-p%phi2(p%n)
	p%jphi2(p%n)=real(NPHISCATT)*p%phi2(p%n)/360d0
	p%jphi2(p%n)=p%jphi2(p%n)+1
	if(p%jphi2(p%n).lt.1) p%jphi2(p%n)=1-p%jphi2(p%n)
	if(p%jphi2(p%n).gt.NPHISCATT/2) p%jphi2(p%n)=NPHISCATT+1-p%jphi2(p%n)

	phot%x=phot%x+phot%vx*v
	phot%y=phot%y+phot%vy*v
	phot%z=phot%z+phot%vz*v

	if(inext.ge.D%nR) then
		return
	endif
	if(inext.lt.0) then
		p%hitstar=.true.
		return
	endif

	phot%i=inext
	phot%j=jnext
	phot%irg=irgnext

	goto 1

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine count2dpath(phot,p)
	use Parameters
	IMPLICIT NONE
	type(Photon) phot
	real*8 tau(nlam),tau_a(nlam),tau_e(nlam),v,fact(nlam)
	real*8 Kext(nlam),Kabs(nlam),emis(nlam),spec(nlam)
	real*8 absfrac(nlam),wT1,wT2,x,y,z,phi
	integer inext,jnext,ntrace,iT,i,j,l1,l2,irgnext
	logical hitstar
	type(path) p

	p%hitstar=.false.

	p%n=0

1	continue

	call Trace2edgeRG(phot,v,inext,jnext,irgnext)
	phot%onEdge=.true.

	p%n=p%n+1
	
	phot%x=phot%x+phot%vx*v
	phot%y=phot%y+phot%vy*v
	phot%z=phot%z+phot%vz*v

	if(inext.ge.D%nR) then
		return
	endif
	if(inext.lt.0) then
		p%hitstar=.true.
		return
	endif

	phot%i=inext
	phot%j=jnext
	phot%irg=irgnext

	goto 1

	return
	end


!c-----------------------------------------------------------------------
!c This subroutine determines the distance a photon has to travel to
!c the next border of the cell. The output variables are v, the distance
!c inext, the indicator for the next radial cell, and jnext, for the
!c next theta cell.
!c This version includes the regridded theta cells for highly anisotropic
!c scattering.
!c-----------------------------------------------------------------------
	subroutine Trace2edgeRG(phot,v,inext,jnext,irgnext)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 b,v,r,R1,R2,T1,T2,vR1,vR2,vT1,vT2,theta
	integer inext,jnext,irgnext
	logical hitR1,hitR2,hitR,hitT1,hitT2,hitT,hitTsame

	r=phot%x**2+phot%y**2+phot%z**2
	R1=C(phot%i,phot%j)%xedge2(1)
	R2=C(phot%i,phot%j)%xedge2(2)
	T1=C(phot%i,phot%j)%thetrg(phot%irg)
	T2=C(phot%i,phot%j)%thetrg(phot%irg+1)

c	if(C(phot%i,phot%j)%nrg.gt.1) then
c		if(phot%irg.eq.1) then
c			T1=C(phot%i,phot%j)%xedge2(3)
c		else
c			theta=D%thet(phot%j)+real(phot%irg-1)*(D%thet(phot%j+1)-D%thet(phot%j))/real(C(phot%i,phot%j)%nrg)
c			T1=cos(theta)**2
c		endif
c		if(phot%irg.eq.C(phot%i,phot%j)%nrg) then
c			T2=C(phot%i,phot%j)%xedge2(4)
c		else
c			theta=D%thet(phot%j)+real(phot%irg)*(D%thet(phot%j+1)-D%thet(phot%j))/real(C(phot%i,phot%j)%nrg)
c			T2=cos(theta)**2
c		endif
c	else
c		T1=C(phot%i,phot%j)%xedge2(3)
c		T2=C(phot%i,phot%j)%xedge2(4)
c	endif

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
		if(phot%j.ne.(D%nTheta-1).or.phot%irg.ne.C(phot%i,phot%j)%nrg) then
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
		if(phot%j.ne.1.or.phot%irg.ne.1) then
			hitT1=hitT(phot,T1,r,b,vT1)
		else
			hitT1=.false.
			vT1=1d200
		endif
		if(phot%j.ne.(D%nTheta-1).or.phot%irg.ne.C(phot%i,phot%j)%nrg) then
			hitT2=hitTsame(phot,T2,r,b,vT2)
		else
			hitT2=.false.
			vT2=1d200
		endif
	endif
	endif

	if(.not.hitR1.and..not.hitR2.and..not.hitT1.and..not.hitT2) then
		print*,'nothing to hit!',vT1,vT2,vR1,vR2
		hitR1=hitR(phot,R1,r,b,vR1)
		hitR2=hitR(phot,R2,r,b,vR2)
		hitT1=hitT(phot,T1,r,b,vT1)
		hitT2=hitT(phot,T2,r,b,vT2)
		if(.not.hitR1.and..not.hitR2.and..not.hitT1.and..not.hitT2) print*,'still nothing'
		print*,phot%z,phot%vz,T1,T2
		print*,phot%nr
		stop
	endif

	v=1d200
	if(hitR1.and.vR1.lt.v) then
		v=vR1
		inext=phot%i-1
		jnext=phot%j
		irgnext=phot%irg
		phot%edgeNr=2
	endif
	if(hitR2.and.vR2.lt.v) then
		v=vR2
		inext=phot%i+1
		jnext=phot%j
		irgnext=phot%irg
		phot%edgeNr=1
	endif
	if(hitT1.and.vT1.lt.v) then
		v=vT1
		inext=phot%i
		if(C(phot%i,phot%j)%nrg.gt.1) then
			irgnext=phot%irg-1
			jnext=phot%j
			if(irgnext.lt.1) then
				jnext=phot%j-1
				irgnext=C(phot%i,jnext)%nrg
			endif
		else
			jnext=phot%j-1
			if(jnext.gt.0) irgnext=C(phot%i,jnext)%nrg
		endif
		phot%edgeNr=4
	endif
	if(hitT2.and.vT2.lt.v) then
		v=vT2
		inext=phot%i
		if(C(phot%i,phot%j)%nrg.gt.1) then
			irgnext=phot%irg+1
			jnext=phot%j
			if(irgnext.gt.C(phot%i,phot%j)%nrg) then
				jnext=phot%j+1
				irgnext=1
				if(jnext.eq.D%nTheta) then
					irgnext=C(phot%i,D%nTheta-1)%nrg
				endif
			endif
		else
			jnext=phot%j+1
			irgnext=1
		endif
		phot%edgeNr=3
	endif


	if(jnext.eq.D%nTheta) then
		jnext=D%nTheta-1
		irgnext=C(phot%i,D%nTheta-1)%nrg
		phot%edgeNr=4
	endif
	if(jnext.eq.0) then
		jnext=1
		irgnext=1
		phot%edgeNr=3
	endif


	return
	end	
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine trace2exit(phot,tau,withscatt)
	use Parameters
	IMPLICIT NONE
	type(Photon) phot
	real*8 tau,v
	integer inext,jnext,ntrace,i,j,ii
	logical withscatt

	tau=0d0

1	continue

	call Trace2edge(phot,v,inext,jnext)
	phot%onEdge=.true.

	phot%x=phot%x+phot%vx*v
	phot%y=phot%y+phot%vy*v
	phot%z=phot%z+phot%vz*v

	if(withscatt) then
		tau=tau+v*C(phot%i,phot%j)%dens*AU*C(phot%i,phot%j)%Kext
	else
		tau=tau+v*C(phot%i,phot%j)%dens*AU*C(phot%i,phot%j)%Kabs
	endif
		
	if(inext.ge.D%nR.or.inext.le.0) return

	phot%i=inext
	phot%j=jnext

	goto 1

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------



	real*8 function ReadMCScatt(angle0,lam0,nf)
	use Parameters
	IMPLICIT NONE
	real*8 angle0,lam0,nf
	character*500 file
	real*8,allocatable :: angle(:),scat1(:),scat2(:)
	real*8 lam1,lam2,wl1,wl2,wa1,wa2
	integer i,j,nl,nangle,ia
	real*8,allocatable :: iscat1(:),iscat2(:)
	
	write(file,'(a,"MCScattered.dat")') outdir(1:len_trim(outdir))
	open(unit=40,file=file,RECL=6000)
	read(40,*) nangle
	read(40,*) nl
	allocate(angle(nangle))
	allocate(scat1(nangle))
	allocate(scat2(nangle))
	allocate(iscat1(nangle))
	allocate(iscat2(nangle))
	do i=1,nangle
		read(40,*) angle(i)
		angle(i)=angle(i)*pi/180d0
	enddo
	ia=1
	do i=1,nangle-1
		if(angle0.ge.angle(i).and.angle0.lt.angle(i+1)) ia=i
	enddo
	if(angle0.gt.angle(nangle)) ia=nangle-1
	
	read(40,*) lam1,(scat1(j),j=1,nangle),(iscat1(j),j=1,nangle)
	if(lam0.le.lam1) then
		wa1=(angle(ia+1)-angle0)/(angle(ia+1)-angle(ia))
		wa2=(angle0-angle(ia))/(angle(ia+1)-angle(ia))
		ReadMCScatt=(wa1*scat1(ia)+wa2*scat1(ia+1))*D%distance**2/1d23
		nf=(wa1*iscat1(ia)+wa2*iscat1(ia+1))
		close(unit=40)
		return
	endif
	do i=1,nl-1
		read(40,*) lam2,(scat2(j),j=1,nangle),(iscat2(j),j=1,nangle)
		if(lam0.ge.lam1.and.lam0.lt.lam2) then
			wl1=(lam2-lam0)/(lam2-lam1)
			wl2=(lam0-lam1)/(lam2-lam1)
			wa1=(angle(ia+1)-angle0)/(angle(ia+1)-angle(ia))
			wa2=(angle0-angle(ia))/(angle(ia+1)-angle(ia))
			ReadMCScatt=(wl1*(wa1*scat1(ia)+wa2*scat1(ia+1))+wl2*(wa1*scat2(ia)+wa2*scat2(ia+1)))*D%distance**2/1d23
			nf=(wl1*(wa1*iscat1(ia)+wa2*iscat1(ia+1))+wl2*(wa1*iscat2(ia)+wa2*iscat2(ia+1)))
			close(unit=40)
			return
		endif
		lam1=lam2
		scat1(1:nangle)=scat2(1:nangle)
		iscat1(1:nangle)=iscat2(1:nangle)
	enddo
	ReadMCScatt=0d0
	nf=0d0
	close(unit=40)
	return
	end

	subroutine MakeImage(image,tel,Rmax)
c	Rmax,IMDIM,nintegrate,width,SNR,Diam)
	use Parameters
	IMPLICIT NONE
	type(Telescope) tel
	type(RPhiImage) image
	integer IMDIM
c	parameter(IMDIM=500)
	real*8 w1,w2,dx,dy,gasdev,SNR
	real*8 wx,wy,wi,coswi,sinwi,cos2pa,sin2pa,Q,U
	real*8 Eu,Ex,Ey,EE
	
	real*8,allocatable :: im(:,:),imQ(:,:),imU(:,:),imV(:,:),imP(:,:)
	real*8 lam0,Rmax,i1,tot,ran2,max,psfdev,nphot,err
	integer ix,iy,i,j,k,i10,nintegrate,inphot
	integer il10000,il1000,il100,il10,jj,number_invalid
	integer ir10000,ir1000,ir100,ir10,ir1,jp1,jp2,ir2
	real*8 x,y,R,phi,flux,il1,min,wp1,wp2,ranR,ranPhi
	real*8 fluxQ,fluxU,fluxV,immin,poidev
	real*8 sint,cost,sin2t,cos2t,z,wr1,wr2
	character*500 filename
	type(photon) phot

	real*8 al11,ar11,i11,c11,p11,s11
	real*8 al12,ar12,i12,c12,p12,s12
	real*8 al21,ar21,i21,c21,p21,s21
	real*8 al22,ar22,i22,c22,p22,s22
	real*8,allocatable :: rrr(:,:)
	
	logical checktellertje

	IMDIM=tel%npixel
	nintegrate=tel%nint

	allocate(im(IMDIM,IMDIM))
	if(scat_how.eq.2.or.nplanets.gt.0) then
		allocate(imQ(IMDIM,IMDIM))
		allocate(imU(IMDIM,IMDIM))
		allocate(imV(IMDIM,IMDIM))
		allocate(imP(IMDIM,IMDIM))
	endif
	lam0=image%lam

c	wi=widthi*pi/180d0
c	coswi=cos(wi)
c	sinwi=sin(wi)
	
	write(*,'("Creating rectangular image at wavelength: ",f10.3)') lam0
	write(9,'("Creating rectangular image at wavelength: ",f10.3)') lam0

	im=0d0
	if(scat_how.eq.2) then
		imQ=0d0
		imU=0d0
		imV=0d0
	endif

	allocate(rrr(2,nintegrate))
	do k=1,nintegrate
		rrr(1,k)=ran2(idum)
		rrr(2,k)=ran2(idum)
	enddo

	call tellertje(1,100)
!$OMP PARALLEL IF(multicore)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(j,k,w2,ranR,ranPhi,R,phi,wp1,jp1,wp2,jp2,wr1,ir1,wr2,ir2,
!$OMP&  i11,i12,i21,i22,p11,p12,p21,p22,al11,al12,al21,al22,ar11,ar12,ar21,ar22,
!$OMP&  c11,c12,c21,c22,s11,s12,s21,s22,flux,fluxQ,fluxU,fluxV,x,ix,y,iy)
!$OMP& SHARED(image,im,imQ,imU,imV,D,nintegrate,Rmax,scat_how,IMDIM,rrr)

!$OMP DO
	do i=1,image%nr-1

	if(checktellertje(i+1,image%nr+1)) then
!$OMP CRITICAL
		write(*,'(".",$)')
		call flush(6)
		write(9,'(".",$)')
		call flush(9)
!$OMP END CRITICAL
	endif
	do j=1,image%nphi
		w2=2d0*pi*AU**2*(image%R(i+1)-image%R(i))/real(image%nphi)
		do k=1,nintegrate
			ranR=rrr(1,k)		!ran2(idum)
			ranPhi=rrr(2,k)		!ran2(idum)

			R=image%R(i)+(image%R(i+1)-image%R(i))*ranR
			phi=pi*(real(j-1)+ranPhi)/real(image%nphi)

			if(ranPhi.lt.0.5d0) then
				wp1=0.5d0-ranPhi
				jp1=j-1
				if(jp1.eq.0) jp1=1
				wp2=ranPhi+0.5d0
				jp2=j
			else
				wp1=ranPhi-0.5d0
				jp1=j+1
				if(jp1.eq.image%nphi+1) jp1=image%nphi
				wp2=1.5d0-ranPhi
				jp2=j
			endif
			
			wr1=(1d0-ranR)
			ir1=i
			wr2=ranR
			ir2=i+1
			
			if(scat_how.eq.2) then
				if(abs(image%R(ir2)-image%R(ir1)).gt.(Rmax/real(IMDIM))) then
				i11=image%image(ir1,jp1)-sqrt(image%imageQ(ir1,jp1)**2+image%imageU(ir1,jp1)**2+image%imageV(ir1,jp1)**2)
				i12=image%image(ir1,jp2)-sqrt(image%imageQ(ir1,jp2)**2+image%imageU(ir1,jp2)**2+image%imageV(ir1,jp2)**2)
				i21=image%image(ir2,jp1)-sqrt(image%imageQ(ir2,jp1)**2+image%imageU(ir2,jp1)**2+image%imageV(ir2,jp1)**2)
				i22=image%image(ir2,jp2)-sqrt(image%imageQ(ir2,jp2)**2+image%imageU(ir2,jp2)**2+image%imageV(ir2,jp2)**2)
				p11=sqrt(image%imageQ(ir1,jp1)**2+image%imageU(ir1,jp1)**2+image%imageV(ir1,jp1)**2)
				p12=sqrt(image%imageQ(ir1,jp2)**2+image%imageU(ir1,jp2)**2+image%imageV(ir1,jp2)**2)
				p21=sqrt(image%imageQ(ir2,jp1)**2+image%imageU(ir2,jp1)**2+image%imageV(ir2,jp1)**2)
				p22=sqrt(image%imageQ(ir2,jp2)**2+image%imageU(ir2,jp2)**2+image%imageV(ir2,jp2)**2)
				al11=sqrt((p11+image%imageQ(ir1,jp1))/2d0)
				al12=sqrt((p12+image%imageQ(ir1,jp2))/2d0)
				al21=sqrt((p21+image%imageQ(ir2,jp1))/2d0)
				al22=sqrt((p22+image%imageQ(ir2,jp2))/2d0)
				ar11=sqrt((p11-image%imageQ(ir1,jp1))/2d0)
				ar12=sqrt((p12-image%imageQ(ir1,jp2))/2d0)
				ar21=sqrt((p21-image%imageQ(ir2,jp1))/2d0)
				ar22=sqrt((p22-image%imageQ(ir2,jp2))/2d0)
				c11=image%imageU(ir1,jp1)/(2d0*al11*ar11)
				c12=image%imageU(ir1,jp2)/(2d0*al12*ar12)
				c21=image%imageU(ir2,jp1)/(2d0*al21*ar21)
				c22=image%imageU(ir2,jp2)/(2d0*al22*ar22)
				s11=image%imageV(ir1,jp1)/(2d0*al11*ar11)
				s12=image%imageV(ir1,jp2)/(2d0*al12*ar12)
				s21=image%imageV(ir2,jp1)/(2d0*al21*ar21)
				s22=image%imageV(ir2,jp2)/(2d0*al22*ar22)

				i11=log10(i11)
				i12=log10(i12)
				i21=log10(i21)
				i22=log10(i22)
				al11=log10(al11)
				al12=log10(al12)
				al21=log10(al21)
				al22=log10(al22)
				ar11=log10(ar11)
				ar12=log10(ar12)
				ar21=log10(ar21)
				ar22=log10(ar22)

				flux=(10d0**(wr1*(wp1*al11+wp2*al12)+wr2*(wp1*al21+wp2*al22)))**2
				flux=flux+(10d0**(wr1*(wp1*ar11+wp2*ar12)+wr2*(wp1*ar21+wp2*ar22)))**2
				flux=flux+10d0**(wr1*(wp1*i11+wp2*i12)+wr2*(wp1*i21+wp2*i22))
				flux=flux*w2*R

				fluxQ=(10d0**(wr1*(wp1*al11+wp2*al12)+wr2*(wp1*al21+wp2*al22)))**2
				fluxQ=fluxQ-(10d0**(wr1*(wp1*ar11+wp2*ar12)+wr2*(wp1*ar21+wp2*ar22)))**2
				fluxQ=fluxQ*w2*R

				fluxU=10d0**(wr1*(wp1*al11+wp2*al12)+wr2*(wp1*al21+wp2*al22))
				fluxU=fluxU*(10d0**(wr1*(wp1*ar11+wp2*ar12)+wr2*(wp1*ar21+wp2*ar22)))
				fluxU=fluxU*2d0*(wr1*(wp1*c11+wp2*c12)+wr2*(wp1*c21+wp2*c22))
				fluxU=fluxU*w2*R

				fluxV=10d0**(wr1*(wp1*al11+wp2*al12)+wr2*(wp1*al21+wp2*al22))
				fluxV=fluxV*(10d0**(wr1*(wp1*ar11+wp2*ar12)+wr2*(wp1*ar21+wp2*ar22)))
				fluxV=fluxV*2d0*(wr1*(wp1*s11+wp2*s12)+wr2*(wp1*s21+wp2*s22))
				fluxV=fluxV*w2*R

				if(number_invalid(flux)+number_invalid(fluxQ)+number_invalid(fluxU)+number_invalid(fluxV).ne.0) then
					flux=(wr1*(wp1*image%image(ir1,jp1)+wp2*image%image(ir1,jp2))
     &					+wr2*(wp1*image%image(ir2,jp1)+wp2*image%image(ir2,jp2)))*w2*R
					fluxQ=(wr1*(wp1*image%imageQ(ir1,jp1)+wp2*image%imageQ(ir1,jp2))
     &					+wr2*(wp1*image%imageQ(ir2,jp1)+wp2*image%imageQ(ir2,jp2)))*w2*R
					fluxU=(wr1*(wp1*image%imageU(ir1,jp1)+wp2*image%imageU(ir1,jp2))
     &					+wr2*(wp1*image%imageU(ir2,jp1)+wp2*image%imageU(ir2,jp2)))*w2*R
					fluxV=(wr1*(wp1*image%imageV(ir1,jp1)+wp2*image%imageV(ir1,jp2))
     &					+wr2*(wp1*image%imageV(ir2,jp1)+wp2*image%imageV(ir2,jp2)))*w2*R
				endif

				else
				flux=(wr1*(wp1*image%image(ir1,jp1)+wp2*image%image(ir1,jp2))
     &				+wr2*(wp1*image%image(ir2,jp1)+wp2*image%image(ir2,jp2)))*w2*R
				fluxQ=(wr1*(wp1*image%imageQ(ir1,jp1)+wp2*image%imageQ(ir1,jp2))
     &				+wr2*(wp1*image%imageQ(ir2,jp1)+wp2*image%imageQ(ir2,jp2)))*w2*R
				fluxU=(wr1*(wp1*image%imageU(ir1,jp1)+wp2*image%imageU(ir1,jp2))
     &				+wr2*(wp1*image%imageU(ir2,jp1)+wp2*image%imageU(ir2,jp2)))*w2*R
				fluxV=(wr1*(wp1*image%imageV(ir1,jp1)+wp2*image%imageV(ir1,jp2))
     &				+wr2*(wp1*image%imageV(ir2,jp1)+wp2*image%imageV(ir2,jp2)))*w2*R
				endif
			else
				flux=(wr1*(wp1*log10(image%image(ir1,jp1))+wp2*log10(image%image(ir1,jp2)))
     &				+wr2*(wp1*log10(image%image(ir2,jp1))+wp2*log10(image%image(ir2,jp2))))
     			flux=(10d0**flux)*w2*R
			endif

c			wx=gasdev(idum)*widthy
c			wy=gasdev(idum)*widthx
			
			x=-R*cos(phi+D%PA*pi/180d0)!+wx*coswi+wy*sinwi
			x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1d0
			ix=x
			y=R*sin(phi+D%PA*pi/180d0)!+wy*coswi+wx*sinwi
			y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1d0
			iy=y
			if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
				im(ix,iy)=im(ix,iy)+flux/real(2*nintegrate)
				if(scat_how.eq.2) then
					imQ(ix,iy)=imQ(ix,iy)+fluxQ/real(2*nintegrate)
					imU(ix,iy)=imU(ix,iy)+fluxU/real(2*nintegrate)
					imV(ix,iy)=imV(ix,iy)+fluxV/real(2*nintegrate)
				endif
			endif

c			wx=gasdev(idum)*widthy
c			wy=gasdev(idum)*widthx

			x=-R*cos(-phi+D%PA*pi/180d0)!+wx*coswi+wy*sinwi
			x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1d0
			ix=x
			y=R*sin(-phi+D%PA*pi/180d0)!+wy*coswi+wx*sinwi
			y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1d0
			iy=y
			if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
				im(ix,iy)=im(ix,iy)+flux/real(2*nintegrate)
				if(scat_how.eq.2) then
					imQ(ix,iy)=imQ(ix,iy)+fluxQ/real(2*nintegrate)
					imU(ix,iy)=imU(ix,iy)-fluxU/real(2*nintegrate)
					imV(ix,iy)=imV(ix,iy)-fluxV/real(2*nintegrate)
				endif
			endif
		enddo
	enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(100,100)
	
	deallocate(rrr)


	call AddPlanet(image,im,imQ,imU,imV,Rmax,IMDIM)

	if(tel%psffile.ne.' ') then
		call ConvolutionFILE(im,imQ,imU,imV,IMDIM,tel%psffile)
	else if(tel%D.gt.0d0.or.tel%width.gt.0d0) then
		call Convolution(im,imQ,imU,imV,IMDIM,lam0,tel%D,tel%D2,tel%spider
     &			,tel%mask,tel%wmask,tel%owa,tel%strehl,Rmax*2d0*parsec/D%distance
     &			,tel%width,tel%snoise)
	endif
	
	do i=1,IMDIM
	do j=1,IMDIM
		im(i,j)=im(i,j)*1d23*image%zscale/D%distance**2
		if(scat_how.eq.2.and.useobspol) then
			imQ(i,j)=imQ(i,j)*1d23*image%zscale/D%distance**2
			imU(i,j)=imU(i,j)*1d23*image%zscale/D%distance**2
			imV(i,j)=imV(i,j)*1d23*image%zscale/D%distance**2
			phot%Q=imQ(i,j)
			phot%U=imU(i,j)
			if(D%PA.ne.0d0.and.im(i,j).ne.0d0) then
				x=1d0
				y=0d0
				z=0d0
				phot%vx=0d0
				phot%vy=0d0
				phot%vz=1d0
				phot%Sx=cos(-D%PA*pi/180d0)
				phot%Sy=sin(-D%PA*pi/180d0)
				phot%Sz=0d0
				phot%E=im(i,j)
				call RotateStokes(phot,x,y,z)
			endif
			imQ(i,j)=phot%Q
			imU(i,j)=phot%U
		endif
	enddo
	enddo

	if(tel%texp.gt.0d0.and.tel%scaletype.eq.2.and.tel%D.gt.0d0) then
	do i=1,IMDIM
	do j=1,IMDIM
		if(scat_how.ne.2) then
			nphot=im(i,j)*1d-3*(Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2)
   		  	nphot=nphot*tel%texp
			if(nphot.lt.100d0) then
				inphot=poidev(nphot,idum)
				nphot=real(inphot)/tel%texp
				im(i,j)=nphot*1d3/((Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2))
			else
				im(i,j)=im(i,j)+im(i,j)*gasdev(idum)/sqrt(nphot)
			endif
			im(i,j)=im(i,j)+gasdev(idum)*tel%RON
		else
			Eu=im(i,j)
			im(i,j)=0d0
		
			Ex=(Eu+imQ(i,j))/2d0
			nphot=Ex*1d-3*(Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2)
			nphot=nphot*tel%texp/4d0
			if(nphot.lt.100d0) then
				inphot=poidev(nphot,idum)
				nphot=real(inphot)*4d0/tel%texp
				Ex=nphot*1d3/((Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2))
			else
				Ex=Ex+Ex*gasdev(idum)/sqrt(nphot)
			endif
			Ex=Ex+gasdev(idum)*tel%RON
			if(Ex.lt.0d0) Ex=0d0
			Ey=(Eu-imQ(i,j))/2d0
			nphot=Ey*1d-3*(Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2)
			nphot=nphot*tel%texp/4d0
			if(nphot.lt.100d0) then
				inphot=poidev(nphot,idum)
				nphot=real(inphot)*4d0/tel%texp
				Ey=nphot*1d3/((Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2))
			else
				Ey=Ey+Ey*gasdev(idum)/sqrt(nphot)
			endif
			Ey=Ey+gasdev(idum)*tel%RON
			if(Ey.lt.0d0) Ey=0d0
			im(i,j)=im(i,j)+(Ex+Ey)/2d0
			imQ(i,j)=Ex-Ey

			Ex=(Eu+imU(i,j))/2d0
			nphot=Ex*1d-3*(Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2)
			nphot=nphot*tel%texp/4d0
			if(nphot.lt.100d0) then
				inphot=poidev(nphot,idum)
				nphot=real(inphot)*4d0/tel%texp
				Ex=nphot*1d3/((Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2))
			else
				Ex=Ex+Ex*gasdev(idum)/sqrt(nphot)
			endif
			Ex=Ex+gasdev(idum)*tel%RON
			if(Ex.lt.0d0) Ex=0d0
			Ey=(Eu-imU(i,j))/2d0
			nphot=Ey*1d-3*(Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2)
			nphot=nphot*tel%texp/4d0
			if(nphot.lt.100d0) then
				inphot=poidev(nphot,idum)
				nphot=real(inphot)*4d0/tel%texp
				Ey=nphot*1d3/((Rmax*2d0*parsec/D%distance)**2*
     &			1.509e7*tel%dlam*pi*((tel%D/2d0)**2-(tel%D2/2d0)**2)/(lam0*real(IMDIM)**2))
			else
				Ey=Ey+Ey*gasdev(idum)/sqrt(nphot)
			endif
			Ey=Ey+gasdev(idum)*tel%RON
			if(Ey.lt.0d0) Ey=0d0
			im(i,j)=im(i,j)+(Ex+Ey)/2d0
			imU(i,j)=Ex-Ey
		endif
	enddo
	enddo
	endif
	if(scat_how.eq.2.and.tel%Ptelescope.ne.0d0) then
c add telescope polarization
		imQ=imQ+tel%Ptelescope*cos(tel%APtelescope*pi/180d0)*im
		imU=imU+tel%Ptelescope*sin(tel%APtelescope*pi/180d0)*im
	endif

	call FITSImage(image,im,IMDIM,Rmax,'AI',tel%flag)
	if(scat_how.eq.2.and.useobspol) then
		call FITSImage(image,imQ,IMDIM,Rmax,'AQ',tel%flag)
		call FITSImage(image,imU,IMDIM,Rmax,'AU',tel%flag)
		call FITSImage(image,imV,IMDIM,Rmax,'AV',tel%flag)
		imP=imQ/im
		call FITSImage(image,imP,IMDIM,Rmax,'RQ',tel%flag)
		imP=imU/im
		call FITSImage(image,imP,IMDIM,Rmax,'RU',tel%flag)
		imP=sqrt((imQ/im)**2+(imU/im)**2)
		call FITSImage(image,imP,IMDIM,Rmax,'LP',tel%flag)
		imP=imV/im
		call FITSImage(image,imP,IMDIM,Rmax,'CP',tel%flag)
		do i=1,IMDIM
		do j=1,IMDIM
			if(i.ne.((IMDIM+1)/2).or.j.ne.((IMDIM+1)/2)) then
				x=real(i)-real(IMDIM+1)/2d0
				y=real(j)-real(IMDIM+1)/2d0
				R=sqrt(x**2+y**2)
				x=x/R
				y=y/R	
				z=0d0
				phot%vx=0d0
				phot%vy=0d0
				phot%vz=1d0
				phot%Sx=1d0
				phot%Sy=0.00001d0
				phot%Sz=0d0
				phot%Q=imQ(i,j)
				phot%U=imU(i,j)
				phot%E=im(i,j)
				call RotateStokes(phot,x,y,z)
				imP(i,j)=phot%Q
			endif
		enddo
		enddo
		imP=imP/im
		call FITSImage(image,imP,IMDIM,Rmax,'DQ',tel%flag)
		do i=1,IMDIM
		do j=1,IMDIM
			if(i.ne.((IMDIM+1)/2).or.j.ne.((IMDIM+1)/2)) then
				x=real(i)-real(IMDIM+1)/2d0
				y=real(j)-real(IMDIM+1)/2d0
				R=sqrt(x**2+y**2)
				x=x/R
				y=y/R
				z=0d0
				phot%vx=0d0
				phot%vy=0d0
				phot%vz=1d0
				phot%Sx=1d0
				phot%Sy=0.00001d0
				phot%Sz=0d0
				phot%Q=imQ(i,j)
				phot%U=imU(i,j)
				phot%E=im(i,j)
				call RotateStokes(phot,x,y,z)
				imP(i,j)=phot%U
			endif
		enddo
		enddo
		imP=imP/im
		call FITSImage(image,imP,IMDIM,Rmax,'DU',tel%flag)
		imP=sqrt(imQ**2+imU**2)
		call FITSImage(image,imP,IMDIM,Rmax,'QU',tel%flag)
		call ExPoImage(image,im,imQ,imU,IMDIM,Rmax,'EX',tel%flag)
		if(tel%iprad.gt.0d0) then
			fluxQ=0d0
			fluxU=0d0
			flux=0d0
			do i=1,IMDIM
			do j=1,IMDIM
				x=real(i)-real(IMDIM+1)/2d0
				y=real(j)-real(IMDIM+1)/2d0
				R=sqrt(x**2+y**2)
				if(R.le.tel%iprad) then
					fluxQ=fluxQ+imQ(i,j)
					fluxU=fluxU+imU(i,j)
					flux=flux+im(i,j)
				endif
			enddo
			enddo
			print*,fluxQ/flux
			print*,fluxU/flux
			call ExPoImage(image,im,imQ-im*fluxQ/flux,imU-im*fluxU/flux,IMDIM,Rmax,'CX',tel%flag)
		endif
	endif
	deallocate(im)
	if(scat_how.eq.2.or.nplanets.gt.0) then
		deallocate(imQ)
		deallocate(imU)
		deallocate(imV)
		deallocate(imP)
	endif

	i10=(180d0*image%angle/pi)/10d0
	i1=(180d0*image%angle/pi)-10d0*i10
	if(i1.ge.9.95d0) then
		i1=i1-9.95d0
		i10=i10+1
	endif
	
	il10000=lam0/10000d0
	il1000=(lam0-10000d0*il10000)/1000d0
	il100=(lam0-1000d0*il1000-10000d0*il10000)/100d0
	il10=(lam0-100d0*il100-1000d0*il1000-10000d0*il10000)/10d0
	il1=lam0-10d0*il10-100d0*il100-1000d0*il1000-10000d0*il10000

	ir10000=Rmax/10000d0
	ir1000=(Rmax-10000d0*ir10000)/1000d0
	ir100=(Rmax-1000d0*ir1000-10000d0*ir10000)/100d0
	ir10=(Rmax-100d0*ir100-1000d0*ir1000-10000d0*ir10000)/10d0
	ir1=Rmax-10d0*ir10-100d0*ir100-1000d0*ir1000-10000d0*ir10000

	write(filename,'(a,"RPhiImage_i",i1,f3.1,"_l",i1,i1,i1,i1,f4.2,"_z",i1,i1,i1,i1,i1,".dat")')
     & outdir(1:len_trim(outdir)),i10,i1,il10000,il1000,il100,il10,il1,ir10000,ir1000,ir100,ir10,ir1
c	open(unit=50,file=filename,RECL=1000)
c	do i=1,image%nr
c	do j=1,image%nphi
c		write(50,*) image%r(i)*sin(image%phi(j))*image%rscale,image%r(i)*cos(image%phi(j))
c     &					*image%rscale,image%image(i,j)*image%zscale
c	enddo
c	enddo
c	close(unit=50)

	return
	end


	subroutine AddPlanet(image,im,imQ,imU,imV,Rmax,IMDIM)
	use Parameters
	IMPLICIT NONE
	type(RPhiImage) image
	integer IMDIM
	real*8 w1,w2,widthx,widthy,widthi,dx,dy,gasdev,SNR
	real*8 wx,wy,wi,coswi,sinwi,cos2pa,sin2pa,Q,U

	real*8 im(IMDIM,IMDIM),imQ(IMDIM,IMDIM),imU(IMDIM,IMDIM),imV(IMDIM,IMDIM)
	real*8 lam0,Rmax,i1,tot,ran2,max,psfdev
	integer ix,iy,i,j,k,i10,nintegrate,iplanets
	integer il10000,il1000,il100,il10,jj
	integer ir10000,ir1000,ir100,ir10,ir1,jp1,jp2
	real*8 x,y,R,phi,flux,il1,min,wp1,wp2,ranR,ranPhi
	real*8 fluxQ,fluxU,fluxV,immin,Tflux,Planck
	real*8 sint,cost,sin2t,cos2t,z,tau,fact,theta
	character*500 filename
	type(photon) phot

c	wi=widthi*pi/180d0
c	coswi=cos(wi)
c	sinwi=sin(wi)
	
	do iplanets=1,nplanets
	Tflux=Planck(Planets(iplanets)%T,image%lam)*pi*Planets(iplanets)%R**2

	R=Planets(iplanets)%d
	phi=Planets(iplanets)%phi*pi/180d0
	theta=Planets(iplanets)%theta

c	x=R*cos(phi)*cos(image%angle)
c	y=R*sin(phi)

	x=R*cos(phi)*sin(theta)*cos(image%angle)-R*cos(theta)*sin(image%angle)
	y=R*sin(phi)*sin(theta)

	R=sqrt(x**2+y**2)
	phi=acos(x/R)
	if(y.lt.0d0) phi=2d0*pi-phi

	do i=1,image%nr-1
		if(R.gt.image%R(i).and.R.le.image%R(i+1)) goto 1
	enddo
1	continue
	if(phi.lt.pi) then
		j=phi*real(image%nphi)/pi+1
	else
		j=(2d0*pi-phi)*real(image%nphi)/pi+1
	endif
	if(j.lt.1) j=1
	if(j.gt.image%nphi) j=image%nphi


	tau=0d0
	do k=1,image%p(i,j)%n
		if(image%p(i,j)%k(k).eq.1) then
			tau=tau+image%p(i,j)%v(k)*C(image%p(i,j)%i(k),image%p(i,j)%j(k))%dens
     &				*C(image%p(i,j)%i(k),image%p(i,j)%j(k))%Kext*AU
		endif
	enddo
	fact=exp(-tau)

	write(*,'("Adding planet: ",a)') Planets(iplanets)%name(1:len_trim(Planets(iplanets)%name))
	write(9,'("Adding planet: ",a)') Planets(iplanets)%name(1:len_trim(Planets(iplanets)%name))

c	do k=1,nintegrate*1000
c		call tellertje(k,nintegrate*1000)
c		wx=gasdev(idum)*widthy
c		wy=gasdev(idum)*widthx

		x=-R*cos(phi+D%PA*pi/180d0)!+wx*coswi+wy*sinwi
		x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1
		ix=x
		y=R*sin(phi+D%PA*pi/180d0)!+wy*coswi+wx*sinwi
		y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1
		iy=y
		if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
			im(ix,iy)=im(ix,iy)+((Planets(iplanets)%E+Tflux)*fact)!/real(nintegrate*1000)
			if(scat_how.eq.2) then
				imQ(ix,iy)=imQ(ix,iy)+(Planets(iplanets)%Q*fact)!/real(nintegrate*1000)
				imU(ix,iy)=imU(ix,iy)+(Planets(iplanets)%U*fact)!/real(nintegrate*1000)
				imV(ix,iy)=imV(ix,iy)+(Planets(iplanets)%V*fact)!/real(nintegrate*1000)
			endif
			tot=((Planets(iplanets)%E+Tflux)*fact)
		endif
c	enddo
	enddo

	print*,tot/sum(im(1:IMDIM,1:IMDIM))

	return
	end
	


	subroutine RImage(image,tel)
	use Parameters
	IMPLICIT NONE
	type(RPhiImage) image
	character*500 filename
	real*8 flux,w1,w2,il1,i1,fluxR(image%nr),scaling,P,Pflux,Qflux,Uflux,PSflux
	integer i,k,i10,il10,il100,il1000,il10000,IMDIM
	type(Telescope) tel

	scaling=(image%zscale*1d23/D%distance**2)*(tel%fov(1)**2*AU**2)/(tel%npixel**2)

	i10=(180d0*image%angle/pi)/10d0
	i1=(180d0*image%angle/pi)-10d0*i10
	if(i1.ge.9.95d0) then
		i1=i1-9.95d0
		i10=i10+1
	endif

	il10000=image%lam/10000d0
	il1000=(image%lam-10000d0*il10000)/1000d0
	il100=(image%lam-1000d0*il1000-10000d0*il10000)/100d0
	il10=(image%lam-100d0*il100-1000d0*il1000-10000d0*il10000)/10d0
	il1=image%lam-10d0*il10-100d0*il100-1000d0*il1000-10000d0*il10000
	write(filename,'(a,"RImage_i",i1,f3.1,"_l",i1,i1,i1,i1,f4.2,".dat")')
     & outdir(1:len_trim(outdir)),i10,i1,il10000,il1000,il100,il10,il1
	open(unit=50,file=filename,RECL=1000)


	flux=0d0
	Pflux=0d0
	Qflux=0d0
	Uflux=0d0
	do i=1,image%nr-1
	PSflux=0d0
	do k=1,image%nPhi
		w1=2d0*pi*abs(image%R(i))*AU**2/real(image%nPhi)
		w2=2d0*pi*abs(image%R(i+1))*AU**2/real(image%nPhi)
		flux=flux+(image%R(i+1)-image%R(i))*
     &		(w1*image%image(i,k)+w2*image%image(i+1,k))/2d0
		if(scat_how.eq.2) then
			Qflux=Qflux+(image%R(i+1)-image%R(i))*
     &		(w1*image%imageQ(i,k)+w2*image%imageQ(i+1,k))/2d0
			Uflux=Uflux+(image%R(i+1)-image%R(i))*
     &		(w1*image%imageU(i,k)+w2*image%imageU(i+1,k))/2d0
			PSflux=PSflux+(image%R(i+1)-image%R(i))*
     &		sqrt(((w1*image%imageU(i,k)+w2*image%imageU(i+1,k))/2d0)**2+
     &			 ((w1*image%imageQ(i,k)+w2*image%imageQ(i+1,k))/2d0)**2)
		endif
		PSflux=PSflux/(2d0*pi*AU**2*(image%R(i+1)**2-image%R(i)**2))
	enddo
	if(scat_how.ne.2) then
		write(50,*) image%R(i)*image%rscale,flux*1d23/D%distance**2,scaling*image%image(i,image%nPhi/2)
	else
		write(50,*) image%R(i)*image%rscale,flux*1d23/D%distance**2,scaling*image%image(i,image%nPhi/2)
     &				,Qflux*1d23/D%distance**2,Uflux*1d23/D%distance**2,scaling*PSflux
	endif
	enddo
	if(scat_how.ne.2) then
		write(50,*) image%R(image%nr)*image%rscale,flux*1d23/D%distance**2,scaling*image%image(image%nr,image%nPhi/2)
	else
		P=sqrt(image%imageQ(image%nr,image%nPhi/2)**2+image%imageU(image%nr,image%nPhi/2)**2)
		write(50,*) image%R(image%nr)*image%rscale,flux*1d23/D%distance**2,scaling*image%image(image%nr,image%nPhi/2),scaling*P
	endif

	close(unit=50)

	return
	end



	subroutine FITSImage(image,im,IMDIM,Rmax,add,flag)
	use Parameters
	IMPLICIT NONE
	integer IMDIM
	type(RPhiImage) image
	character*500 filename,cmd
	character*2 add
	character*20 flag
	real*8 flux,w1,w2,il1,i1,ir1,im(IMDIM,IMDIM)
	real*8 Rmax,cornerpix,pixscale
	real*8,allocatable :: array(:,:)
	integer i,k,i10,il10,il100,il1000,il10000
	integer ir10,ir100,ir1000,ir10000

      integer status,unit,blocksize,bitpix,naxis,naxes(2)
      integer j,group,fpixel,nelements
      logical simple,extend,truefalse

	allocate(array(IMDIM,IMDIM))


	i10=(180d0*image%angle/pi)/10d0
	i1=(180d0*image%angle/pi)-10d0*i10
	if(i1.ge.9.95d0) then
		i1=i1-9.95d0
		i10=i10+1
	endif

	il10000=image%lam/10000d0
	il1000=(image%lam-10000d0*il10000)/1000d0
	il100=(image%lam-1000d0*il1000-10000d0*il10000)/100d0
	il10=(image%lam-100d0*il100-1000d0*il1000-10000d0*il10000)/10d0
	il1=image%lam-10d0*il10-100d0*il100-1000d0*il1000-10000d0*il10000

	ir10000=2d0*Rmax/10000d0
	ir1000=(2d0*Rmax-10000d0*ir10000)/1000d0
	ir100=(2d0*Rmax-1000d0*ir1000-10000d0*ir10000)/100d0
	ir10=(2d0*Rmax-100d0*ir100-1000d0*ir1000-10000d0*ir10000)/10d0
	ir1=2d0*Rmax-10d0*ir10-100d0*ir100-1000d0*ir1000-10000d0*ir10000

	if(ir1.ge.9.95d0) then
		ir10=ir10+1
		ir1=ir1-9.95d0
	endif
	if(ir10.ge.10) then
		ir100=ir100+1
		ir10=ir10-10
	endif
	if(ir100.ge.10) then
		ir1000=ir1000+1
		ir100=ir100-10
	endif
	if(ir1000.ge.10) then
		ir10000=ir10000+1
		ir1000=ir1000-10
	endif

	write(filename,'(a,"Image",a2,"_i",i1,f3.1,"_l",i1,i1,i1,i1,f4.2,"_fov",i1,i1,i1,i1,f3.1,a,".fits")')
     & outdir(1:len_trim(outdir)),add,i10,i1,il10000,il1000,il100,il10,il1,ir10000,ir1000,ir100,ir10,ir1,flag(1:len_trim(flag))

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
      naxes(1)=IMDIM
      naxes(2)=IMDIM
      extend=.true.

C     write the required header keywords
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do i=1,IMDIM
	do j=1,IMDIM
		array(j,i)=im(i,j)
	enddo
	enddo

C     write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)
      call ftpprd(unit,group,fpixel,nelements,array,status)

	if(image%scaletype.eq.1) then
C     write another optional keyword to the header
	    call ftpkys(unit,'CTYPE1','X [AU]','X-coordinate',status)
      	call ftpkys(unit,'CTYPE2','Y [AU]','Y-coordinate',status)
	else if(image%scaletype.eq.2) then
      	call ftpkys(unit,'CTYPE1','X [arcsec]','X-coordinate',status)
      	call ftpkys(unit,'CTYPE2','Y [arcsec]','Y-coordinate',status)
	endif
      call ftpkyd(unit,'CRPIX1',0.5d0,10,'Corner of the image',status)
      call ftpkyd(unit,'CRPIX2',0.5d0,10,'Corner of the image',status)

		cornerpix=-Rmax*image%rscale
      call ftpkyd(unit,'CRVAL1',cornerpix,10,'Center of the image',status)
      call ftpkyd(unit,'CRVAL2',cornerpix,10,'Center of the image',status)

		pixscale=2d0*Rmax*image%rscale/real(IMDIM)
	if(image%scaletype.eq.1) then
      	call ftpkyd(unit,'CDELT1',pixscale,10,'Pixel scale (in AU)',status)
      	call ftpkyd(unit,'CDELT2',pixscale,10,'Pixel scale (in AU)',status)
	else if(image%scaletype.eq.2) then
      	call ftpkyd(unit,'CDELT1',pixscale,10,'Pixel scale (in arcsec)',status)
      	call ftpkyd(unit,'CDELT2',pixscale,10,'Pixel scale (in arcsec)',status)
	endif
	

C     close the file and free the unit number
      call ftclos(unit, status)
      call ftfiou(unit, status)

	deallocate(array)

	return
	end





	subroutine ExPoImage(image,im,imQ,imU,IMDIM,Rmax,add,flag)
	use Parameters
	IMPLICIT NONE
	integer IMDIM
	type(RPhiImage) image
	character*500 filename,cmd
	character*2 add
	character*20 flag
	real*8 flux,w1,w2,il1,i1,ir1,im(IMDIM,IMDIM),imQ(IMDIM,IMDIM),imU(IMDIM,IMDIM)
	real*8 Rmax,cornerpix,pixscale
	real*8,allocatable :: array(:,:,:)
	integer i,k,i10,il10,il100,il1000,il10000
	integer ir10,ir100,ir1000,ir10000

      integer status,unit,blocksize,bitpix,naxis,naxes(3)
      integer j,group,fpixel,nelements
      logical simple,extend,truefalse

	allocate(array(IMDIM,IMDIM,4))


	i10=(180d0*image%angle/pi)/10d0
	i1=(180d0*image%angle/pi)-10d0*i10
	if(i1.ge.9.95d0) then
		i1=i1-9.95d0
		i10=i10+1
	endif

	il10000=image%lam/10000d0
	il1000=(image%lam-10000d0*il10000)/1000d0
	il100=(image%lam-1000d0*il1000-10000d0*il10000)/100d0
	il10=(image%lam-100d0*il100-1000d0*il1000-10000d0*il10000)/10d0
	il1=image%lam-10d0*il10-100d0*il100-1000d0*il1000-10000d0*il10000

	ir10000=2d0*Rmax/10000d0
	ir1000=(2d0*Rmax-10000d0*ir10000)/1000d0
	ir100=(2d0*Rmax-1000d0*ir1000-10000d0*ir10000)/100d0
	ir10=(2d0*Rmax-100d0*ir100-1000d0*ir1000-10000d0*ir10000)/10d0
	ir1=2d0*Rmax-10d0*ir10-100d0*ir100-1000d0*ir1000-10000d0*ir10000

	if(ir1.ge.9.95d0) then
		ir10=ir10+1
		ir1=ir1-9.95d0
	endif
	if(ir10.ge.10) then
		ir100=ir100+1
		ir10=ir10-10
	endif
	if(ir100.ge.10) then
		ir1000=ir1000+1
		ir100=ir100-10
	endif
	if(ir1000.ge.10) then
		ir10000=ir10000+1
		ir1000=ir1000-10
	endif

	write(filename,'(a,"Image",a2,"_i",i1,f3.1,"_l",i1,i1,i1,i1,f4.2,"_fov",i1,i1,i1,i1,f3.1,a,".fits")')
     & outdir(1:len_trim(outdir)),add,i10,i1,il10000,il1000,il100,il10,il1,ir10000,ir1000,ir100,ir10,ir1,flag(1:len_trim(flag))
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
      naxis=3
      naxes(1)=IMDIM
      naxes(2)=IMDIM
	  naxes(3)=4
      extend=.true.

C     write the required header keywords
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do i=1,IMDIM
	do j=1,IMDIM
		array(j,i,1)=im(i,j)
		array(j,i,2)=imQ(i,j)
		array(j,i,3)=imU(i,j)
		array(j,i,4)=sqrt(imQ(i,j)**2+imU(i,j)**2)
	enddo
	enddo

C     write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)*naxes(3)
      call ftpprd(unit,group,fpixel,nelements,array,status)

	if(image%scaletype.eq.1) then
C     write another optional keyword to the header
	    call ftpkys(unit,'CTYPE1','X [AU]','X-coordinate',status)
      	call ftpkys(unit,'CTYPE2','Y [AU]','Y-coordinate',status)
	else if(image%scaletype.eq.2) then
      	call ftpkys(unit,'CTYPE1','X [arcsec]','X-coordinate',status)
      	call ftpkys(unit,'CTYPE2','Y [arcsec]','Y-coordinate',status)
	endif
      call ftpkyd(unit,'CRPIX1',0.5d0,10,'Corner of the image',status)
      call ftpkyd(unit,'CRPIX2',0.5d0,10,'Corner of the image',status)

		cornerpix=-Rmax*image%rscale
      call ftpkyd(unit,'CRVAL1',cornerpix,10,'Center of the image',status)
      call ftpkyd(unit,'CRVAL2',cornerpix,10,'Center of the image',status)

		pixscale=2d0*Rmax*image%rscale/real(IMDIM)
	if(image%scaletype.eq.1) then
      	call ftpkyd(unit,'CDELT1',pixscale,10,'Pixel scale (in AU)',status)
      	call ftpkyd(unit,'CDELT2',pixscale,10,'Pixel scale (in AU)',status)
	else if(image%scaletype.eq.2) then
      	call ftpkyd(unit,'CDELT1',pixscale,10,'Pixel scale (in arcsec)',status)
      	call ftpkyd(unit,'CDELT2',pixscale,10,'Pixel scale (in arcsec)',status)
	endif
	

C     close the file and free the unit number
      call ftclos(unit, status)
      call ftfiou(unit, status)

	deallocate(array)

	return
	end


	subroutine MakeSimpleImage(image,tel,Rmax,addedname)
c	Rmax,IMDIM,nintegrate,width,SNR,Diam)
	use Parameters
	IMPLICIT NONE
	type(Telescope) tel
	type(RPhiImage) image
	integer IMDIM
c	parameter(IMDIM=500)
	real*8 w1,w2,dx,dy,gasdev,SNR
	real*8 wx,wy,wi,coswi,sinwi,cos2pa,sin2pa,Q,U
	real*8 Eu,Ex,Ey,EE
	
	character*2,optional :: addedname

	real*8,allocatable :: im(:,:),imQ(:,:),imU(:,:),imV(:,:),imP(:,:)
	real*8 lam0,Rmax,i1,tot,ran2,max,psfdev,nphot,err
	integer ix,iy,i,j,k,i10,nintegrate,inphot
	integer il10000,il1000,il100,il10,jj,number_invalid
	integer ir10000,ir1000,ir100,ir10,ir1,jp1,jp2,ir2
	real*8 x,y,R,phi,flux,il1,min,wp1,wp2,ranR,ranPhi
	real*8 fluxQ,fluxU,fluxV,immin,poidev
	real*8 sint,cost,sin2t,cos2t,z,wr1,wr2
	character*500 filename
	type(photon) phot
	integer, allocatable :: count(:,:)

	real*8 al11,ar11,i11,c11,p11,s11
	real*8 al12,ar12,i12,c12,p12,s12
	real*8 al21,ar21,i21,c21,p21,s21
	real*8 al22,ar22,i22,c22,p22,s22
	
	logical checktellertje

	IMDIM=tel%npixel
	nintegrate=tel%nint

	allocate(im(IMDIM,IMDIM))
	allocate(count(IMDIM,IMDIM))
	lam0=image%lam

c	wi=widthi*pi/180d0
c	coswi=cos(wi)
c	sinwi=sin(wi)
	
	write(*,'("Creating rectangular image at wavelength: ",f10.3)') lam0
	write(9,'("Creating rectangular image at wavelength: ",f10.3)') lam0

	im=0d0
	count=0

	do i=1,image%nr-1
	call tellertje(i,image%nr-1)
	do j=1,image%nphi
		do k=1,nintegrate
			ranR=ran2(idum)
			ranPhi=ran2(idum)

			R=image%R(i)+(image%R(i+1)-image%R(i))*ranR
			phi=pi*(real(j-1)+ranPhi)/real(image%nphi)
			
			flux=image%image(i,j)

			x=-R*cos(phi)
			x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1
			ix=x
			y=R*sin(phi)	
			y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1
			iy=y
			if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
				im(ix,iy)=im(ix,iy)+flux
				count(ix,iy)=count(ix,iy)+1
				y=-R*sin(phi)	
				y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1
				iy=y
				im(ix,iy)=im(ix,iy)+flux
				count(ix,iy)=count(ix,iy)+1
			endif
		enddo
	enddo
	enddo
	im=im/real(count)

	call FITSImage(image,im,IMDIM,Rmax,addedname,tel%flag)

	deallocate(im)
	deallocate(count)

	return
	end



	subroutine MakeDensImage(image,Rmax,IMDIM,nintegrate)
	use Parameters
	IMPLICIT NONE
	type(RPhiImage) image
	integer IMDIM
	real*8 im(IMDIM,IMDIM),w1,w2,width,dx,dy,gasdev
	real*8 lam0,Rmax,i1,tot,ran2
	integer ix,iy,i,j,k,i10,nintegrate
	integer il10000,il1000,il100,il10,jj,ii
	integer ir10000,ir1000,ir100,ir10,ir1,jp1,jp2
	real*8 x,y,R,phi,flux,il1,min,wp1,wp2,ranR,ranPhi
	real*8 fluxQ,fluxU,fluxV,immin
	character*500 filename
	logical hit(IMDIM,IMDIM)
	
	write(*,'("Creating density image")')
	write(9,'("Creating density image")')
	
	hit=.false.

	image%nr=D%nR
	image%nphi=(D%nTheta-1)*2d0+1
	if(.not.allocated(image%image)) then
		allocate(image%image(image%nr,image%nphi))
		allocate(image%Phi(image%nphi))
		allocate(image%R(image%nr))
	endif
	
	do i=1,image%nr
		image%R(i)=D%R(i)
	enddo
	do i=1,D%nTheta-1
		image%Phi(i)=D%thet(i)
		image%Phi(image%nphi-i+1)=pi-D%thet(i)
	enddo
	image%Phi(D%nTheta)=pi/2d0
	image%nphi=image%nphi-1

	do ii=1,ngrains

	do i=1,image%nr
	do j=1,D%nTheta-1
		image%image(i,j)=C(i,j)%mass*C(i,j)%w(ii)
		image%image(i,image%nphi+1-j)=C(i,j)%mass*C(i,j)%w(ii)
	enddo
	enddo
	
	im=0d0

	do i=1,image%nr-1
	call tellertje(i,image%nr-1)
	do j=1,image%nphi
		w2=2d0*pi*AU**2*(image%R(i+1)-image%R(i))/real(image%nphi)
		do k=1,nintegrate
			ranR=ran2(idum)
			ranPhi=ran2(idum)

			R=image%R(i)+(image%R(i+1)-image%R(i))*ranR
			phi=pi*(real(j-1)+ranPhi)/real(image%nphi)
			if(j.ne.image%nphi) then
				phi=image%Phi(j)+ranPhi*(image%Phi(j+1)-image%Phi(j))
			else
				phi=image%Phi(j)+ranPhi*(pi-image%Phi(j))
			endif
			
			flux=image%image(i,j)

			x=-R*cos(phi)
			x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1
			ix=x
			y=R*sin(phi)	
			y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1
			iy=y
			if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
				im(ix,iy)=im(ix,iy)+flux/real(2*nintegrate)
				hit(ix,iy)=.true.
				y=-R*sin(phi)	
				y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1
				iy=y
				im(ix,iy)=im(ix,iy)+flux/real(2*nintegrate)
				hit(ix,iy)=.true.
			endif
		enddo
	enddo
	enddo

	call FITSImageDens(image,im,IMDIM,Rmax,ii)

	enddo

	return
	end


	subroutine FITSImageDens(image,im,IMDIM,Rmax,iadd)
	use Parameters
	IMPLICIT NONE
	integer IMDIM,iadd
	type(RPhiImage) image
	character*500 filename,cmd
	real*8 flux,w1,w2,il1,ir1,im(IMDIM,IMDIM)
	real*8 array(IMDIM,IMDIM),Rmax,cornerpix,pixscale
	integer i,k,i10,il10,il100,il1000,il10000,i1
	integer ir10,ir100,ir1000,ir10000

      integer status,unit,blocksize,bitpix,naxis,naxes(2)
      integer j,group,fpixel,nelements
      logical simple,extend,truefalse




	i10=(iadd)/10d0
	i1=(iadd)-10d0*i10

	ir10000=2d0*Rmax/10000d0
	ir1000=(2d0*Rmax-10000d0*ir10000)/1000d0
	ir100=(2d0*Rmax-1000d0*ir1000-10000d0*ir10000)/100d0
	ir10=(2d0*Rmax-100d0*ir100-1000d0*ir1000-10000d0*ir10000)/10d0
	ir1=2d0*Rmax-10d0*ir10-100d0*ir100-1000d0*ir1000-10000d0*ir10000

	if(ir1.ge.9.95d0) then
		ir10=ir10+1
		ir1=ir1-9.95d0
	endif
	if(ir10.ge.10) then
		ir100=ir100+1
		ir10=ir10-10
	endif
	if(ir100.ge.10) then
		ir1000=ir1000+1
		ir100=ir100-10
	endif
	if(ir1000.ge.10) then
		ir10000=ir10000+1
		ir1000=ir1000-10
	endif


	write(filename,'(a,"DensImage",i1,i1,"_fov",i1,i1,i1,i1,f3.1,".fits")')
     & outdir(1:len_trim(outdir)),i10,i1,ir10000,ir1000,ir100,ir10,ir1
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
      naxes(1)=IMDIM
      naxes(2)=IMDIM
      extend=.true.

C     write the required header keywords
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do i=1,IMDIM
	do j=1,IMDIM
		array(j,i)=im(i,j)
	enddo
	enddo

C     write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)
      call ftpprd(unit,group,fpixel,nelements,array,status)

C     write another optional keyword to the header
      call ftpkys(unit,'CTYPE1','X [AU]','X-coordinate',status)
      call ftpkys(unit,'CTYPE2','Y [AU]','Y-coordinate',status)

      call ftpkyd(unit,'CRPIX1',0.5d0,10,'Corner of the image',status)
      call ftpkyd(unit,'CRPIX2',0.5d0,10,'Corner of the image',status)

		cornerpix=-Rmax
      call ftpkyd(unit,'CRVAL1',cornerpix,10,'Center of the image',status)
      call ftpkyd(unit,'CRVAL2',cornerpix,10,'Center of the image',status)

		pixscale=2d0*Rmax/real(IMDIM)
      call ftpkyd(unit,'CDELT1',pixscale,10,'Pixel scale (in AU)',status)
      call ftpkyd(unit,'CDELT2',pixscale,10,'Pixel scale (in AU)',status)

	

C     close the file and free the unit number
      call ftclos(unit, status)
      call ftfiou(unit, status)


	return
	end





	subroutine MakeFluxImage(fluxcontr,Rmax,IMDIM,nintegrate,angle,lam0,basename)
	use Parameters
	IMPLICIT NONE
	type(RPhiImage) image
	integer IMDIM
	real*8 im(IMDIM,IMDIM),w1,w2,width,dx,dy,gasdev,angle
	real*8 lam0,Rmax,i1,tot,ran2,fluxcontr(0:D%nR,0:D%nTheta,2)
	integer ix,iy,i,j,k,i10,nintegrate
	integer il10000,il1000,il100,il10,jj,ii
	integer ir10000,ir1000,ir100,ir10,ir1,jp1,jp2
	real*8 x,y,R,phi,flux,il1,min,wp1,wp2,ranR,ranPhi
	real*8 fluxQ,fluxU,fluxV,immin
	character*500 filename
	character(len=*) :: basename
	logical hit(IMDIM,IMDIM)
	
	write(*,'("Creating flux image")')
	write(9,'("Creating flux image")')
	
	hit=.false.

	image%nr=D%nR
	image%nphi=(D%nTheta-1)*2d0+1
	if(.not.allocated(image%image)) then
		allocate(image%image(image%nr,image%nphi))
		allocate(image%Phi(image%nphi))
		allocate(image%R(image%nr))
	endif

	image%angle=angle
	image%lam=lam0
	
	do i=1,image%nr
		image%R(i)=D%R(i)
	enddo
	do i=1,D%nTheta-1
		image%Phi(i)=D%thet(i)
		image%Phi(image%nphi-i+1)=pi-D%thet(i)
	enddo
	image%Phi(D%nTheta)=pi/2d0
	image%nphi=image%nphi-1

	do i=1,image%nr
	do j=1,D%nTheta-1
		image%image(i,j)=fluxcontr(i,j,1)
		image%image(i,image%nphi+1-j)=fluxcontr(i,j,2)
	enddo
	enddo
	
	im=0d0

	do i=1,image%nr-1
	call tellertje(i,image%nr-1)
	do j=1,image%nphi
		w2=2d0*pi*AU**2*(image%R(i+1)-image%R(i))/real(image%nphi)
		do k=1,nintegrate
			ranR=ran2(idum)
			ranPhi=ran2(idum)

			R=image%R(i)+(image%R(i+1)-image%R(i))*ranR
			phi=pi*(real(j-1)+ranPhi)/real(image%nphi)
			if(j.ne.image%nphi) then
				phi=image%Phi(j)+ranPhi*(image%Phi(j+1)-image%Phi(j))
			else
				phi=image%Phi(j)+ranPhi*(pi-image%Phi(j))
			endif
			
			flux=image%image(i,j)

			x=-R*cos(phi)
			x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1
			ix=x
			y=R*sin(phi)	
			y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1
			iy=y
			if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
				im(ix,iy)=im(ix,iy)+flux/real(2*nintegrate)
				hit(ix,iy)=.true.
				y=-R*sin(phi)	
				y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1
				iy=y
				im(ix,iy)=im(ix,iy)+flux/real(2*nintegrate)
				hit(ix,iy)=.true.
			endif
		enddo
	enddo
	enddo

	call FITSImageFlux(image,im,IMDIM,Rmax,basename)

	return
	end


	subroutine FITSImageFlux(image,im,IMDIM,Rmax,basename)
	use Parameters
	IMPLICIT NONE
	integer IMDIM,iadd
	type(RPhiImage) image
	character*500 filename,cmd
	character(len=*) basename
	real*8 flux,w1,w2,il1,im(IMDIM,IMDIM),i1
	real*8 array(IMDIM,IMDIM),Rmax,cornerpix,pixscale
	integer i,k,i10,il10,il100,il1000,il10000

      integer status,unit,blocksize,bitpix,naxis,naxes(2)
      integer j,group,fpixel,nelements
      logical simple,extend,truefalse


	i10=(180d0*image%angle/pi)/10d0
	i1=(180d0*image%angle/pi)-10d0*i10
	if(i1.ge.9.95d0) then
		i1=i1-9.95d0
		i10=i10+1
	endif
	il10000=image%lam/10000d0
	il1000=(image%lam-10000d0*il10000)/1000d0
	il100=(image%lam-1000d0*il1000-10000d0*il10000)/100d0
	il10=(image%lam-100d0*il100-1000d0*il1000-10000d0*il10000)/10d0
	il1=image%lam-10d0*il10-100d0*il100-1000d0*il1000-10000d0*il10000

	write(filename,'(a,a,"_i",i1,f3.1,"_l",i1,i1,i1,i1,f4.2,".fits")')
     & outdir(1:len_trim(outdir)),trim(basename),i10,i1,il10000,il1000,il100,il10,il1
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
      naxes(1)=IMDIM
      naxes(2)=IMDIM
      extend=.true.

C     write the required header keywords
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do i=1,IMDIM
	do j=1,IMDIM
		array(j,i)=im(i,j)
	enddo
	enddo

C     write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)
      call ftpprd(unit,group,fpixel,nelements,array,status)

C     write another optional keyword to the header
      call ftpkys(unit,'CTYPE1','X [AU]','X-coordinate',status)
      call ftpkys(unit,'CTYPE2','Y [AU]','Y-coordinate',status)

      call ftpkyd(unit,'CRPIX1',0.5d0,10,'Corner of the image',status)
      call ftpkyd(unit,'CRPIX2',0.5d0,10,'Corner of the image',status)

		cornerpix=-Rmax
      call ftpkyd(unit,'CRVAL1',cornerpix,10,'Center of the image',status)
      call ftpkyd(unit,'CRVAL2',cornerpix,10,'Center of the image',status)

		pixscale=2d0*Rmax/real(IMDIM)
      call ftpkyd(unit,'CDELT1',pixscale,10,'Pixel scale (in AU)',status)
      call ftpkyd(unit,'CDELT2',pixscale,10,'Pixel scale (in AU)',status)

	

C     close the file and free the unit number
      call ftclos(unit, status)
      call ftfiou(unit, status)


	return
	end




	subroutine TELFWHM(image,tel,Rmax,FWHM1,FWHM2)
c	Rmax,IMDIM,nintegrate,width,SNR,Diam)
	use Parameters
	IMPLICIT NONE
	type(Telescope) tel
	type(RPhiImage) image
	integer IMDIM
c	parameter(IMDIM=500)
	real*8 w1,w2,dx,dy,gasdev,SNR
	real*8 wx,wy,wi,coswi,sinwi,cos2pa,sin2pa,Q,U
	real*8 Eu,Ex,Ey,EE

	real*8,allocatable :: im(:,:),imQ(:,:),imU(:,:),imV(:,:),imP(:,:)
	real*8,allocatable :: im1(:),im2(:),xx(:)
	real*8 lam0,Rmax,i1,tot,ran2,max,psfdev,nphot,err
	integer ix,iy,i,j,k,i10,nintegrate,inphot,scat_temp
	integer il10000,il1000,il100,il10,jj,number_invalid
	integer ir10000,ir1000,ir100,ir10,ir1,jp1,jp2,ir2
	real*8 x,y,R,phi,flux,il1,min,wp1,wp2,ranR,ranPhi
	real*8 fluxQ,fluxU,fluxV,immin,poidev
	real*8 sint,cost,sin2t,cos2t,z,wr1,wr2,FWHM1,FWHM2
	character*500 filename
	type(photon) phot

	real*8 al11,ar11,i11,c11,p11,s11
	real*8 al12,ar12,i12,c12,p12,s12
	real*8 al21,ar21,i21,c21,p21,s21
	real*8 al22,ar22,i22,c22,p22,s22

	IMDIM=tel%npixel
	nintegrate=tel%nint

	allocate(im(IMDIM,IMDIM))
	allocate(im1(IMDIM))
	allocate(im2(IMDIM))
	allocate(xx(IMDIM))
	lam0=image%lam

c	wi=widthi*pi/180d0
c	coswi=cos(wi)
c	sinwi=sin(wi)
	
	write(*,'("Creating rectangular image at wavelength: ",f10.3)') lam0
	write(9,'("Creating rectangular image at wavelength: ",f10.3)') lam0

	im=0d0

	do i=1,image%nr-1
	call tellertje(i,image%nr-1)
	do j=1,image%nphi
		w2=2d0*pi*AU**2*(image%R(i+1)-image%R(i))/real(image%nphi)
		do k=1,nintegrate
			ranR=ran2(idum)
			ranPhi=ran2(idum)

			R=image%R(i)+(image%R(i+1)-image%R(i))*ranR
			phi=pi*(real(j-1)+ranPhi)/real(image%nphi)

			if(ranPhi.lt.0.5d0) then
				wp1=0.5d0-ranPhi
				jp1=j-1
				if(jp1.eq.0) jp1=1
				wp2=ranPhi+0.5d0
				jp2=j
			else
				wp1=ranPhi-0.5d0
				jp1=j+1
				if(jp1.eq.image%nphi+1) jp1=image%nphi
				wp2=1.5d0-ranPhi
				jp2=j
			endif
			
			wr1=(1d0-ranR)
			ir1=i
			wr2=ranR
			ir2=i+1
			
			flux=(wr1*(wp1*log10(image%image(ir1,jp1))+wp2*log10(image%image(ir1,jp2)))
     &				+wr2*(wp1*log10(image%image(ir2,jp1))+wp2*log10(image%image(ir2,jp2))))
 			flux=(10d0**flux)*w2*R

c			wx=gasdev(idum)*widthy
c			wy=gasdev(idum)*widthx
			
			x=-R*cos(phi+D%PA*pi/180d0)!+wx*coswi+wy*sinwi
			x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1d0
			ix=x
			y=R*sin(phi+D%PA*pi/180d0)!+wy*coswi+wx*sinwi
			y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1d0
			iy=y
			if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
				im(ix,iy)=im(ix,iy)+flux/real(2*nintegrate)
			endif

c			wx=gasdev(idum)*widthy
c			wy=gasdev(idum)*widthx

			x=-R*cos(-phi+D%PA*pi/180d0)!+wx*coswi+wy*sinwi
			x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1d0
			ix=x
			y=R*sin(-phi+D%PA*pi/180d0)!+wy*coswi+wx*sinwi
			y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1d0
			iy=y
			if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
				im(ix,iy)=im(ix,iy)+flux/real(2*nintegrate)
			endif
		enddo
	enddo
	enddo

	call AddPlanet(image,im,imQ,imU,imV,Rmax,IMDIM)

	scat_temp=scat_how
	scat_how=0

	if(tel%psffile.ne.' ') then
		call ConvolutionFILE(im,imQ,imU,imV,IMDIM,tel%psffile)
	else if(tel%D.gt.0d0.or.tel%width.gt.0d0) then
		call Convolution(im,imQ,imU,imV,IMDIM,lam0,tel%D,tel%D2,tel%spider
     &			,tel%mask,tel%wmask,tel%owa,tel%strehl,Rmax*2d0*parsec/D%distance
     &			,tel%width,tel%snoise)
	endif
	
	scat_how=scat_temp

	do i=1,IMDIM
	do j=1,IMDIM
		im(i,j)=im(i,j)*1d23*image%zscale/D%distance**2
	enddo
	enddo

	do i=1,IMDIM
		im1(i)=sum(im(i,1:IMDIM))
		im2(i)=sum(im(1:IMDIM,i))
	enddo

	do i=1,IMDIM
		if(im1(i).lt.1d-30) im1(i)=1d-30
		if(im2(i).lt.1d-30) im2(i)=1d-30
	enddo	

	do i=1,IMDIM
		xx(i)=real(i)*Rmax*2d0/real(IMDIM)-Rmax
	enddo
	call FitGauss(xx,im1,IMDIM,FWHM1)
	FWHM1=dabs(2d0*sqrt(2d0*log(2d0))*FWHM1)
	call FitGauss(xx,im2,IMDIM,FWHM2)
	FWHM2=dabs(2d0*sqrt(2d0*log(2d0))*FWHM2)
	
	il10000=image%lam/10000d0
	il1000=(image%lam-10000d0*il10000)/1000d0
	il100=(image%lam-1000d0*il1000-10000d0*il10000)/100d0
	il10=(image%lam-100d0*il100-1000d0*il1000-10000d0*il10000)/10d0
	il1=image%lam-10d0*il10-100d0*il100-1000d0*il1000-10000d0*il10000
	write(filename,'(a,"FWHM_l",i1,i1,i1,i1,f4.2,".dat")')
     & outdir(1:len_trim(outdir)),il10000,il1000,il100,il10,il1

	open(unit=90,file=filename,RECL=1000)
	do i=1,IMDIM
		write(90,*) xx(i)*image%rscale,im2(i),im1(i)
	enddo
	close(unit=90)
	
	deallocate(im)
	deallocate(im1)
	deallocate(im2)
	deallocate(xx)

	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine FWHM(image,IMDIM,nintegrate,FWHM1,FWHM2,width)
	use Parameters
	IMPLICIT NONE
	type(RPhiImage) image
	integer IMDIM
c	parameter(IMDIM=500)
	real*8 w1,w2,dx,dy,gasdev,width,Diam
	real*8 lam0,Rmax,i1,tot,ran2,wx,wy,widthtot
	integer ix,iy,i,j,k,i10,nintegrate
	integer il10000,il1000,il100,il10,jj
	integer ir10000,ir1000,ir100,ir10,ir1,jp1,jp2,ir2
	real*8 x,y,R,phi,flux,il1,min,wp1,wp2,ranR,ranPhi,xx(IMDIM)
	real*8 im1(IMDIM),im2(IMDIM),FWHM1,FWHM2
	real*8 chi2,chi0,wr1,wr2
	character*500 filename
	
	lam0=image%lam

	write(*,'("Creating rectangular image at wavelength: ",f10.3)') lam0
	write(9,'("Creating rectangular image at wavelength: ",f10.3)') lam0

	im1=0d0
	im2=0d0
	
	Rmax=D%R(D%nR)

	widthtot=width!sqrt(width**2+Diam/image%lam)

	do i=1,image%nr-1
	call tellertje(i,image%nr-1)
	do j=1,image%nphi
		w2=2d0*pi*AU**2*(image%R(i+1)-image%R(i))/real(image%nphi)
		do k=1,nintegrate
			ranR=ran2(idum)
			ranPhi=ran2(idum)

			R=image%R(i)+(image%R(i+1)-image%R(i))*ranR
			phi=pi*(real(j-1)+ranPhi)/real(image%nphi)

			if(ranPhi.lt.0.5d0) then
				wp1=0.5d0-ranPhi
				jp1=j-1
				if(jp1.eq.0) jp1=1
				wp2=ranPhi+0.5d0
				jp2=j
			else
				wp1=ranPhi-0.5d0
				jp1=j+1
				if(jp1.eq.image%nphi+1) jp1=image%nphi
				wp2=1.5d0-ranPhi
				jp2=j
			endif
			
			wr1=(1d0-ranR)
			ir1=i
			wr2=ranR
			ir2=i+1
			
			flux=(wr1*(wp1*image%image(ir1,jp1)+wp2*image%image(ir1,jp2))
     &				+wr2*(wp1*image%image(ir2,jp1)+wp2*image%image(ir2,jp2)))*w2*R

			wx=gasdev(idum)*widthtot
			wy=gasdev(idum)*widthtot
			x=-R*cos(phi+D%PA*pi/180d0)+wx
			x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1
			ix=x
			y=R*sin(phi+D%PA*pi/180d0)+wy
			y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1
			iy=y
			if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
				im1(ix)=im1(ix)+flux/real(2*nintegrate)
				im2(iy)=im2(iy)+flux/real(2*nintegrate)
			endif
			wx=gasdev(idum)*widthtot
			wy=gasdev(idum)*widthtot
			x=-R*cos(-phi+D%PA*pi/180d0)+wx
			x=(real(IMDIM)*(x+Rmax)/(2d0*Rmax))+1
			ix=x
			y=R*sin(-phi+D%PA*pi/180d0)+wy
			y=(real(IMDIM)*(y+Rmax)/(2d0*Rmax))+1
			iy=y
			if(ix.le.IMDIM.and.iy.le.IMDIM.and.ix.gt.0.and.iy.gt.0) then
				im1(ix)=im1(ix)+flux/real(2*nintegrate)
				im2(iy)=im2(iy)+flux/real(2*nintegrate)
			endif
		enddo
	enddo
	enddo

	do i=1,IMDIM
		if(im1(i).lt.1d-30) im1(i)=1d-30
		if(im2(i).lt.1d-30) im2(i)=1d-30
	enddo	

	do i=1,IMDIM
		xx(i)=real(i)*D%R(D%nR)*2d0/real(IMDIM)-D%R(D%nR)
	enddo
	call FitGauss(xx,im1,IMDIM,FWHM1)
	FWHM1=2d0*sqrt(2d0*log(2d0))*FWHM1
	call FitGauss(xx,im2,IMDIM,FWHM2)
	FWHM2=2d0*sqrt(2d0*log(2d0))*FWHM2
	
	il10000=image%lam/10000d0
	il1000=(image%lam-10000d0*il10000)/1000d0
	il100=(image%lam-1000d0*il1000-10000d0*il10000)/100d0
	il10=(image%lam-100d0*il100-1000d0*il1000-10000d0*il10000)/10d0
	il1=image%lam-10d0*il10-100d0*il100-1000d0*il1000-10000d0*il10000
	write(filename,'(a,"FWHM_l",i1,i1,i1,i1,f4.2,".dat")')
     & outdir(1:len_trim(outdir)),il10000,il1000,il100,il10,il1

	open(unit=90,file=filename,RECL=1000)
	do i=1,IMDIM
		write(90,*) xx(i)*image%rscale,im2(i),im1(i)
	enddo
	close(unit=90)

	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine TraceDens(image)
	use Parameters
	IMPLICIT NONE
	real*8 phi
	integer i,j,k,ii,jj
	type(RPhiImage) image
	real*8 wT1,wT2,emis(0:D%nR,D%nTheta),frac,wl1,wl2
	integer ilam1,ilam2,iT,ip,jp,kp,jj1,jj2,djj,njj

	image%lam=0d0

	do i=1,image%nr
		call tellertje(i,image%nr)
		do j=1,image%nphi
			image%image(i,j)=0d0
			do k=1,image%p(i,j)%n
				ip=image%p(i,j)%i(k)
				jp=image%p(i,j)%j(k)
				kp=image%p(i,j)%k(k)
				image%image(i,j)=image%image(i,j)+image%p(i,j)%v(k)*C(ip,jp)%dens
			enddo
		enddo
	enddo

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

