c New version of TraceMono that actually does the line radiative integration
	subroutine TraceLine(image,iTr,flux,Nphot,NphotStar,linefile,popfile,velo,dvelo,nvelo,abun,poptype,i_up,i_low)
	use Parameters
	IMPLICIT NONE
	integer i,j,k,ii,jj,Nphot,NphotStar,nvelo,iv,iTr
	real*8 lam0,tau(nvelo),phi,flux(nvelo),velo(nvelo),dvelo,abun
	type(RPhiImage) image
	real*8 tau_e,tau_s,tau_a,w1,w2,nu0,exptau_e(nvelo),velo0
	real*8 wT1,wT2,emis(0:D%nR,D%nTheta),frac,wl1,wl2,Ksca
	integer ilam1,ilam2,iT,ip,jp,kp,jj1,jj2,djj,njj,irg,iopac
	real*8 scat(2,0:D%nR,D%nTheta),fact(nvelo)
	real*8 scatQ(2,0:D%nR,D%nTheta),scatU(2,0:D%nR,D%nTheta)
	real*8 scatV(2,0:D%nR,D%nTheta),nf,sf,fracirg(D%nTheta,360)
	character*500 fluxfile,linefile,popfile
	character*10 poptype
	integer il10,il100,il1000,il10000,i10,kk,nkk,i_low,i_up
	real*8 i1,il1,Planck,velo_emis(nvelo),velo_scat(nvelo),velo_tau(nvelo)
	real*8 velo_scatQ(nvelo),velo_scatU(nvelo),velo_scatV(nvelo),add
	real*8 im(nvelo),imQ(nvelo),imU(nvelo),imV(nvelo),gu,gl,x,h,T
	real*8 Aul,Bul,Blu,mol_mass,n_up(0:D%nR-1,D%nTheta-1),n_low(0:D%nR-1,D%nTheta-1)
	real*8 y1(nvelo),y2(nvelo),y3(nvelo),y4(nvelo),Resolution,clight,nm(0:D%nR-1,D%nTheta-1)
	real*8,allocatable :: velo_imageI(:,:,:),velo_imageQ(:,:,:),velo_imageU(:,:,:),velo_imageV(:,:,:)
	parameter(h=6.626068e-27) ! cm^2 g/s
	parameter(clight=2.9979d5) !km/s
	logical aniso
	
	aniso=.false.
	if(scat_how.eq.2) aniso=.true.

	tau_max=1d2

	poptype='LTE'
	call ComputeLTE(linefile,iTr,i_low,i_up,Aul,Bul,Blu,nu0,n_low,n_up,mol_mass)
	nm=mol_mass*abun
	lam0=clight*1d9/nu0

	if(popfile.ne.' ') call PopulationsProdimo(popfile,i_low,i_up,n_low,n_up,nm,mol_mass,poptype)

	do i=1,D%nR-1
	do j=1,D%nTheta-1
		n_up(i,j)=n_up(i,j)*nm(i,j)
		n_low(i,j)=n_low(i,j)*nm(i,j)
	enddo
	enddo

	if(lam0.lt.lam(1).or.lam0.gt.lam(nlam)) then
		write(*,'("Wavelength not in the grid: ",f10.3)') lam0
		return
	endif

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

	if(storescatt) then
		if(.not.scattcomputed(ilam1)) then
			do i=0,D%nR
			do j=0,D%nTheta
				C(i,j)%scattfield(1,0,1:2)=0d0
			enddo
			enddo
			call TraceMono(lam(ilam1),Nphot,image%angle,NphotStar)
			do i=0,D%nR
			do j=0,D%nTheta
				C(i,j)%scattfield(1,ilam1,1)=C(i,j)%scattfield(1,0,1)
				C(i,j)%scattfield(1,ilam1,2)=C(i,j)%scattfield(1,0,2)
			enddo
			enddo
			scattcomputed(ilam1)=.true.
			nscattcomputed(ilam1)=Nphot+NphotStar
		else if(nscattcomputed(ilam1).lt.Nphot) then
			do i=0,D%nR
			do j=0,D%nTheta
				C(i,j)%scattfield(1,0,1:2)=0d0
			enddo
			enddo
			call TraceMono(lam(ilam1),Nphot-nscattcomputed(ilam1),image%angle,NphotStar)
			do i=0,D%nR
			do j=0,D%nTheta
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
			do j=0,D%nTheta
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
			do j=0,D%nTheta
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
	else if(scattering.and.(Nphot.ne.0.or.NphotStar.ne.0)) then
		do i=0,D%nR
		do j=0,D%nTheta
			C(i,j)%scattfield=0d0
			if(aniso) then
				C(i,j)%scattQ=0d0
				C(i,j)%scattU=0d0
				C(i,j)%scattV=0d0
			endif
		enddo
		enddo
		call TraceMono(lam0,Nphot,image%angle,NphotStar)
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
				if(iT.gt.TMAX) iT=TMAX
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
			if(iT.gt.TMAX) iT=TMAX
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
	enddo
	enddo
	
	call setup_velo(linefile,velo,nvelo,dvelo,lam0,Aul,Blu,n_up,n_low,mol_mass)

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

	allocate(velo_imageI(nvelo,image%nr,image%nphi))
	if(aniso) then
		allocate(velo_imageQ(nvelo,image%nr,image%nphi))
		allocate(velo_imageU(nvelo,image%nr,image%nphi))
		allocate(velo_imageV(nvelo,image%nr,image%nphi))
	endif
	
	do i=1,image%nr
	call tellertje(i,image%nr)
	do j=1,image%nphi
		im=0d0
		if(aniso) then
			imQ=0d0
			imU=0d0
			imV=0d0
		endif
		tau=0d0
		fact=1d0
		
		do k=1,image%p(i,j)%n

		ip=image%p(i,j)%i(k)
		if(ip.ne.0) then
		
		jp=image%p(i,j)%j(k)
		kp=image%p(i,j)%k(k)
		irg=image%p(i,j)%irg(k)
		
		if(.not.storescatt.and.scattering.and.(Nphot+NphotStar).ne.0) then
		jj1=image%p(i,j)%jphi1(k)
		jj2=image%p(i,j)%jphi2(k)
		scat(kp,ip,jp)=0d0
		if(aniso) then
			scatQ(kp,ip,jp)=0d0
			scatU(kp,ip,jp)=0d0
			scatV(kp,ip,jp)=0d0
		endif
		djj=-1
		if(jj1.le.jj2) djj=1
		njj=0
		jj=jj1
1		scat(kp,ip,jp)=scat(kp,ip,jp)+C(ip,jp)%Albedo*C(ip,jp)%scattfield(irg,jj,kp)/(C(ip,jp)%V*fracirg(jp,irg))
		if(aniso) then
			scatQ(kp,ip,jp)=scatQ(kp,ip,jp)+C(ip,jp)%Albedo*C(ip,jp)%scattQ(irg,jj,kp)/(C(ip,jp)%V*fracirg(jp,irg))
			scatU(kp,ip,jp)=scatU(kp,ip,jp)+C(ip,jp)%Albedo*C(ip,jp)%scattU(irg,jj,kp)/(C(ip,jp)%V*fracirg(jp,irg))
			scatV(kp,ip,jp)=scatV(kp,ip,jp)+C(ip,jp)%Albedo*C(ip,jp)%scattV(irg,jj,kp)/(C(ip,jp)%V*fracirg(jp,irg))
		endif
		njj=njj+1
		if(jj.eq.jj2) goto 2
		jj=jj+djj
		goto 1
2		scat(kp,ip,jp)=scat(kp,ip,jp)/real(njj)
		if(aniso) then
			scatQ(kp,ip,jp)=scatQ(kp,ip,jp)/real(njj)
			scatU(kp,ip,jp)=scatU(kp,ip,jp)/real(njj)
			scatV(kp,ip,jp)=scatV(kp,ip,jp)/real(njj)
		endif
		endif

		tau_e=image%p(i,j)%v(k)*C(ip,jp)%dens*C(ip,jp)%Kext*AU

		nkk=2d0*abs(image%p(i,j)%velo1(k)-image%p(i,j)%velo2(k))/dvelo+1

		do kk=1,nkk
		
		velo0=image%p(i,j)%velo1(k)+(image%p(i,j)%velo2(k)-image%p(i,j)%velo1(k))*(real(kk)-0.5)/real(nkk)

		call velo_abs_ext(velo0,ip,jp,velo_emis,velo_tau,velo,nvelo,dvelo)

		velo_emis=velo_emis*1d-2+emis(ip,jp)*C(ip,jp)%Kext*C(ip,jp)%dens/(C(ip,jp)%gasdens*gas2dust)
		velo_tau=tau_e/real(nkk)+velo_tau*C(ip,jp)%gasdens*gas2dust*image%p(i,j)%v(k)*AU/real(nkk)

		velo_emis=velo_emis*C(ip,jp)%gasdens*gas2dust*(image%p(i,j)%v(k)*AU/real(nkk))/velo_tau
		
		velo_scat=scat(kp,ip,jp)*tau_e/(velo_tau*real(nkk))
		if(aniso) then
			velo_scatQ=scatQ(kp,ip,jp)*tau_e/(real(nkk)*velo_tau)
			velo_scatU=scatU(kp,ip,jp)*tau_e/(real(nkk)*velo_tau)
			velo_scatV=scatV(kp,ip,jp)*tau_e/(real(nkk)*velo_tau)
		endif

		if(sum(velo_tau)/real(nvelo).lt.1d-5) then
			im=im+(2d0*velo_scat+velo_emis)*velo_tau*fact
			if(aniso) then
				imQ=imQ+2d0*velo_scatQ*velo_tau*fact
				imU=imU+2d0*velo_scatU*velo_tau*fact
				imV=imV+2d0*velo_scatV*velo_tau*fact
			endif
			fact=fact*(1d0-velo_tau)
		else
			exptau_e=exp(-velo_tau)
			im=im+(2d0*velo_scat+velo_emis)*(1d0-exptau_e)*fact
			if(aniso) then
				imQ=imQ+2d0*velo_scatQ*(1d0-exptau_e)*fact
				imU=imU+2d0*velo_scatU*(1d0-exptau_e)*fact
				imV=imV+2d0*velo_scatV*(1d0-exptau_e)*fact
			endif
			fact=fact*exptau_e
		endif
		tau=tau+velo_tau
		
		enddo
		
		if(tau(1).gt.tau_max) goto 10
		endif
		enddo
		if(image%p(i,j)%hitstar.and.tracestar) then
			im=im+(D%Fstar(ilam1)*wl1+D%Fstar(ilam2)*wl2)*fact*dimstar/(pi*(D%R(0)*AU)**2)
		endif
10	continue
		velo_imageI(1:nvelo,i,j)=im
		if(aniso) then
			velo_imageQ(1:nvelo,i,j)=imQ
			velo_imageU(1:nvelo,i,j)=imU
			velo_imageV(1:nvelo,i,j)=imV
		endif
	enddo
	enddo

	flux=0d0
	do i=1,image%nr-1
	do k=1,image%nPhi
		w1=2d0*pi*abs(image%R(i))*AU**2/real(image%nPhi)
		w2=2d0*pi*abs(image%R(i+1))*AU**2/real(image%nPhi)
		if(k.ne.image%nPhi) then
			do iv=1,nvelo
				y1(iv)=velo_imageI(iv,i,k)
				y2(iv)=velo_imageI(iv,i,k+1)
			enddo
			call addcurves(y1,y2,velo,y3,nvelo)
			do iv=1,nvelo
				y1(iv)=velo_imageI(iv,i+1,k)
				y2(iv)=velo_imageI(iv,i+1,k+1)
			enddo
			call addcurves(y1,y2,velo,y4,nvelo)
			call addcurves(w1*y3,w2*y4,velo,im,nvelo)
			do iv=1,nvelo
				flux(iv)=flux(iv)+(image%R(i+1)-image%R(i))*im(iv)/4d0
			enddo
		else
			do iv=1,nvelo
				y1(iv)=velo_imageI(iv,i,k)
				y2(iv)=velo_imageI(iv,i,1)
			enddo
			call addcurves(y1,y2,velo,y3,nvelo)
			do iv=1,nvelo
				y1(iv)=velo_imageI(iv,i+1,k)
				y2(iv)=velo_imageI(iv,i+1,1)
			enddo
			call addcurves(y1,y2,velo,y4,nvelo)
			call addcurves(w1*y3,w2*y4,velo,im,nvelo)
			do iv=1,nvelo
				flux(iv)=flux(iv)+(image%R(i+1)-image%R(i))*im(iv)/4d0
			enddo
		endif
	enddo
	enddo

c	Resolution=600d0
c	Resolution=clight/Resolution
c	do iv=1,nvelo
c		y1(iv)=exp(-(velo(iv)/Resolution)**2)
c	enddo
c	y2=y1/sum(y1)
c	im=flux
c	flux=0d0
c	do iv=1,nvelo
c		do i=1,nvelo
c			ii=i-iv+nvelo/2
c			if(ii.lt.1) ii=1
c			if(ii.gt.nvelo) ii=nvelo
c			flux(i)=flux(i)+im(iv)*y2(ii)
c		enddo
c	enddo
		

c symemtric line profiles
c	im=flux
c	do iv=1,nvelo
c		flux(iv)=(im(iv)+im(nvelo+1-iv))/2d0
c	enddo

	image%flux=flux(nvelo/2)

	write(*,'("Total flux:     ",e15.3," Jy")') flux(1)*1e23/D%distance**2
	write(9,'("Total flux:     ",e15.3," Jy")') flux(1)*1e23/D%distance**2

	deallocate(velo_imageI)
	if(aniso) then
		deallocate(velo_imageQ)
		deallocate(velo_imageU)
		deallocate(velo_imageV)
	endif

	do i=1,D%nR-1
	do j=1,D%nTheta-1
		deallocate(C(i,j)%line_emis)
		deallocate(C(i,j)%line_abs)
	enddo
	enddo

	return
	end
	
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine TracePathLine(image,angle,nphi,nj0,s_theta,nstar,lam0)
	use Parameters
	IMPLICIT NONE
	character*500 input,tmp,specfile
	integer i,j,k,l,t1,t2,nstar,s_theta
	real*8 tau,angle,spec(nlam),distance,tottime,ct,lam0
	real*8 ph,inp,determineT,lam1,lam2,v,dvelo
	integer starttime,stoptime,starttrace,cr,l1,l2,nj0,ip,jp
	type(Photon) phot,photcount
	integer ilam,MinPhot,nabs,nj,nj2,nphi
	real*8 x,y,z,phi,theta,Albedo,w1,w2,fact(nlam),Rad,r,ran2
	logical hitstar
	real*8 extstar(nlam),R01(D%nTheta),t_offset
	type(RPhiImage) image

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
		image%nr=nj0*(D%nR-1)+(D%nTheta-1)*2*(D%nRfix+2)/s_theta+1+nstar
	else
		image%nr=(D%nR-2)/(-nj0)+1+(D%nTheta-1)*2*(D%nRfix+2)/s_theta+1+nstar
	endif

	nj=0
	do j=1,D%nTheta-1,s_theta
		nj=nj+1
		nj=nj+1
	enddo
	do j=1,D%nTheta-1,s_theta
		nj=nj+1
		nj=nj+1
	enddo
	do i=1,D%nRfix
	do j=1,D%nTheta-1,s_theta
		nj=nj+1
		nj=nj+1
	enddo
	enddo
	if(nj0.gt.0) then
		do i=1,D%nR-1
		do j=1,nj0
			nj=nj+1
		enddo
		enddo
	else
		do i=1,D%nR-1,-nj0
			nj=nj+1
		enddo
	endif
	do i=1,nstar-1
		nj=nj+1
	enddo
	nj=nj+1
	image%nr=nj+1

	
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
	do j=1,D%nTheta-1,s_theta
		nj=nj+1
		image%R(nj)=R01(D%nTheta-1)*abs(cos(pi/2d0-D%theta_av(j)+angle))
		nj=nj+1
		image%R(nj)=abs(R01(D%nTheta-1)*cos(D%theta_av(j)-pi/2d0+angle))
	enddo
	do j=1,D%nTheta-1,s_theta
		nj=nj+1
		image%R(nj)=R01(D%nTheta-1)*abs(cos(pi/2d0-D%theta_av(j)))
		nj=nj+1
		image%R(nj)=abs(R01(D%nTheta-1)*cos(D%theta_av(j)-pi/2d0))
	enddo
	do i=1,D%nRfix
	do j=1,D%nTheta-1,s_theta
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

	image%nr=nj+1

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
	
	do i=1,image%nr
	call tellertje(i,image%nr)
	t_offset=ran2(idum)-0.5d0
	do k=1,image%nPhi
		image%phi(k)=2d0*pi*(real(k)-0.5)/real(image%nPhi)+2d0*pi*t_offset/real(image%nPhi)
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
		allocate(image%p(i,k)%velo1(image%p(i,k)%n))
		allocate(image%p(i,k)%velo2(image%p(i,k)%n))

		call trace2dpathLine(phot,image%p(i,k))
	enddo
	enddo


	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine trace2dpathLine(phot,p)
	use Parameters
	IMPLICIT NONE
	type(Photon) phot
	real*8 tau(nlam),tau_a(nlam),tau_e(nlam),v,fact(nlam)
	real*8 Kext(nlam),Kabs(nlam),emis(nlam),spec(nlam)
	real*8 absfrac(nlam),wT1,wT2,x,y,z,phi,G,r
	integer inext,jnext,ntrace,iT,i,j,l1,l2,irgnext
	logical hitstar
	type(path) p
	parameter(G=6.67300d-8) ! in cm^3/g/s^2

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

c v=29.7859 * (1/a^(1/2)) * [y,-x,0]/sqrt(x^2+y^2)       in km/s

	r=sqrt(x*x+y*y+z*z)
	p%velo1(p%n)=0D0
c Keppler velocity
	p%velo1(p%n)=p%velo1(p%n)+1d-5*sqrt(G*D%Mstar*sin(D%theta_av(phot%j))/(AU*r))*(phot%vx*y-phot%vy*x)/sqrt(x**2+y**2)
c 20 km/s outflow
c	p%velo1(p%n)=p%velo1(p%n)-20d0*(phot%vx*x+phot%vy*y+phot%vz*z)/r


	x=phot%x+phot%vx*v
	y=phot%y+phot%vy*v
	z=phot%z+phot%vz*v
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

	r=sqrt(x*x+y*y+z*z)
	p%velo2(p%n)=0d0
c Keppler velocity
	p%velo2(p%n)=p%velo2(p%n)+1d-5*sqrt(G*D%Mstar*sin(D%theta_av(phot%j))/(AU*r))*(phot%vx*y-phot%vy*x)/sqrt(x**2+y**2)
c 20 km/s outflow
c	p%velo2(p%n)=p%velo2(p%n)-20d0*(phot%vx*x+phot%vy*y+phot%vz*z)/r

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


	subroutine setup_velo(linefile,velo,nvelo,dvelo,lam0,Aul,Blu,n_up,n_low,mol_mass)
c Aul, Bul, Blu	: Einstein coefficients
c n_up, n_low	: number of molecules per unit volume that are in the upper/lower level
c mol_mass		: molecule mass in proton masses
	use Parameters
	IMPLICIT NONE
	character*500 linefile
	integer nvelo,i,j,iv
	real*8 velo(nvelo),dvelo,lam0,velo_T,mp,tot,clight,mol_mass,nu0,gu,gl,x
	real*8 phi(nvelo),fact,Aul,Bul,Blu,n_up(0:D%nR-1,D%nTheta-1),n_low(0:D%nR-1,D%nTheta-1),T,h
	parameter(mp=1.67262158d-24) ! the proton mass in gram
	parameter(clight=2.9979d10) ! cm/s
	parameter(h=6.626068e-27) ! cm^2 g/s

	do iv=-(nvelo-1)/2,(nvelo-1)/2
		velo(iv+(nvelo-1)/2+1)=real(iv)*dvelo
	enddo

	nu0=clight/(lam0*1d-4)
	Bul=Aul*2d0*h*nu0**3/(clight**2)
	
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		allocate(C(i,j)%line_emis(nvelo))
		allocate(C(i,j)%line_abs(nvelo))

		if(useTgas) then
			T=C(i,j)%Tgas
		else
			T=C(i,j)%T
		endif

		velo_T=1d-5*sqrt(2d0*kb*T/(mp*mol_mass)) ! factor 1e-5 is to convert from cm/s to km/s
c		velo_T=velo_T+1d0
c		if(velo_T.lt.dvelo/2d0) velo_T=dvelo/2d0
		do iv=1,nvelo
			phi(iv)=(clight/(sqrt(pi)*velo_T*nu0*1d5))*exp(-(velo(iv)/velo_T)**2)
		enddo
		tot=sum(phi(1:nvelo))*dvelo/(lam0*1d-9)
		phi=phi/tot

		fact=h*nu0/(4d0*pi*mol_mass*mp)
		C(i,j)%line_emis=fact*n_up(i,j)*Aul*phi
		C(i,j)%line_abs=fact*(n_low(i,j)*Blu-n_up(i,j)*Bul)*phi
	enddo
	enddo
	
	return
	end
	

	subroutine velo_abs_ext(v,i,j,velo_emis,velo_tau,velo,nvelo,dvelo)
	use Parameters
	IMPLICIT NONE
	integer i,j,nvelo,iv,shift1,shift2,ii
	real*8 v,velo_emis(nvelo),velo_tau(nvelo),velo(nvelo),ran2,dvelo,w1,w2

	shift1=(v/dvelo)
	shift2=(v/dvelo+sign(1d0,v))
	w1=abs(real(shift2)*dvelo-v)/dvelo
	w2=1d0-w1

c	call hunt(velo,nvelo,v,shift)
c	if(ran2(idum).lt.abs((v-velo(shift))/(velo(2)-velo(1)))) shift=shift+1
c	shift=(nvelo+1)/2-shift
	do iv=1,nvelo
		ii=iv-shift1
		if(ii.lt.1) ii=1
		if(ii.gt.nvelo) ii=nvelo
		velo_tau(iv)=C(i,j)%line_abs(ii)*w1
		velo_emis(iv)=C(i,j)%line_emis(ii)*w1
		ii=iv-shift2
		if(ii.lt.1) ii=1
		if(ii.gt.nvelo) ii=nvelo
		velo_tau(iv)=velo_tau(iv)+C(i,j)%line_abs(ii)*w2
		velo_emis(iv)=velo_emis(iv)+C(i,j)%line_emis(ii)*w2
	enddo

	return
	end
	


	subroutine addcurves(y1in,y2in,x,y,n)
	IMPLICIT NONE
	integer n,i,imax1,imax2,j,nj,ii,s1,s2
	real*8 x(n),y1(n),y2(n),y(n)
	real*8 w1,w2,tot,cont_flux1,cont_flux2,y1in(n),y2in(n),tot1,tot2

	cont_flux1=(y1in(1)+y1in(n))/2d0
	y1=y1in-cont_flux1
	cont_flux2=(y2in(1)+y2in(n))/2d0	
	y2=y2in-cont_flux2

	imax1=(n-1)/2
	imax2=(n-1)/2
	do i=1,n
		if(abs(y1(i)).gt.y1(imax1)) imax1=i
		if(abs(y2(i)).gt.y2(imax2)) imax2=i
	enddo

	nj=abs(imax1-imax2)-1

	tot1=sum(y1)
	tot2=sum(y2)
	if((tot1/tot2).lt.1e-6.or.(tot2/tot1).lt.1e-6) nj=1

	if(nj.lt.1) then
		y=y1in+y2in
	else
		y=y1+y2
		do j=1,nj
			s1=sign(j,imax1-imax2)
			s2=s1-(imax1-imax2)
			w1=real(abs(s2))/real(abs(s1)+abs(s2))
			w2=real(abs(s1))/real(abs(s1)+abs(s2))

			do i=1,n
				ii=i-s1
				if(ii.lt.1) ii=1
				if(ii.gt.n) ii=n
				y(ii)=y(ii)+y1(i)*w1
				ii=i-s2
				if(ii.lt.1) ii=1
				if(ii.gt.n) ii=n
				y(ii)=y(ii)+y2(i)*w2
			enddo
		enddo
		y=2d0*y/real(nj+2)+cont_flux1+cont_flux2

		tot=sum(y1in(1:n)+y2in(1:n))

		y=y*tot/sum(y(1:n))
	endif

	return
	end
	

