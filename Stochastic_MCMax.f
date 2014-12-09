	module PAHModule
	IMPLICIT NONE
	real*8,allocatable :: BBQHP(:,:)
	end module PAHModule

	subroutine Stochastic(niter)
	use Parameters
	use InputOutput
	IMPLICIT NONE
	integer niter,ii,i,j,iopac
	type(Photon) phot
	real*8 determineTP
	logical outputionization
	character*500 filename
	
	if(niter.ne.0.and.NphotUV.gt.0.and.UV_PAH) call CreateLRF(NphotUV,10000)
	
	if(qhp_solver.eq.0) then
c use the PAH module from Michiel Min
		call StochasticMC(niter)
	else if(qhp_solver.eq.1) then
c use the PAH module from Kees Dullemond
		call PAHMCMax(niter)
	else if(qhp_solver.eq.2) then
c use the PAH module from Kees Dullemond
		call PAHequilibrium(niter)
	else
		write(*,'("QHP solver unknown")')
		write(9,'("QHP solver unknown")')
		stop
	endif

	do i=1,D%nR-1
		do j=1,D%nTheta-1
			do ii=1,ngrains
				if(Grain(ii)%qhp) then
					phot%i=i
					phot%j=j
					if(niter.eq.0) then
						phot%E=0d0
						do iopac=1,Grain(ii)%nopac
							phot%E=phot%E+Grain(ii)%Kpabsstar(iopac)*C(i,j)%wopac(ii,iopac)
						enddo
						phot%E=0.5d0*(1d0-sqrt(1d0-(D%Rstar/D%R_av(i))**2))*phot%E*D%Lstar/(pi*D%Rstar**2)
					else
						phot%E=C(i,j)%EJvQHP(Grain(ii)%qhpnr)/C(i,j)%w(ii)
					endif
					C(i,j)%Tqhp(Grain(ii)%qhpnr)=determineTP(phot,ii)
					if(C(i,j)%Tqhp(Grain(ii)%qhpnr).lt.2.8) C(i,j)%Tqhp(Grain(ii)%qhpnr)=2.8
				endif
			enddo
		enddo
	enddo
	do j=1,D%nTheta-1
		C(0,j)%Tqhp(1:nqhp)=C(1,j)%Tqhp(1:nqhp)
	enddo

	outputionization=.false.
	do ii=1,ngrains
		if(Grain(ii)%qhp.and.Grain(ii)%nopac.gt.1) outputionization=.true.
	enddo
	if(outputionization) then
		if(outputfits) then
			write(filename,'(a,"/ionization.fits.gz")') outdir(1:len_trim(outdir))
		else
			write(filename,'(a,"/ionization.dat")') outdir(1:len_trim(outdir))
		endif
		call outputstruct(filename,(/'G0     ','NE     '/),2,0)
	endif

	return
	end

	subroutine StochasticMC(niter)
	use Parameters
	use PAHModule
	IMPLICIT NONE
	integer i,j,ii,niter,itemp,ir,l,iopac,iT
	real*8 tot,tau,temp0,temp1,Planck
	real*8,allocatable :: Kabs(:)

	allocate(Kabs(nlam))
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Computing multi photon PAH emissivity")')
	write(9,'("Computing multi photon PAH emissivity")')
	write(*,'("Using Monte Carlo solver")')
	write(9,'("Using Monte Carlo solver")')

	temp0=10d0
	temp1=100000d0

	do itemp=1,NTQHP
		tgridqhp(itemp) = temp0 * (temp1/temp0)**
     %                      ((itemp-1.d0)/(NTQHP-1.d0))
	enddo

	if(.not.allocated(BBQHP)) then
		allocate(BBQHP(NTQHP,nlam))
		do itemp=1,NTQHP
			do i=1,nlam
				BBQHP(itemp,i)=Planck(tgridqhp(itemp),lam(i))
			enddo
		enddo
	endif

	do i=1,D%nR-1
		call tellertje(i,D%nR-1)
		do j=1,D%nTheta-1
			tau=0d0
			if(C(i,j)%Ni.le.0.or.niter.eq.0) then
				do ir=0,i
					do ii=1,ngrains
						tau=tau+C(ir,j)%dens*(D%R(ir+1)-D%R(ir))*AU*C(ir,j)%w(ii)*Grain(ii)%Kpstar(1)
					enddo
				enddo
				if(niter.eq.0) tau=0d0
				C(i,j)%LRF(1:nlam)=D%Fstar(1:nlam)/(4d0*pi*D%R_av(i)**2)*exp(-tau)
			endif
			do ii=1,ngrains
				if(Grain(ii)%qhp.and.Grain(ii)%nopac.gt.1) call PAHionization(ii,i,j)
			enddo
			if(niter.ne.0.or.((tau.lt.10d0.and.tau.gt.1d0).or.j.eq.1)) then
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
						do l=1,nlam
							Kabs(l)=0d0
							do iopac=1,Grain(ii)%nopac
								Kabs(l)=Kabs(l)+C(i,j)%wopac(ii,iopac)*Grain(ii)%Kabs(iopac,l)
							enddo
						enddo
						call Stochastic_MCMax(nu,nlam,Grain(ii)%Nc,Grain(ii)%Mc,Kabs(1:nlam)
     &		,Grain(ii)%Td_qhp,C(i,j)%LRF(1:nlam),C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)
     &		,tgridqhp,C(i,j)%tdistr(Grain(ii)%qhpnr,1:NTQHP),NTQHP,idum)
						call integrate(C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam),tot)
						if(tot.gt.1d-200) then
							C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)/tot
						else
							iT=C(i,j)%T
							if(iT.lt.3) iT=3
							if(iT.gt.TMAX) iT=TMAX
							C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=0d0
							do iopac=1,Grain(ii)%nopac
								C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)+
     &	C(i,j)%wopac(ii,iopac)*Grain(ii)%Kabs(iopac,1:nlam)*BB(1:nlam,iT)
							enddo
							call integrate(C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam),tot)
							C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)/tot
						endif
 	 				endif
				enddo
			else if(tau.le.1d0.and.j.ne.1) then
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
						C(i,j)%tdistr(Grain(ii)%qhpnr,1:NTQHP)=C(i,1)%tdistr(Grain(ii)%qhpnr,1:NTQHP)
						C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(i,1)%QHP(Grain(ii)%qhpnr,1:nlam)
					endif
				enddo
			else
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
						C(i,j)%tdistr(Grain(ii)%qhpnr,1:NTQHP)=0d0
						C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(1,1)%QHP(Grain(ii)%qhpnr,1:nlam)
					endif
				enddo
			endif
		enddo
	enddo


	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	deallocate(BBQHP)
	deallocate(Kabs)

	return
	end
	

	subroutine Stochastic_MCMax(nu,nlam,Na,Ma0,Kabs,Td,LRF,QHPemis,Tgrid,Tdistr,Nt,idum)
	use PAHModule
	IMPLICIT NONE
	integer nlam,i,j,idum,iT,ip,Ne,Nt,tel,Ntel,iTnew,ilam,iT0
	real*8 LRF(nlam),Kabs(nlam),nu(nlam),dnu(nlam),Td,Ma,Na,QHPemis(nlam)
	real*8 U(Nt),pi,c,h,Planck,T,Mass,Nphot(nlam),k,tot,ran2,Ma0
	real*8 Cool(Nt),Tdistr(Nt),r,dU,time,time0,wphot,T0,NphotSum(nlam)
	real*8 spec2(nlam),tot2,w1,w2,E,Tgrid(Nt),temp0,temp1,mp
	parameter(pi=3.1415926536)
	parameter(c=2.9979d10)		! cm/s
	parameter(h=6.6261d-27)	! cm^2 g/s
	parameter(k=1.3807d-16)	! cm^2 g / s^2 / K
	parameter(mp=1.67e-24)	! g

	dnu(1)=(nu(1)-nu(2))/2d0
	do i=2,nlam-1
		dnu(i)=(nu(i-1)-nu(i+1))/2d0
	enddo
	dnu(nlam)=(nu(nlam-1)-nu(nlam))/2d0

	Ma=Ma0*mp
	
	Mass=Na*Ma
	Nphot=LRF*Kabs*Mass*dnu/(h*nu)

	do i=1,Nt
		T=Tgrid(i)
		U(i)=3d0*Na*k*(T-Td*atan(T/Td))
		do j=1,nlam
			spec2(j)=BBQHP(i,j)*kabs(j)*Mass
		enddo
		call integrate(spec2,Cool(i))
	enddo

	spec2(1)=LRF(1)*Kabs(1)
	do i=2,nlam
		spec2(i)=spec2(i-1)+LRF(i)*Kabs(i)*nu(i)
	enddo

	NphotSum=Nphot
	do i=2,nlam
		NphotSum(i)=NphotSum(i-1)+NphotSum(i)
	enddo

	tot=sum(Nphot(1:nlam))

	Ne=(Na)**0.125
	if(Ne.lt.25) Ne=25
	if(Ne.gt.100) Ne=100
	Ne=Ne/2
	Ne=Ne*2

	Tdistr=0d0

	call integrate(LRF*Kabs,tot2)
	tot2=tot2*Mass
	call hunt(Cool,Nt,tot2,iT0)
	if(iT0.lt.2) then
		iT0=2
		T0=Tgrid(iT0)
	else
		w1=Cool(iT0+1)-tot2
		w2=tot2-Cool(iT0)
		T0=(w1*Tgrid(iT0)+w2*Tgrid(iT0+1))/(w1+w2)
	endif

	do j=1,nlam

	r=ran2(idum)*spec2(nlam)
	call hunt(spec2,nlam,r,ilam)
	ilam=ilam+1
	wphot=0d0

	iT=iT0
	T=T0
	time=0d0

	do i=1,Ne
		if(i.gt.2) wphot=Nphot(ilam)/(LRF(ilam)*Kabs(ilam)*nu(ilam))
		if(i.ne.Ne/2) then
			r=ran2(idum)*tot
			call hunt(NphotSum,nlam,r,ip)
			w1=NphotSum(ip+1)-r
			w2=r-NphotSum(ip)
			E=h*(w1*nu(ip)+w2*nu(ip+1))/(w1+w2)
		else
			ip=ilam
			E=h*nu(ip)
		endif
		
		if(T.lt.Tgrid(Nt)) then
			call hunt(Tgrid,Nt,T,iT)
			w1=Tgrid(iT+1)-T
			w2=T-Tgrid(iT)
		else
			iT=Nt-1
			w1=0d0
			w2=1d0
		endif
		tot2=E+(w1*U(iT)+w2*U(iT+1))/(w1+w2)
		call hunt(U,Nt,tot2,iT)
		if(iT.lt.1) then
			iT=1
			T=Tgrid(1)
		else if(iT.ge.(Nt-1)) then
			iT=Nt-1
			T=Tgrid(Nt-1)
		else
			w1=U(iT+1)-tot2
			w2=tot2-U(iT)
			T=(w1*Tgrid(iT)+w2*Tgrid(iT+1))/(w1+w2)
		endif

		time0=time-log(ran2(idum))/tot
		if(T.lt.Tgrid(Nt)) then
			call hunt(Tgrid,Nt,T,iT)
			w1=Tgrid(iT+1)-T
			w2=T-Tgrid(iT)
		else
			iT=Nt-1
			w1=0d0
			w2=1d0
		endif
		iTnew=iT
		if(iTnew.lt.1) iTnew=1
		dU=(w1*U(iT)+w2*U(iT+1))/(w1+w2)-U(iTnew)

		do while((time+dU/Cool(iT)).lt.time0.and.iT.gt.1)
			time=time+dU/Cool(iT)
			Tdistr(iT)=Tdistr(iT)+(dU/Cool(iT))*wphot
			T=Tgrid(iT)
			iTnew=iT-1
			if(iTnew.lt.1) iTnew=1
			dU=U(iT)-U(iTnew)
			iT=iTnew
		enddo
		dU=(time0-time)*Cool(iT)

		if(T.lt.Tgrid(Nt)) then
			call hunt(Tgrid,Nt,T,iT)
			w1=Tgrid(iT+1)-T
			w2=T-Tgrid(iT)
		else
			iT=Nt-1
			w1=0d0
			w2=1d0
		endif
		tot2=(w1*U(iT)+w2*U(iT+1))/(w1+w2)-dU
		call hunt(U,Nt,tot2,iT)
		if(iT.lt.1) then
			iT=1
			T=Tgrid(1)
		else if(iT.ge.(Nt-1)) then
			iT=Nt-1
			T=Tgrid(Nt-1)
		else
			w1=U(iT+1)-tot2
			w2=tot2-U(iT)
			T=(w1*Tgrid(iT)+w2*Tgrid(iT+1))/(w1+w2)
		endif
		time=time0

		Tdistr(iT)=Tdistr(iT)+dU/Cool(iT)*wphot!/2d0
	enddo
1	continue
	enddo
	
	tot=sum(Tdistr(1:Nt))
	Tdistr=Tdistr/tot
	QHPemis=0d0
	do i=2,Nt
		do j=1,nlam
			QHPemis(j)=QHPemis(j)+BBQHP(i,j)*kabs(j)*Tdistr(i)
		enddo
	enddo
	
	return
	end
	
	
	


	subroutine increaseQHPeq(phot,ii)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 E1,kp0,kp1,tot,T0,T1,Planck,spec(nlam)
	integer i,j,ii,iopac,l,iqhp,iT0,iT1,epsT0,epsT1

	iqhp=Grain(ii)%qhpnr

	E1=C(phot%i,phot%j)%EabsQHP(iqhp)/(C(phot%i,phot%j)%mass*C(phot%i,phot%j)%w(ii))
	T0=C(phot%i,phot%j)%Tqhp(iqhp)
	j=int(T0/dT)

	kp0=0d0
	do iopac=1,Grain(ii)%nopac
		kp0=kp0+Grain(ii)%Kp(iopac,j)*C(phot%i,phot%j)%wopac(ii,iopac)
	enddo
	do i=j,TMAX-1
		kp1=0d0
		do iopac=1,Grain(ii)%nopac
			kp1=kp1+Grain(ii)%Kp(iopac,i+1)*C(phot%i,phot%j)%wopac(ii,iopac)
		enddo
		if(kp0.le.E1.and.kp1.ge.E1) then
			C(phot%i,phot%j)%Tqhp(iqhp)=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dT
			goto 1
		endif
		kp0=kp1
	enddo
c not found, starting from 1 K
	j=1
	kp0=0d0
	do iopac=1,Grain(ii)%nopac
		kp0=kp0+Grain(ii)%Kp(iopac,j)*C(phot%i,phot%j)%wopac(ii,iopac)
	enddo
	do i=j,TMAX-1
		kp1=0d0
		do iopac=1,Grain(ii)%nopac
			kp1=kp1+Grain(ii)%Kp(iopac,i+1)*C(phot%i,phot%j)%wopac(ii,iopac)
		enddo
		if(kp0.le.E1.and.kp1.ge.E1) then
			C(phot%i,phot%j)%Tqhp(iqhp)=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dT
			goto 1
		endif
		kp0=kp1
	enddo
	C(phot%i,phot%j)%Tqhp(iqhp)=real(TMAX-1)*dT

1	continue

	T1=C(phot%i,phot%j)%Tqhp(iqhp)

	iT0=int(T0/dT)
	iT1=int(T1/dT)

	epsT0=T0-real(iT0)*dT
	epsT1=T1-real(iT1)*dT
	if(iT0.ge.TMAX-1) iT0=TMAX-2
	if(iT1.ge.TMAX-1) iT1=TMAX-1
	if(iT0.lt.1) iT0=1
	if(iT1.lt.2) iT1=2

	do l=1,nlam
		spec(l)=0d0
		do iopac=1,Grain(ii)%nopac
			spec(l)=spec(l)+C(phot%i,phot%j)%wopac(ii,iopac)*Grain(ii)%Kabs(iopac,l)
		enddo
	enddo

	if(iT0.eq.iT1) then
		do l=1,nlam
			spec(l)=spec(l)*(BB(l,iT0+1)-BB(l,iT0))
		enddo
	else
		do l=1,nlam
			spec(l)=spec(l)*(epsT1*BB(l,iT1+1)+(1d0-epsT1)*BB(l,iT1)-epsT0*BB(l,iT0+1)-(1d0-epsT0)*BB(l,iT0))
		enddo
	endif

	call integrate(spec,tot)
	C(phot%i,phot%j)%QHP(Grain(ii)%qhpnr,1:nlam)=spec(1:nlam)/tot

	return
	end


	
	

	subroutine PAHequilibrium(niter)
	use Parameters
	IMPLICIT NONE
	integer i,j,ii,niter,iqhp,iopac,l
	real*8 determineTP,tot,Planck
	type(Photon) phot

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Computing PAH emissivity assuming equilibrium")')
	write(9,'("Computing PAH emissivity assuming equilibrium")')

	do i=1,D%nR-1
		call tellertje(i,D%nR-1)
		do j=1,D%nTheta-1
			do ii=1,ngrains
				if(Grain(ii)%qhp) then
					if(niter.eq.0) C(i,j)%LRF(1:nlam)=D%Fstar(1:nlam)/(4d0*pi*D%R_av(i)**2)
					if(Grain(ii)%nopac.gt.1) call PAHionization(ii,i,j)
					iqhp=Grain(ii)%qhpnr
					phot%i=i
					phot%j=j
					if(niter.eq.0) then
						phot%E=0d0
						do iopac=1,Grain(ii)%nopac
							phot%E=phot%E+Grain(ii)%Kpabsstar(iopac)*C(i,j)%wopac(ii,iopac)
						enddo
						phot%E=0.5d0*(1d0-sqrt(1d0-(D%Rstar/D%R_av(i))**2))*phot%E*D%Lstar/(pi*D%Rstar**2)
					else
						phot%E=C(i,j)%EJvQHP(iqhp)/C(i,j)%w(ii)
					endif
					C(i,j)%Tqhp(iqhp)=determineTP(phot,ii)
					if(C(i,j)%Tqhp(iqhp).lt.2.8) C(i,j)%Tqhp(iqhp)=2.8
					do l=1,nlam
						C(i,j)%QHP(iqhp,l)=0d0
						do iopac=1,Grain(ii)%nopac
							C(i,j)%QHP(iqhp,l)=C(i,j)%QHP(iqhp,l)+Planck(C(i,j)%Tqhp(iqhp),lam(l))*
     &								C(i,j)%wopac(ii,iopac)*Grain(ii)%Kabs(iopac,l)
						enddo
					enddo
					call integrate(C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam),tot)
					C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)/tot
				endif
			enddo
		enddo
	enddo


	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	return
	end
	

