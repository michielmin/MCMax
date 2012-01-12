	module PAHModule
	IMPLICIT NONE
	real*8,allocatable :: BBQHP(:,:)
	end module PAHModule

	subroutine Stochastic(niter)
	use Parameters
	use PAHModule
	IMPLICIT NONE
	integer i,j,ii,niter,itemp,ir
	real*8 tot,tau,temp0,temp1,Planck
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Computing multi photon PAH emissivity")')
	write(9,'("Computing multi photon PAH emissivity")')

	open(unit=90,file='LRF.dat',RECL=6000)
	do i=1,nlam
		write(90,*) lam(i),(C(j,1)%LRF(i),j=1,D%nR-1)
	enddo
	close(unit=90)

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
			if(niter.ne.0.or.((tau.lt.10d0.and.tau.gt.1d0).or.j.eq.1)) then
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
						call Stochastic_MCMax(nu,nlam,Grain(ii)%Nc,Grain(ii)%Mc,Grain(ii)%Kabs(1,1:nlam)
     &		,Grain(ii)%Td_qhp,C(i,j)%LRF(1:nlam),C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)
     &		,tgridqhp,C(i,j)%tdistr(Grain(ii)%qhpnr,1:NTQHP),NTQHP,idum)
						call integrate(C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam),tot)
						if(tot.gt.1d-200) then
							C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)/tot
						else
							C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=0d0
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

	if(i.eq.1.and.j.eq.1) then
		open(unit=42,file='QHP.dat',RECL=6000)
		do ii=1,nlam
			write(42,*) lam(ii),C(i,j)%QHP(1,ii)
		enddo
		close(unit=42)
	endif
	
		enddo
		do j=1,D%nTheta-1
			do ii=1,nqhp
				C(i,j)%Tqhp(ii)=0d0
				do itemp=1,NTQHP
					C(i,j)%Tqhp(ii)=C(i,j)%Tqhp(ii)+C(i,j)%tdistr(ii,itemp)*tgridqhp(itemp)**4
				enddo
				C(i,j)%Tqhp(ii)=C(i,j)%Tqhp(ii)**0.25d0
			enddo
		enddo
	enddo
	do j=1,D%nTheta-1
		C(0,j)%Tqhp(1:nqhp)=C(1,j)%Tqhp(1:nqhp)
	enddo


	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')


	deallocate(BBQHP)

	return
	end
	

	subroutine Stochastic_MCMax(nu,nlam,Na,Ma0,Kabs,Td,LRF,QHPemis,Tgrid,Tdistr,Nt,idum)
	use PAHModule
	IMPLICIT NONE
	integer nlam,i,j,idum,iT,ip,Ne,Nt,tel,Ntel,iTnew,ilam,iT0
	real*8 LRF(nlam),Kabs(nlam),nu(nlam),dnu(nlam),Td,Ma,Na,QHPemis(nlam)
	real*8 U(Nt),pi,c,h,Planck,T,Mass,Nphot(nlam),k,tot,ran2,Ma0
	real*8 Cool(Nt),Tdistr(Nt),r,dU,time,time0,wphot,T0
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

	tot=sum(Nphot(1:nlam))

	do i=1,Nt
		T=Tgrid(i)
		U(i)=3d0*Na*k*(T-Td*atan(T/Td))
		do j=1,nlam
			spec2(j)=BBQHP(i,j)*kabs(j)*Mass
		enddo
		call integrate(spec2,Cool(i))
	enddo

	do i=2,nlam
		Nphot(i)=Nphot(i-1)+Nphot(i)
	enddo
	spec2(1)=LRF(1)*Kabs(1)
	do i=2,nlam
		spec2(i)=spec2(i-1)+LRF(i)*Kabs(i)*nu(i)
	enddo

	Ne=(Na)**0.125
	if(Ne.lt.20) Ne=20
	Ne=20

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
			call hunt(Nphot,nlam,r,ip)
			w1=Nphot(ip+1)-r
			w2=r-Nphot(ip)
			E=h*(w1*nu(ip)+w2*nu(ip+1))/(w1+w2)
		else
			ip=ilam
			E=h*nu(ip)
		endif
		
		call hunt(Tgrid,Nt,T,iT)
		w1=Tgrid(iT+1)-T
		w2=T-Tgrid(iT)
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
		call hunt(Tgrid,Nt,T,iT)
		w1=Tgrid(iT+1)-T
		w2=T-Tgrid(iT)
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

		call hunt(Tgrid,Nt,T,iT)
		w1=Tgrid(iT+1)-T
		w2=T-Tgrid(iT)
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

		Tdistr(iT)=Tdistr(iT)+dU/Cool(iT)*wphot/2d0
	enddo
1	continue
	enddo
	
	tot=sum(Tdistr(2:Nt))
	Tdistr=Tdistr/tot
	QHPemis=0d0
	do i=2,Nt
		do j=1,nlam
			QHPemis(j)=QHPemis(j)+BBQHP(i,j)*kabs(j)*Tdistr(i)
		enddo
	enddo
	
	return
	end
	
