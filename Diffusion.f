c 20071022 MM: Smoothing before diffusion is turned off due to instabilities
c 20080317 MM: Changed the diffusion to use the symetry of the system
c-----------------------------------------------------------------------
c This subroutine computes temperature diffusion in the optically thick
c approximation. It uses the Roseland approximation for energy diffusion
c NphotMin is the minimum number of photons that have passed through a
c cell. For cells with less photons, energy diffusion is used.
c-----------------------------------------------------------------------
	subroutine Diffuse(NPhotMin,dTMax)
	use Parameters
	IMPLICIT NONE
	logical assigned(0:D%nR,0:D%nTheta)
	integer number(0:D%nR,0:D%nTheta),N,i,j,jm,jp,j0,INFO,l,ia
	integer ii,iT,NPhotMin,NRHS,ncor,iter,niter,iopac,number_invalid
	integer,allocatable :: celi(:),celj(:),IWORK(:)
	real*8,allocatable :: Dc(:),T4(:),M(:,:),r(:),t(:),cooling(:)
	real*8 dx0,dx1,dx2,f1,f2,Rr,Rt,Rabs,dTMax,determineT,Etot,Dav,pow,gasdev,Kext,dens
	real*8 Krjm,rhojm,rjm,tjm,Krjp,rhojp,rjp,tjp,T4jm,T4jp,T4j0,prgopt(4),ErrorE,Error0
	logical,allocatable :: ok(:)
	type(photon) phot
	real*8 Rrad,Rthet,Ra,kp,ComputeCoolingFraction

	niter=3
	if(FLD) niter=10
	
	N=0
	ncor=0
	do i=0,D%nR
		do j=0,D%nTheta
			assigned(i,j)=.false.
			number(i,j)=0
		enddo
	enddo


	
c first number all cells that have to be included
c number(i,j) is the number in the array of cell i,j
c celi(i),celj(i) are the cell-coordinates of array element i
	C(1:D%nR-1,1:D%nTheta-1)%diff=.false.
	do i=2,D%nR-2
		do j=2,D%nTheta-1
			if(C(i,j)%Ni.lt.NphotMin.or.(C(i,j)%dT/C(i,j)%T).gt.dTMax) then
			if(C(i,j)%thick.or.C(i,j)%Ni.le.1) then
				ncor=ncor+1
				C(i,j)%diff=.true.
				if(.not.assigned(i,j)) then
					N=N+1
					number(i,j)=N
					assigned(i,j)=.true.
				endif
				if(.not.assigned(i-1,j)) then
					N=N+1
					number(i-1,j)=N
					assigned(i-1,j)=.true.
				endif
				if(.not.assigned(i+1,j)) then
					N=N+1
					number(i+1,j)=N
					assigned(i+1,j)=.true.
				endif
				if(.not.assigned(i,j-1)) then
					N=N+1
					number(i,j-1)=N
					assigned(i,j-1)=.true.
				endif
				if(.not.assigned(i,j+1).and.j.ne.(D%nTheta-1)) then
					N=N+1
					number(i,j+1)=N
					assigned(i,j+1)=.true.
				endif
			endif
			endif
		enddo
	enddo
	
	if(N.lt.1) return
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Diffusive approximation when low photon statistics")')
	write(*,'("Number of photons below: ",i8)') NPhotMin
	write(*,'("Error on T above:        ",f5.3)') dTMax
	write(*,'("Number of cells:         ",i8)') ncor
	write(9,'("Diffusive approximation when low photon statistics")')
	write(9,'("Number of photons below: ",i8)') NPhotMin
	write(9,'("Error on T above:        ",f5.3)') dTMax
	write(9,'("Number of cells:         ",i8)') ncor

	allocate(celi(N))
	allocate(celj(N))
	allocate(r(N))
	allocate(t(N))
	allocate(IWORK(10*N*N))
	allocate(T4(N))
	allocate(M(N,N))
	allocate(ok(0:N))
	allocate(Dc(N))
	allocate(cooling(N))


	do i=1,D%nR-1
		do j=1,D%nTheta-1
			if(assigned(i,j)) then
				celi(number(i,j))=i
				celj(number(i,j))=j
				r(number(i,j))=D%R_av(i)
				t(number(i,j))=D%theta_av(j)
			endif
		enddo
	enddo


c*****************************************************
c*****************************************************
c logaritmic gridding in r
c*****************************************************

	do iter=1,niter

c now compute the roseland mean opacities Kr
c also store the density and temperature information in the arrays
c label those cells that have sufficient photon statistics
c fill the T4 vector with zeros, except when suff. photon statistics
	do i=1,N
		iT=C(celi(i),celj(i))%T/dT
		if(iT.le.1) iT=1
		if(iT.ge.(TMAX-1)) iT=TMAX-1
		call DiffCoeffCell(celi(i),celj(i),iT)
		Kext=C(celi(i),celj(i))%KDext-C(celi(i),celj(i))%KDQHP
		dens=0d0
		do ii=1,ngrains
			if(.not.Grain(ii)%qhp) dens=dens+C(celi(i),celj(i))%w(ii)*C(celi(i),celj(i))%dens
		enddo
		Dc(i)=1d0/(Kext*dens)
		if(((D%R(D%nR)-D%R(celi(i)))*AU/Dc(i)).lt.2d0) C(celi(i),celj(i))%thick=.false.
		if(C(celi(i),celj(i))%Ni.lt.NPhotMin
     &		.or.(C(celi(i),celj(i))%dT/C(celi(i),celj(i))%T).gt.dTMax) then
			ok(i)=.false.
			if(celi(i).eq.1.or.celi(i).eq.D%nR-1.or.celj(i).eq.1) then
				ok(i)=.true.
			endif
			if(.not.C(celi(i),celj(i))%thick.and.C(celi(i),celj(i))%Ni.gt.1) then
				ok(i)=.true.
			endif
		else
			ok(i)=.true.
		endif
c		cooling(i)=0d0
c		if(iter.ne.1) then
c			call tellertje(i,N)
c			cooling(i)=ComputeCoolingFraction(celi,celj,ok,i,N)
c		endif
	enddo

	if(FLD) then
	do i=1,N
c the radial part
		if(celi(i).gt.1) then
			rjm=D%R_av(celi(i)-1)
			T4jm=C(celi(i)-1,celj(i))%T**4
		else
			rjm=D%R_av(celi(i))
			T4jm=C(celi(i),celj(i))%T**4
		endif
		if(celi(i).lt.(D%nR-1)) then
			rjp=D%R_av(celi(i)+1)
			T4jp=C(celi(i)+1,celj(i))%T**4
		else
			rjp=D%R_av(celi(i))
			T4jp=C(celi(i),celj(i))%T**4
		endif

		Rrad=(T4jp-T4jm)/(rjp-rjm)
c the theta part
		if(celj(i).ne.(D%nTheta-1)) then
			tjm=D%theta_av(celj(i)+1)
			T4jm=C(celi(i),celj(i)+1)%T**4
		else
			tjm=pi-D%theta_av(celj(i))
			T4jm=C(celi(i),celj(i))%T**4
		endif
		if(celj(i).gt.1) then
			tjp=D%theta_av(celj(i)-1)
			T4jp=C(celi(i),celj(i)-1)%T**4
		else
			tjp=D%theta_av(celj(i))
			T4jp=C(celi(i),celj(i))%T**4
		endif

		Rthet=(T4jp-T4jm)/((tjp-tjm)*r(i))
c absolute value
		Etot=C(celi(i),celj(i))%T**4
		if(C(celi(i),celj(i))%T.lt.dT) Etot=dT**4

		Ra=sqrt(Rrad**2+Rthet**2)*Dc(i)/Etot
		Dc(i)=Dc(i)*(2d0+Ra)/(6d0+3d0*Ra+Ra**2)
	enddo
	endif

	ok(0)=.false.

	do i=1,N
		if(ok(i)) then
			M(i,1:N)=0d0
			M(i,i)=1d0
			T4(i)=C(celi(i),celj(i))%T**4
		else
c the radial part
			M(i,1:N)=0d0
			j0=number(celi(i),celj(i))
			jm=number(celi(i)-1,celj(i))
			jp=number(celi(i)+1,celj(i))

			dx0=log(r(jp)/r(jm))
			dx1=log(r(jp)/r(j0))
			dx2=log(r(j0)/r(jm))

			f1=Dc(jp)*r(jp)
			f2=Dc(jm)*r(jm)

			M(i,jm)=M(i,jm)+f2/(dx2*dx0*r(j0)**3)
			M(i,j0)=M(i,j0)-f2/(dx2*dx0*r(j0)**3)

			M(i,jp)=M(i,jp)+f1/(dx1*dx0*r(j0)**3)
			M(i,j0)=M(i,j0)-f1/(dx1*dx0*r(j0)**3)

c the theta part
			j0=number(celi(i),celj(i))
			if(celj(i).ne.(D%nTheta-1)) then
				jm=number(celi(i),celj(i)+1)
				tjm=t(jm)
			else
				jm=number(celi(i),celj(i))
				tjm=pi-t(jm)
			endif
			jp=number(celi(i),celj(i)-1)

			dx0=t(jp)-tjm
			dx1=t(jp)-t(j0)
			dx2=t(j0)-tjm

			f1=Dc(jp)*sin(t(jp))
			f2=Dc(jm)*sin(tjm)

			M(i,jm)=M(i,jm)+f2/(dx2*dx0*r(j0)**2*sin(t(j0)))
			M(i,j0)=M(i,j0)-f2/(dx2*dx0*r(j0)**2*sin(t(j0)))

			M(i,jp)=M(i,jp)+f1/(dx1*dx0*r(j0)**2*sin(t(j0)))
			M(i,j0)=M(i,j0)-f1/(dx1*dx0*r(j0)**2*sin(t(j0)))

			T4(i)=0d0
		endif
	enddo

	NRHS=1
	call DGESV( N, NRHS, M, N, IWORK, T4, N, INFO )
	if(INFO.ne.0) then
		write(*,'("Matrix solution failed!",i5)') INFO
		write(9,'("Matrix solution failed!",i5)') INFO
	endif

	do i=1,N
		if(.not.ok(i)) then
			if(T4(i).gt.0d0) C(celi(i),celj(i))%T=T4(i)**0.25
			if(T4(i).lt.0d0) C(celi(i),celj(i))%T=dT
			C(celi(i),celj(i))%FradR=0d0
			C(celi(i),celj(i))%FradZ=0d0

			iT=(C(celi(i),celj(i))%T)/dT
			if(iT.le.1) iT=1
			if(iT.ge.(TMAX-1)) iT=TMAX-1 
			Etot=0d0
			do ii=1,ngrains
			   do iopac=1,Grain(ii)%nopac
			      Etot=Etot+C(celi(i),celj(i))%w(ii)*C(celi(i),celj(i))%wopac(ii,iopac)*Grain(ii)%Kp(iopac,iT)
			   enddo
			enddo
			C(celi(i),celj(i))%EJv=Etot

			if(C(celi(i),celj(i))%T.gt.(real(TMAX-1)*dT)) C(celi(i),celj(i))%T=real(TMAX-1)*dT
			if(C(celi(i),celj(i))%T.lt.1d0) then
				C(celi(i),celj(i))%T=1d0
			endif
			if(number_invalid(C(celi(i),celj(i))%T).ne.0) C(celi(i),celj(i))%T=1d0
			if(.not.tcontact) C(celi(i),celj(i))%TP(1:ngrains)=C(celi(i),celj(i))%T
			if(computeLRF) then
				iT=(C(celi(i),celj(i))%T/dT+0.5d0)
				if(iT.le.1) iT=1
				if(iT.ge.(TMAX-1)) iT=TMAX-1
				kp=0d0
				do ii=1,ngrains
					if(.not.Grain(ii)%qhp) then
						do iopac=1,Grain(ii)%nopac
							kp=kp+(Grain(ii)%Kp(iopac,iT)*C(celi(i),celj(i))%w(ii)*C(celi(i),celj(i))%wopac(ii,iopac))
						enddo
					endif
				enddo
				C(celi(i),celj(i))%LRF(1:nlam)=BB(1:nlam,iT)*C(celi(i),celj(i))%EJv/kp
				C(celi(i),celj(i))%nLRF(1:nlam)=100
			endif
		endif
	enddo
	
	enddo

	deallocate(celi)
	deallocate(celj)
	deallocate(r)
	deallocate(t)
	deallocate(IWORK)
	deallocate(T4)
	deallocate(M)
	deallocate(ok)
	deallocate(Dc)
	deallocate(cooling)

	return
	end



	real*8 function ComputeCoolingFraction(celi,celj,ok,i0,N)
	use Parameters
	IMPLICIT NONE
	integer N
	integer celi(N),celj(N),i,ilam,nexit,iexit
	integer inext,jnext,iT,i0,ii,iopac,j
	logical ok(N)
	type(photon) phot
	real*8 ran2,Rad,Theta,phi,tau,v,f,Kext,tot
	real*8,allocatable :: spec(:)
	
	nexit=4
	allocate(spec(nlam))
	spec=0d0
	do iexit=1,nexit
		do ilam=1,nlam
			phot%i=celi(i0)
			phot%j=celj(i0)
			Rad=ran2(idum)
			Rad=sqrt(D%R(phot%i)**2*Rad+D%R(phot%i+1)**2*(1d0-Rad))
			Theta=ran2(idum)
			Theta=D%Theta(phot%j)*Theta+D%Theta(phot%j+1)*(1d0-Theta)
			phot%z=Rad*Theta
			phi=1d-4
			phot%x=Rad*sqrt(1d0-Theta**2)*sin(phi)
			phot%y=Rad*sqrt(1d0-Theta**2)*cos(phi)
			phi=pi*ran2(idum)
			phot%vx=0d0
			phot%vy=cos(phi)
			phot%vz=sin(phi)
			phot%onEdge=.false.
			i=i0
			tau=0d0
1			do while(phot%i.gt.0.and.phot%i.lt.D%nR)
				call Trace2edge(phot,v,inext,jnext)
				Kext=0d0
				do ii=1,ngrains
					do iopac=1,Grain(ii)%nopac
						Kext=Kext+Grain(ii)%Kext(iopac,ilam)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
					enddo
				enddo
				tau=tau+v*Kext*AU*C(phot%i,phot%j)%dens
				phot%i=inext
				phot%j=jnext
				phot%x=phot%x+phot%vx*v
				phot%y=phot%y+phot%vy*v
				phot%z=phot%z+phot%vz*v
				phot%onEdge=.true.
			enddo
			spec(ilam)=spec(ilam)+exp(-tau)/real(nexit)
		enddo
	enddo

	iT=C(celi(i0),celj(i0))%T/dT
	if(iT.le.1) iT=1
	if(iT.ge.(TMAX-1)) iT=TMAX-1

	spec(1:nlam)=spec(1:nlam)*BB(1:nlam,iT)
	call integrate(spec,f)
	spec(1:nlam)=BB(1:nlam,iT)
	call integrate(spec,tot)
	f=f/tot

	ComputeCoolingFraction=f

	deallocate(spec)
	
	return
	end
	
	