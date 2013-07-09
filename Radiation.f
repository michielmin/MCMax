	subroutine integrate(spec,L)
	use Parameters
	IMPLICIT NONE
	real*8 spec(nlam),L
	integer i
	
	L=0d0
	do i=1,nlam-1
		L=L+(spec(i)+(spec(i+1)-spec(i))/2d0)*dnu(i)
	enddo

	return
	end



	subroutine makePdis(spec,Pspec)
	use Parameters
	IMPLICIT NONE
	real*8 spec(nlam),Pspec(nlam),scale
	integer i
	
	Pspec(1)=(spec(1)+(spec(2)-spec(1))/2d0)*dnu(1)
	do i=2,nlam-1
		Pspec(i)=Pspec(i-1)+(spec(i)+(spec(i+1)-spec(i))/2d0)*dnu(i)
	enddo
	scale=Pspec(nlam-1)
	do i=1,nlam-1
		Pspec(i)=Pspec(i)/scale
	enddo
	return
	end


	
	subroutine interact(phot)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 spec(nlam),Albedo,ran2,T0,T1,increaseT,increaseTP
	real*8 x,y,z,phi,theta,r,Ksca,Kext,kp,Kabs,Edust,Egas
	real*8 epsT1,epsT0,KabsQHP,KabsTot,Ks(nlam),w1,w2,KscaR
	integer i,iT,iT0,iT1,l,ii,iopac
	type(Mueller) M

	if(C(phot%i,phot%j)%useFE.and.computeTgas.and.g2d_heat) then
		Egas=phot%E*C(phot%i,phot%j)%FE(0)
		Edust=phot%E*(1d0-C(phot%i,phot%j)%FE(0))
	else
		Egas=0d0
		Edust=phot%E
	endif

	Ksca=0d0
	Kext=0d0
	KabsQHP=0d0
	Ksca=(C(phot%i,phot%j)%KscaTot(phot%ilam1)*phot%wl1+C(phot%i,phot%j)%KscaTot(phot%ilam2)*phot%wl2)
	Kext=Ksca+(C(phot%i,phot%j)%KabsTot(phot%ilam1)*phot%wl1+C(phot%i,phot%j)%KabsTot(phot%ilam2)*phot%wl2)
	do i=1,ngrains
	do iopac=1,Grain(i)%nopac
		if(Grain(i)%qhp) KabsQHP=KabsQHP+(Grain(i)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(i)%Kabs(iopac,phot%ilam2)*phot%wl2)*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
	enddo
	enddo
	Kext=Kext+KabsQHP
	KabsTot=Kext-Ksca
	
	Albedo=Ksca/Kext

	if(.not.scattering.or.ran2(idum).gt.Albedo) then
c absorption
	phot%scatt=.false.
	phot%pol=.false.
	phot%Q=0d0
	phot%U=0d0
	phot%V=0d0

c normal absorption and reemission
	if(.not.tcontact) then
c No thermal contact
	spec(1:nlam)=0d0
	kp=0d0
	C(phot%i,phot%j)%Eabs=C(phot%i,phot%j)%Eabs+Edust
	do i=1,ngrains
	if(.not.Grain(i)%qhp) then
	do iopac=1,Grain(i)%nopac
		C(phot%i,phot%j)%EabsP(i)=C(phot%i,phot%j)%EabsP(i)+Edust*((Grain(i)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(i)%Kabs(iopac,phot%ilam2)*phot%wl2)*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac))/(Kext-Ksca)
	enddo
	T0=C(phot%i,phot%j)%TP(i)
	T1=increaseTP(phot,i)

	iT0=int(T0/dT)
	iT1=int(T1/dT)
	epsT0=T0-real(iT0)*dT
	epsT1=T1-real(iT1)*dT
	if(iT0.ge.TMAX-1) iT0=TMAX-2
	if(iT1.ge.TMAX-1) iT1=TMAX-1
	if(iT0.eq.iT1) then
		do iopac=1,Grain(i)%nopac
			do l=1,nlam
				spec(l)=spec(l)+(BB(l,iT0+1)-BB(l,iT0))*Grain(i)%Kabs(iopac,l)
     &					*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
			enddo
			kp=kp+(Grain(i)%Kp(iopac,iT0+1)-Grain(i)%Kp(iopac,iT0))*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
		enddo
	else
		do iopac=1,Grain(i)%nopac
			do l=1,nlam
				spec(l)=spec(l)+(epsT1*BB(l,iT1+1)+(1d0-epsT1)*BB(l,iT1)-epsT0*BB(l,iT0+1)-(1d0-epsT0)*BB(l,iT0))*
     &		Grain(i)%Kabs(iopac,l)*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
			enddo
			kp=kp+(epsT1*Grain(i)%Kp(iopac,iT1+1)+(1d0-epsT1)*Grain(i)%Kp(iopac,iT1)
     &			-epsT0*Grain(i)%Kp(iopac,iT0+1)-(1d0-epsT0)*Grain(i)%Kp(iopac,iT0))
     &			*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
     	enddo
	endif

	C(phot%i,phot%j)%TP(i)=T1
	endif
	enddo
	if(kp.ne.0d0.and.use_qhp) then
		spec(1:nlam)=spec(1:nlam)*(KabsTot-KabsQHP)/(kp*KabsTot)
		kp=1d0
	endif
	do i=1,ngrains
		do iopac=1,Grain(i)%nopac
		if(Grain(i)%qhp) then
			KabsQHP=(Grain(i)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(i)%Kabs(iopac,phot%ilam2)*phot%wl2)*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
			if(use_qhp) then
				spec(1:nlam)=spec(1:nlam)
     &	+KabsQHP*C(phot%i,phot%j)%QHP(Grain(i)%qhpnr,1:nlam)/KabsTot
			endif
		endif
		enddo
	enddo

	if(C(phot%i,phot%j)%useFE.and.computeTgas.and.g2d_heat) then
		spec=spec/kp
		kp=1d0
		T0=C(phot%i,phot%j)%Tgas
		iT0=int(T0/dT)
		spec(1:nlam)=spec(1:nlam)*Edust+Egas*BB(1:nlam,iT0)/BBint(iT0)
		kp=phot%E
	endif

	call emit(phot,spec,kp)
	C(phot%i,phot%j)%T=increaseT(phot)
	
	else
c Thermal contact (single temperature)
	C(phot%i,phot%j)%Eabs=C(phot%i,phot%j)%Eabs+Edust*(KabsTot-KabsQHP)/KabsTot
	T0=C(phot%i,phot%j)%T
	T1=increaseT(phot)


	iT0=int(T0/dT)
	iT1=int(T1/dT)

	epsT0=T0-real(iT0)*dT
	epsT1=T1-real(iT1)*dT
	if(iT0.ge.TMAX-1) iT0=TMAX-2
	if(iT1.ge.TMAX-1) iT1=TMAX-1

	if(iT0.eq.iT1) then
		do l=1,nlam
			spec(l)=(BB(l,iT0+1)-BB(l,iT0))*C(phot%i,phot%j)%KabsTot(l)
		enddo
		kp=0d0
		do i=1,ngrains
			if(.not.Grain(i)%qhp) then
				do iopac=1,Grain(i)%nopac
					kp=kp+(Grain(i)%Kp(iopac,iT0+1)-Grain(i)%Kp(iopac,iT0))*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
				enddo
			endif
		enddo
	else
		do l=1,nlam
			spec(l)=(epsT1*BB(l,iT1+1)+(1d0-epsT1)*BB(l,iT1)-epsT0*BB(l,iT0+1)-(1d0-epsT0)*BB(l,iT0))*C(phot%i,phot%j)%KabsTot(l)
		enddo
		kp=0d0
		do i=1,ngrains
			if(.not.Grain(i)%qhp) then
				do iopac=1,Grain(i)%nopac
					kp=kp+(epsT1*Grain(i)%Kp(iopac,iT1+1)+(1d0-epsT1)*Grain(i)%Kp(iopac,iT1)
     &			-epsT0*Grain(i)%Kp(iopac,iT0+1)-(1d0-epsT0)*Grain(i)%Kp(iopac,iT0))
     &			*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
				enddo
			endif
		enddo
	endif

	if(kp.ne.0d0.and.use_qhp) then
		spec(1:nlam)=spec(1:nlam)*(KabsTot-KabsQHP)/(kp*KabsTot)
		kp=1d0
	endif
	do i=1,ngrains
		if(Grain(i)%qhp) then
		do iopac=1,Grain(i)%nopac
			KabsQHP=(Grain(i)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(i)%Kabs(iopac,phot%ilam2)*phot%wl2)*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
			if(use_qhp) then
				spec(1:nlam)=spec(1:nlam)
     &	+KabsQHP*C(phot%i,phot%j)%QHP(Grain(i)%qhpnr,1:nlam)/KabsTot
			endif
		enddo
		endif
	enddo

	if(C(phot%i,phot%j)%useFE.and.computeTgas.and.g2d_heat) then
		spec=spec/kp
		kp=1d0
		T0=C(phot%i,phot%j)%Tgas
		iT0=int(T0/dT)
		spec(1:nlam)=spec(1:nlam)*Edust+Egas*BB(1:nlam,iT0)/BBint(iT0)
		kp=phot%E
	endif

	call emit(phot,spec,kp)
	
	C(phot%i,phot%j)%T=T1

	endif

	else
c scattering
	phot%scatt=.true.
	if(scat_how.eq.1) then
		call randomdirection(phot%vx,phot%vy,phot%vz)
		if(multiwav) then
			specemit(1:nlam)=specemit(1:nlam)*C(phot%i,phot%j)%KscaTot(1:nlam)
		endif
	else
		KscaR=Ksca*ran2(idum)
		do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			w1=C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
			w2=w1*phot%wl2*Grain(ii)%Ksca(iopac,phot%ilam2)
			KscaR=KscaR-w2
			if(KscaR.lt.0d0) then
				M=Grain(ii)%F(iopac,phot%ilam2)
				goto 1
			endif
			w1=w1*phot%wl1*Grain(ii)%Ksca(iopac,phot%ilam1)
			KscaR=KscaR-w1
			if(KscaR.lt.0d0) then
				M=Grain(ii)%F(iopac,phot%ilam1)
				goto 1
			endif
		enddo
		enddo
c Randomly picking a grain and wavelength went wrong for some reason
c use the slow backup determination of average scattering matrix instead
		M%IF11=0d0
		M%IF12=0d0
		M%F11=0d0
		M%F12=0d0
		M%F22=0d0
		M%F33=0d0
		M%F34=0d0
		M%F44=0d0
		do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			w1=C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)
			w2=w1*phot%wl2*Grain(ii)%Ksca(iopac,phot%ilam2)
			w1=w1*phot%wl1*Grain(ii)%Ksca(iopac,phot%ilam1)
			M%IF11=M%IF11+w1*Grain(ii)%F(iopac,phot%ilam1)%IF11+w2*Grain(ii)%F(iopac,phot%ilam2)%IF11
			M%IF12=M%IF12+w1*Grain(ii)%F(iopac,phot%ilam1)%IF12+w2*Grain(ii)%F(iopac,phot%ilam2)%IF12
			M%F11=M%F11+w1*Grain(ii)%F(iopac,phot%ilam1)%F11+w2*Grain(ii)%F(iopac,phot%ilam2)%F11
			M%F12=M%F12+w1*Grain(ii)%F(iopac,phot%ilam1)%F12+w2*Grain(ii)%F(iopac,phot%ilam2)%F12
			M%F22=M%F22+w1*Grain(ii)%F(iopac,phot%ilam1)%F22+w2*Grain(ii)%F(iopac,phot%ilam2)%F22
			M%F33=M%F33+w1*Grain(ii)%F(iopac,phot%ilam1)%F33+w2*Grain(ii)%F(iopac,phot%ilam2)%F33
			M%F34=M%F34+w1*Grain(ii)%F(iopac,phot%ilam1)%F34+w2*Grain(ii)%F(iopac,phot%ilam2)%F34
			M%F44=M%F44+w1*Grain(ii)%F(iopac,phot%ilam1)%F44+w2*Grain(ii)%F(iopac,phot%ilam2)%F44
		enddo
		enddo
		
1		call scatangle(phot,M)
		if(multiwav) then
			specemit(1:nlam)=0d0
			specemit(phot%ilam1)=phot%wl1
			specemit(phot%ilam2)=phot%wl2
			column(1:ngrains,1:ngrains2)=0d0
			if(exportprodimo) Kext_column=0d0
		endif
	endif
	
	endif
	
	return
	end


c------------------------------------------------------------------------
c This subroutine computes a scattering event for photon phot.
c
c The photon should be fully initialized. Used properties are:
c phot%i 			: radial cell number
c phot%j 			: theta cell number
c phot%ilam1		: wavelength number.
c phot%ilam2		: wavelength number.
c phot%wl1			: weight of lam(phot%ilam1)
c phot%wl2			: weight of lam(phot%ilam2)
c			The wavelength of the photon is:
c				phot%lam=phot%wl1*lam(phot%ilam1)+phot%wl2*lam(phot%ilam2)
c phot%pol			: logical determining if the photon is polarized
c phot%vx,vx,vz		: direction of propagation (changed on output)
c phot%Sx,Sy,Sz		: direction of the reference Q vector (changed on output)
c phot%E,Q,U,V		: Stokes vector of the photon (changed on output)
c
c The grain properties used are:
c M%F(i)	: contains all info on the Mueller matrix integrated over grains at wavelength point i
c	F(i)%F11,F12,F22,F33,F34,F44	: elements of the Mueller matrix (normalized)
c	F(i)%IF11						: integral of the F11 and F12 elements weighted with sin(theta)
c 
c------------------------------------------------------------------------
	subroutine scatangle(phot,M)
	use Parameters
	IMPLICIT NONE
	real*8 theta,Fr,ran2,Fi,FiOld,w1,w2,s1,s2,x,y,z
	real*8 E,Q,U,V,F,IF11,IF12,I0,I1,Itot,phi,wp1,P,P2,r,thet
	type(photon) phot
	integer i,ii,iopac
	type(Mueller) M

	IF11=M%IF11*180d0
	IF12=M%IF12*180d0

	if(phot%pol) then
c Integrate the total scattered intensity over all phi angles
c sin2phi and cos2phi are initialized already to be sin(2*phi) and cos(2*phi)
	Itot=0d0
	do i=1,180
		Itot=Itot+IF11*phot%E+IF12*(phot%Q*cos2phi(i)+phot%U*sin2phi(i))
	enddo
	Itot=ran2(idum)*Itot
c Determine the phi angle for scattering
	I1=0d0
	I0=0d0
	do i=1,180
		I0=I0+IF11*phot%E+IF12*(phot%Q*cos2phi(i)+phot%U*sin2phi(i))
		if(I0.gt.Itot) then
			wp1=(I0-Itot)/(I0-I1)
			phi=pi*(real(i)-wp1)/180d0
			goto 2
		endif
		I1=I0
	enddo
	i=180
	wp1=0d0
	phi=pi
2	continue
	else
c If incoming photon is unpolarized, phi angle is random.
	phi=ran2(idum)*pi
	endif

c Scattering plane is defined. Rotate the Stokes vector to the new plane.
	x=phot%Sx
	y=phot%Sy
	z=phot%Sz
c First rotate the axis
	call rotate(x,y,z,phot%vx,phot%vy,phot%vz,phi)
c Then rotate the Stokes vector to the new reference axis
	if(phot%pol) call RotateStokes(phot,x,y,z)
	phot%Sx=x
	phot%Sy=y
	phot%Sz=z

	Fr=180d0*(phot%E*M%IF11+phot%Q*M%IF12)

c Now determine the scattering angle theta
	Fr=ran2(idum)*Fr
	Fi=0d0
	FiOld=0d0
	do i=1,180
		thet=pi*(real(i)-0.5d0)/180d0
		Fi=Fi+pi*sin(thet)*(M%F11(i)*phot%E+M%F12(i)*phot%Q)
		if(Fi.gt.Fr) then
			theta=real(i)-(Fi-Fr)/(Fi-FiOld)
			theta=theta*pi/180d0
			goto 1
		endif
		FiOld=Fi
	enddo
	theta=pi
	i=180
1	continue
c Due to symmetry, could also have been -theta
	if(ran2(idum).lt.0.5d0) theta=-theta
	
c Now determine the scattered stokes vector.
	E=0d0
	Q=0d0
	U=0d0
	V=0d0
	E=E+M%F11(i)*phot%E
	E=E+M%F12(i)*phot%Q
	Q=Q+M%F12(i)*phot%E
	Q=Q+M%F22(i)*phot%Q
	U=U+M%F33(i)*phot%U
c if rotation is the counterclockwise F=-F
	F=sign(M%F34(i),theta)
	U=U+F*phot%V
	V=V-F*phot%U
	V=V+M%F44(i)*phot%V

c Now renormalize the photon package to conserve energy.
	if(E.ne.0d0) then
		phot%Q=phot%E*Q/E
		phot%U=phot%E*U/E
		phot%V=phot%E*V/E
	endif

c Rotate the propagation vector.
	call rotate(phot%vx,phot%vy,phot%vz,phot%Sx,phot%Sy,phot%Sz,theta)
	r=sqrt(phot%vx**2+phot%vy**2+phot%vz**2)
	phot%vx=phot%vx/r
	phot%vy=phot%vy/r
	phot%vz=phot%vz/r
	r=sqrt(phot%Sx**2+phot%Sy**2+phot%Sz**2)
	phot%Sx=phot%Sx/r
	phot%Sy=phot%Sy/r
	phot%Sz=phot%Sz/r

c Photon is now polarized !
	if(i.ne.1.and.i.ne.180) phot%pol=.true.

	return
	end

c------------------------------------------------------------------------
c This subroutine rotates the Stokes vector from 
c one reference plane to another.
c
c The initial reference plane is stored in the phot%Sx,Sy,Sz
c The new reference plane is x,y,z
c The Stokes vector phot%Q,U is changed, the phot%Sx,Sy,Sz are unaltered.
c Note that the vector x,y,z should have unit length!
c------------------------------------------------------------------------
	subroutine RotateStokes(phot,x,y,z)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 x,y,z,cost,sint,Q,U,P,P2,cos2t,sin2t,sPP2
	
	P=(phot%U/phot%E)**2+(phot%Q/phot%E)**2

	if(P.eq.0d0) return

	cost=phot%Sx*x+phot%Sy*y+phot%Sz*z
	sint=(x-
     &	(phot%Sx*(phot%vy**2+phot%vz**2)-phot%vx*(phot%vy*phot%Sy+phot%vz*phot%Sz))*cost)/
     &	(phot%vy*phot%Sz-phot%vz*phot%Sy)
	
	cos2t=cost**2-sint**2
	sin2t=2d0*sint*cost
	
	Q=phot%Q*cos2t+phot%U*sin2t
	U=-phot%Q*sin2t+phot%U*cos2t

	P2=(U/phot%E)**2+(Q/phot%E)**2

	sPP2=sqrt(P/P2)
	phot%U=sPP2*U
	phot%Q=sPP2*Q
	
	return
	end


	real*8 function increaseT(phot)
	use Parameters
	IMPLICIT NONE
	type(photon) phot,phot2
	real*8 E1,kp0,kp1,determineT,computeE
	integer i,j,ii,iopac

	E1=C(phot%i,phot%j)%Eabs/C(phot%i,phot%j)%mass
	j=int(C(phot%i,phot%j)%T/dT)
	kp0=computeE(phot%i,phot%j,j)

	do i=j,TMAX-1
		kp1=computeE(phot%i,phot%j,i+1)
		if(kp0.le.E1.and.kp1.ge.E1) then
			increaseT=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dT
			return
		endif
		kp0=kp1
	enddo
c not found, starting from 1 K
	phot2=phot
	phot2%E=E1
	increaseT=determineT(phot2)

	return
	end


	real*8 function increaseTP(phot,ii)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 E1,kp0,kp1
	integer i,j,ii,iopac

	E1=C(phot%i,phot%j)%EabsP(ii)/(C(phot%i,phot%j)%mass*C(phot%i,phot%j)%w(ii))
	j=int(C(phot%i,phot%j)%TP(ii)/dT)
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
			increaseTP=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dT
			return
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
			increaseTP=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dT
			return
		endif
		kp0=kp1
	enddo
	increaseTP=real(TMAX-1)*dT
	return
	end


	real*8 function computeE(i,j,iT)
	use Parameters
	IMPLICIT NONE
	integer iT,ii,iopac,i,j

	computeE=0d0
	do ii=1,ngrains
		if(.not.Grain(ii)%qhp) then
			do iopac=1,Grain(ii)%nopac
				computeE=computeE+(Grain(ii)%Kp(iopac,iT)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac))
			enddo
		endif
	enddo
	return
	end

	real*8 function determineT(phot)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 E1,E,T,Emin,Emax,computeE
	integer i,ii,iopac,iTmin,iTmax,iT0,iT

	E1=phot%E
	iTmin=1
	iTmax=TMAX

	iT=C(phot%i,phot%j)%T/dT
	if(iT.lt.1) iT=1
	if(iT.gt.TMAX-1) iT=TMAX-1

	E=computeE(phot%i,phot%j,iT)
	Emax=computeE(phot%i,phot%j,iTmax)
	Emin=computeE(phot%i,phot%j,1)

	if(IsNaN(E)) then
		determineT=3d0
		return
	endif

	if(E.gt.Emax) then
		iT=TMAX-1
		determineT=real(iT)*dT
		return
	endif

	if(E.lt.Emin) then
		iT=1
		determineT=real(iT)*dT
		return
	endif

	iT0=iT
	do while(abs(iTmax-iTmin).gt.1)
		iT=(E1/E)**(0.25)*iT
		if(iT.eq.iT0) then
			if(E1.lt.E) iT=iT0-1
			if(E1.gt.E) iT=iT0+1
		endif
1		continue
		if(iT.le.iTmin) then
			iT=iTmin+1
			goto 1
		endif
		if(iT.ge.iTmax) then
			iT=iTmax-1
			goto 1
		endif
		E=computeE(phot%i,phot%j,iT)
		if(E.ge.E1) then
			iTmax=iT
			Emax=E
		endif
		if(E.le.E1) then
			iTmin=iT
			Emin=E
		endif
		iT0=iT
	enddo

	determineT=(real(iTmin)**4+(real(iTmax)**4-real(iTmin)**4)*(E1-Emin)/(Emax-Emin))**(0.25d0)*dT

	return
	end


	real*8 function determineTslow(phot)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 E1,kp0,kp1
	integer i,ii,iopac

	E1=phot%E
	kp0=0d0
	do ii=1,ngrains
		if(.not.Grain(ii)%qhp) then
			do iopac=1,Grain(ii)%nopac
				kp0=kp0+(Grain(ii)%Kp(iopac,0)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac))
			enddo
		endif
	enddo
	do i=0,TMAX-1
		kp1=0d0
		do ii=1,ngrains
			if(.not.Grain(ii)%qhp) then
				do iopac=1,Grain(ii)%nopac
					kp1=kp1+(Grain(ii)%Kp(iopac,i+1)*C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac))
				enddo
			endif
		enddo
		if(kp0.le.E1.and.kp1.ge.E1) then
			determineTslow=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dT
			return
		endif
		kp0=kp1
	enddo
	determineTslow=real(TMAX-1)*dT
	return
	end



	real*8 function computeEP(i,j,iT,ii)
	use Parameters
	IMPLICIT NONE
	integer iT,ii,iopac,i,j

	computeEP=0d0
	do iopac=1,Grain(ii)%nopac
		computeEP=computeEP+(Grain(ii)%Kp(iopac,iT)*C(i,j)%wopac(ii,iopac))
	enddo

	return
	end


	real*8 function determineTP(phot,ii)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 E1,E,T,Emin,Emax,computeEP
	integer i,ii,iopac,iTmin,iTmax,iT0,iT

	E1=phot%E
	iTmin=0
	iTmax=TMAX

	iT=C(phot%i,phot%j)%TP(ii)/dT
	if(iT.lt.1) iT=1
	if(iT.gt.TMAX-1) iT=TMAX-1

	E=computeEP(phot%i,phot%j,iT,ii)
	Emax=computeEP(phot%i,phot%j,iTmax,ii)
	Emin=0d0

	if(E.gt.Emax) then
		iT=TMAX-1
		determineTP=real(iT)*dT
		return
	endif

	iT0=iT
	do while(abs(iTmax-iTmin).gt.1)
		iT=(E1/E)**(0.25)*iT
		if(iT.eq.iT0) then
			if(E1.lt.E) iT=iT0-1
			if(E1.gt.E) iT=iT0+1
		endif
1		continue
		if(iT.le.iTmin) then
			iT=iTmin+1
			goto 1
		endif
		if(iT.ge.iTmax) then
			iT=iTmax-1
			goto 1
		endif
		E=computeEP(phot%i,phot%j,iT,ii)
		if(E.ge.E1) then
			iTmax=iT
			Emax=E
		endif
		if(E.le.E1) then
			iTmin=iT
			Emin=E
		endif
		iT0=iT
	enddo

	determineTP=(real(iTmin)**4+(real(iTmax)**4-real(iTmin)**4)*(E1-Emin)/(Emax-Emin))**(0.25d0)*dT

	return
	end



	real*8 function determineTPslow(phot,ii)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 E1,kp0,kp1
	integer i,ii,iopac

	E1=phot%E
	kp0=0d0
	do iopac=1,Grain(ii)%nopac
		kp0=kp0+Grain(ii)%Kp(iopac,0)*C(phot%i,phot%j)%wopac(ii,iopac)
	enddo
	do i=0,TMAX-1
		kp1=0d0
		do iopac=1,Grain(ii)%nopac
			kp1=kp1+Grain(ii)%Kp(iopac,i+1)*C(phot%i,phot%j)%wopac(ii,iopac)
		enddo
		if(kp0.le.E1.and.kp1.ge.E1) then
			determineTPslow=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dT
			return
		endif
		kp0=kp1
	enddo
	determineTPslow=real(TMAX-1)*dT
	return
	end




	subroutine emit(phot,spec,Lttot)
	use Parameters
	IMPLICIT NONE
	real*8 spec(nlam),Lr,Lt,ran2,Ltold,Lttot,x,y,z,r
	integer i,ii,iopac
	type(Photon) phot

	nemit=nemit+1
	
	if(multiwav) then
		specemit=spec
		column(1:ngrains,1:ngrains2)=0d0
		if(exportprodimo) Kext_column=0d0
	endif

	call randomdirection(phot%vx,phot%vy,phot%vz)
	phot%scatt=.false.
	phot%pol=.false.

	x=-phot%vy
	y=phot%vx
	z=0d0
	r=sqrt(x**2+y**2+z**2)
	phot%Sx=x/r
	phot%Sy=y/r
	phot%Sz=z/r

	Lr=ran2(idum)*Lttot

	Ltold=0d0
	Lt=0d0
	do i=1,nlam-1
		Lt=Lt+(spec(i)+(spec(i+1)-spec(i))/2d0)*dnu(i)
		if(Lt.ge.Lr) then
			phot%wl1=(Lr-Ltold)/(Lt-Ltold)
			phot%wl2=1d0-phot%wl1
			phot%nu=nu(i)*phot%wl1+nu(i+1)*phot%wl2
			phot%lam=2.9979d14/phot%nu
			phot%ilam1=i
			phot%ilam2=i+1
			goto 1
		endif
		Ltold=Lt
	enddo
	phot%wl1=0d0
	phot%wl2=1d0
	phot%nu=nu(nlam)
	phot%lam=2.9979d14/phot%nu
	phot%ilam1=nlam-1
	phot%ilam2=nlam

1	continue
	do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			Grain(ii)%KabsL(iopac)=Grain(ii)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kabs(iopac,phot%ilam2)*phot%wl2
			Grain(ii)%KextL(iopac)=Grain(ii)%Kext(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kext(iopac,phot%ilam2)*phot%wl2
		enddo
	enddo

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine randomdirection(x,y,z)
	use Parameters
	IMPLICIT NONE
	real*8 x,y,z,r,ran2
	
1	continue
	x=2d0*ran2(idum)-1d0
	y=2d0*ran2(idum)-1d0
	z=2d0*ran2(idum)-1d0
	r=x**2+y**2+z**2
	if(r.gt.1d0) goto 1
	r=sqrt(r)
	x=x/r
	y=y/r
	z=z/r
	
	return
	end


	subroutine ReadTemperatures()
	use Parameters
	IMPLICIT NONE
	character*500 file
	integer i,j,nr,nt,nl,ii,l,iopac
	logical truefalse
	real*8 tot,KappaGas

	if(use_obs_TMC.and.(NphotDiffuse.gt.0.or.dTDiffuse.lt.1d0).and..not.forcediff) then
		if(outputfits) then
			write(file,'(a,"denstempND.fits.gz")') outdir(1:len_trim(outdir))
		else
			write(file,'(a,"denstempND.dat")') outdir(1:len_trim(outdir))
		endif
	else
		if(outputfits) then
			write(file,'(a,"denstemp.fits.gz")') outdir(1:len_trim(outdir))
		else
			write(file,'(a,"denstemp.dat")') outdir(1:len_trim(outdir))
		endif
	endif
	write(*,'("Reading temperature structure from: ",a)') file(1:len_trim(file))
	write(9,'("Reading temperature structure from: ",a)') file(1:len_trim(file))

	call readstruct(file,(/'DENS   ','TEMP   ','COMP   ','GASDENS','DENS0  '/),5,0,.false.)
	do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%TMC=C(i,j)%T
		enddo
	enddo

	if(forcediff) then
		if(outputfits) then
			write(file,'(a,"denstempND.fits.gz")') outdir(1:len_trim(outdir))
		else
			write(file,'(a,"denstempND.dat")') outdir(1:len_trim(outdir))
		endif
		write(*,'("Reading temperature structure from: ",a)') file(1:len_trim(file))
		write(9,'("Reading temperature structure from: ",a)') file(1:len_trim(file))
		call readstruct(file,(/'SKIP   ','TEMPMC '/),2,0,.false.)
	endif


	if(.not.tcontact) then
		do ii=1,ngrains
			if(outputfits) then
				write(file,'(a,"denstempP",i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),ii/10,ii-10*(ii/10)
			else
				write(file,'(a,"denstempP",i1,i1,".dat")') outdir(1:len_trim(outdir)),ii/10,ii-10*(ii/10)
			endif
			write(*,'("Reading temperature structure from: ",a)') file(1:len_trim(file))
			write(9,'("Reading temperature structure from: ",a)') file(1:len_trim(file))
			call readstruct(file,(/'DENSP  ','TEMPP  '/),2,ii,.false.)
		enddo

		do i=1,D%nR-1
			do j=1,D%nTheta-1
				tot=sum(C(i,j)%w(1:ngrains))
				C(i,j)%w(1:ngrains)=C(i,j)%w(1:ngrains)/tot
				C(i,j)%w0(1:ngrains)=C(i,j)%w(1:ngrains)
			enddo
		enddo
	endif

	if(outputfits) then
		write(file,'(a,"Nphotons.fits.gz")') outdir(1:len_trim(outdir))
	else
		write(file,'(a,"Nphotons.dat")') outdir(1:len_trim(outdir))
	endif
	inquire(file=file,exist=truefalse)
	if(.not.truefalse) then
		write(*,'("File ",a," not found")') file(1:len_trim(file))
		write(*,'("Assuming high foton statistics")')
		write(9,'("File ",a," not found")') file(1:len_trim(file))
		write(9,'("Assuming high foton statistics")')
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%Ni=1000000
		enddo
		enddo
	else
		write(*,'("Reading photon statistics from:     ",a)') file(1:len_trim(file))
		write(9,'("Reading photon statistics from:     ",a)') file(1:len_trim(file))
		call readstruct(file,(/'NPHOT  '/),1,0,.false.)
	endif

3	if(outputfits) then
		write(file,'(a,"errors.fits.gz")') outdir(1:len_trim(outdir))
	else
		write(file,'(a,"errors.dat")') outdir(1:len_trim(outdir))
	endif
	inquire(file=file,exist=truefalse)
	if(.not.truefalse) then
		write(*,'("File ",a," not found")') file(1:len_trim(file))
		write(*,'("Assuming perfect temperatures")')
		write(9,'("File ",a," not found")') file(1:len_trim(file))
		write(9,'("Assuming perfect temperatures")')
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%dT=1d-8*C(i,j)%T
		enddo
		enddo
	else
		write(*,'("Reading errors from:                ",a)') file(1:len_trim(file))
		write(9,'("Reading errors from:                ",a)') file(1:len_trim(file))
		call readstruct(file,(/'dTEMP  '/),1,0,.false.)
	endif
	
100	continue

	if(use_qhp) then
		do ii=1,ngrains
			if(Grain(ii)%qhp) then
				if(outputfits) then
					write(file,'(a,"QHPemis",i1,i1,".fits.gz")') outdir(1:len_trim(outdir)),ii/10,ii-10*(ii/10)
				else
					write(file,'(a,"QHPemis",i1,i1,".dat")') outdir(1:len_trim(outdir)),ii/10,ii-10*(ii/10)
				endif
				write(*,'("Reading QHP emission from:          ",a)') file(1:len_trim(file))
				write(9,'("Reading QHP emission from:          ",a)') file(1:len_trim(file))
				call readstruct(file,(/'LAM    ','QHPEJv '/),2,ii,.false.)

				do i=1,D%nR-1
					do j=1,D%nTheta-1
						call integrate(C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam),tot)
						if(tot.gt.1d-200) then
							C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)/tot
						else
							C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=0d0
						endif
					enddo
				enddo
			endif
		enddo
	endif


	if(useTgas) then
		if(outputfits) then
			write(file,'(a,"denstempGas.fits.gz")') outdir(1:len_trim(outdir))
		else
			write(file,'(a,"denstempGas.dat")') outdir(1:len_trim(outdir))
		endif
		write(*,'("Reading gas temperatures from:      ",a)') file(1:len_trim(file))
		write(9,'("Reading gas temperatures from:      ",a)') file(1:len_trim(file))
		call readstruct(file,(/'GASDENS','GASTEMP'/),2,0,.false.)
		do i=1,D%nR-1
			do j=1,D%nTheta-1
				C(i,j)%KappaGas=KappaGas(C(i,j)%gasdens*gas2dust,C(i,j)%Tgas)
			enddo
		enddo
	endif

	return
	end


	subroutine OpticallyThinOld(BW)
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	integer ii,iopac
	real*8 spec(nlam),determineT,determineTP,W,r
	type(photon) phot
c The parameter BW determines if backwarming is assumed.
c When backwarming is assumed, the cooling is a factor of 2 less effective.
	logical BW,BBGrains

	do i=0,D%nR-1
		if(i.ne.0) then
			r=D%R(i)*AU
		else
			r=D%R(1)*AU
		endif
		W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))
		do j=1,D%nTheta-1
			phot%i=i
			phot%j=j
			spec(1:nlam)=0d0
			BBGrains=.false.
			do ii=1,ngrains
				if(.not.Grain(ii)%qhp.and.C(i,j)%w(ii).gt.0d0) then
					do iopac=1,Grain(ii)%nopac
						spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)
     &								*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)/(pi*D%Rstar**2)
					enddo
					BBGrains=.true.
				endif
			enddo
			if(BBGrains) then
				call integrate(spec,phot%E)
				phot%E=W*phot%E
				if(BW) phot%E=phot%E*fBW(1,j)
				C(i,j)%EJv=phot%E
				C(i,j)%T=determineT(phot)
				C(i,j)%TMC=C(i,j)%T
			else
				C(i,j)%T=D%Tstar*W**(0.25)
				C(i,j)%TMC=C(i,j)%T
			endif
			if(.not.tcontact.or.tdes_iter) then
			do ii=1,ngrains
				spec=0d0
				do iopac=1,Grain(ii)%nopac
					spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)/(pi*D%Rstar**2)
				enddo
				call integrate(spec,phot%E)
				phot%E=W*phot%E
				if(BW) phot%E=phot%E*fBW(ii,j)
				C(i,j)%EJvP(ii)=phot%E
				C(i,j)%TP(ii)=determineTP(phot,ii)
			enddo
			endif
			C(i,j)%dT=dT
		enddo
	enddo

	return
	end




	subroutine OpticallyThin(BW)
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	integer ii,iopac
	real*8 determineT,determineTP,W,r,f
	type(photon) phot
c The parameter BW determines if backwarming is assumed.
c When backwarming is assumed, the cooling is a factor of 2 less effective.
	logical BW,BBGrains
	integer iter,iT
	real*8 Tevap,maxT,minT,determinegasfrac,f1,f2,eps
	real*8 Emax,Efrac,A(ngrains),Er,Sig,wtot
	character*500 filename
	real*8 G,mu,ShakuraSunyaevIJ,T
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
	parameter(G=6.67300d-8) ! in cm^3/g/s

	if(.not.allocated(iRfirst)) then
		allocate(iRfirst(1:ngrains,0:D%nTheta))
		iRfirst(1:ngrains,0:D%nTheta)=D%nR-1
	endif

	filename=''
	call tau1heightR(filename)

	write(*,'("Determining optically thin temperatures")')
	write(9,'("Determining optically thin temperatures")')

	do i=0,D%nR-1
		call tellertje(i+1,D%nR)
		if(i.ne.0) then
			r=D%R(i)*AU
		else
			r=D%R(1)*AU
		endif
		W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))
		do j=1,D%nTheta-1
			do ii=1,ngrains
				fBW(ii,j)=1d0
c				Tevap=1500d0
c				minT=0d0
c				maxT=real(TMAX)*dT
c				do iter=1,10
c					if(determinegasfrac(Tevap,i,j,ii).gt.1d0) then
c						maxT=Tevap
c						Tevap=(Tevap+minT)/2d0
c					else
c						minT=Tevap
c						Tevap=(Tevap+maxT)/2d0
c					endif
c				enddo
c				iT=int(Tevap/dT)
c				if(iT.lt.1) iT=1
c				if(iT.gt.TMAX-1)iT=TMAX-1
c				eps=Grain(ii)%Kp(1,iT)/(BBint(iT)*Grain(ii)%Kpabsstar(1))
c				f1=abs((2d0*muRad(j)+1d0/eps)*eps)
c				f2=abs(muRad(j)*(2d0+3d0*muRad(j)*eps)*eps)
c				if(f1.gt.f2) then
c					fBW(ii,j)=f1
c				else
c					fBW(ii,j)=f2
c				endif
c				if(fBW(ii,j).lt.1d0) fBW(ii,j)=1d0
			enddo
				
			do iter=1,3
			phot%i=i
			phot%j=j
			BBGrains=.false.
			phot%E=0d0
			wtot=0d0
			do ii=1,ngrains
				if(.not.Grain(ii)%qhp.and.C(i,j)%w(ii).gt.0d0) then
					do iopac=1,Grain(ii)%nopac
						if(BW) then
							phot%E=phot%E+fBW(ii,j)*Grain(ii)%Kpabsstar(iopac)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
						else
							phot%E=phot%E+Grain(ii)%Kpabsstar(iopac)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
						endif
					enddo
					wtot=wtot+C(i,j)%w(ii)
					BBGrains=.true.
				endif
			enddo
			if(BBGrains) then
				phot%E=W*phot%E*D%Lstar/(pi*D%Rstar**2)/wtot
				if(viscous) then
					T=ShakuraSunyaevIJ(i,j)
					iT=T/dT
					if(iT.gt.TMAX-1)iT=TMAX-1
					if(iT.lt.1)iT=1
					do ii=1,ngrains
						do iopac=1,Grain(ii)%nopac
							phot%E=phot%E+Grain(ii)%Kp(iopac,iT)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
						enddo
					enddo
				endif
				C(i,j)%EJv=phot%E
				C(i,j)%T=determineT(phot)
				C(i,j)%TMC=C(i,j)%T
			else
				C(i,j)%T=D%Tstar*W**(0.25)
				C(i,j)%TMC=C(i,j)%T
			endif
			if(.not.tcontact.or.tdes_iter) then
				do ii=1,ngrains
					phot%E=0d0
					do iopac=1,Grain(ii)%nopac
						phot%E=phot%E+W*Grain(ii)%Kpabsstar(iopac)*C(i,j)%wopac(ii,iopac)*D%Lstar/(pi*D%Rstar**2)
					enddo
					if(BW) phot%E=phot%E*fBW(ii,j)
					if(viscous) then
						do iopac=1,Grain(ii)%nopac
							phot%E=phot%E+Grain(ii)%Kp(iopac,iT)*C(i,j)%wopac(ii,iopac)
						enddo
					endif
					C(i,j)%EJvP(ii)=phot%E
					C(i,j)%TP(ii)=determineTP(phot,ii)
				enddo
			endif
			C(i,j)%dT=dT

			do ii=1,ngrains
				if(.not.tcontact.or.tdes_iter) then
					iT=int(C(i,j)%TP(ii)/dT)
				else
					iT=int(C(i,j)%T/dT)
				endif
				if(iT.lt.1) iT=1
				if(iT.gt.TMAX-1)iT=TMAX-1
				eps=Grain(ii)%Kp(1,iT)/(BBint(iT)*Grain(ii)%Kpabsstar(1))
				f1=abs((2d0*muRad(j)+1d0/eps)*eps)
				f2=abs(muRad(j)*(2d0+3d0*muRad(j)*eps)*eps)
				if(f1.gt.f2) then
					fBW(ii,j)=f1
				else
					fBW(ii,j)=f2
				endif
				if(fBW(ii,j).lt.1d0) fBW(ii,j)=1d0
				if(fBW(ii,j).gt.8d0) fBW(ii,j)=8d0
			enddo

			enddo
			if(computeTgas.or.viscous) C(i,j)%Tgas=C(i,j)%T
		enddo
	enddo

	return
	end



	subroutine BackWarming(fBW0)
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	integer ii,iii,number_invalid,iopac
	real*8 spec(nlam),determineT,determineTP,W,r,f,fBW0,E
	logical BBGrains

	if(.not.allocated(fBW)) allocate(fBW(ngrains,0:D%nTheta))
	if(.not.allocated(iRfirst)) then
		allocate(iRfirst(1:ngrains,0:D%nTheta))
		iRfirst(1:ngrains,0:D%nTheta)=1
	endif

	fBW(1:ngrains,0:D%nTheta)=fBW0
	do j=D%nTheta-1,1,-1
		do iii=1,ngrains
			i=iRfirst(iii,j)
			r=D%R(i)*AU
			W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))
			if(tcontact) then
				spec(1:nlam)=0d0
				BBgrains=.false.
				do ii=1,ngrains
					if(.not.Grain(ii)%qhp.and.C(i,j)%w(ii).gt.0d0) then
						do iopac=1,Grain(ii)%nopac
							spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)
     &									*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)/(pi*D%Rstar**2)
						enddo
						BBGrains=.true.
					endif
				enddo
				if(BBGrains) then
					call integrate(spec,E)
					E=W*E
					f=(C(i,j)%EJv-C(i,j)%dEJv)/E
					if(f.gt.fBW(iii,j).and.f.lt.4d0.and.C(i,j)%Ni.gt.50) then
						fBW(iii,j)=f
					endif
				endif
			else
				if(.not.Grain(iii)%qhp.and.C(i,j)%w(iii).gt.0d0) then
					spec=0d0
					do iopac=1,Grain(iii)%nopac
						spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(iii)%Kabs(iopac,1:nlam)
     &								*C(i,j)%wopac(iii,iopac)/(pi*D%Rstar**2)
					enddo
					call integrate(spec,E)
					E=W*E
					f=(C(i,j)%EJvP(iii)-C(i,j)%dEJv)/E
					if(f.gt.fBW(iii,j).and.f.lt.4d0.and.C(i,j)%Ni.gt.50) then
						fBW(iii,j)=f
					endif
				endif
			endif
			if(j.ne.D%nTheta-1) then
				if(fBW(iii,j).gt.fBW(iii,j+1)) fBW(iii,j)=fBW(iii,j+1)
			endif
		enddo
	enddo

	return
	end


	subroutine RadiationPressure(phot,v,r,z,Kext)
	use Parameters
	IMPLICIT NONE
	real*8 v,r,z,clight,F,Kext
	integer ii
	type(photon) phot
	parameter(clight=2.9979d10) ! cm/s

	if(.not.radpress) return
	
	F=4d0*pi*v*Kext*C(phot%i,phot%j)%dens*AU
	C(phot%i,phot%j)%FradR=C(phot%i,phot%j)%FradR+F*r*phot%E/clight
	C(phot%i,phot%j)%FradZ=C(phot%i,phot%j)%FradZ+F*z*phot%E/clight
	
	return
	end
	
	

	subroutine InteriorIN05()
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	integer ii,iopac
	real*8 spec(nlam),determineT,determineTP,W,r,f1,Kp,Kpabsstar
	type(photon) phot
c The parameter BW determines if backwarming is assumed.
c When backwarming is assumed, the cooling is a factor of 2 less effective.
	logical BW,BBGrains
	integer iter,iT
	real*8 Tevap,maxT,minT,determinegasfrac,f,eps,tau(D%nTheta)
	character*500 filename

	if(.not.allocated(iRfirst)) then
		allocate(iRfirst(1:ngrains,0:D%nTheta))
		iRfirst(1:ngrains,0:D%nTheta)=D%nR-1
	endif

	filename=''
	call tau1heightR(filename)

	tau(1:D%nTheta)=1000d0
	do i=0,D%nR-1
		if(i.ne.0) then
			r=D%R(i)*AU
		else
			r=D%R(1)*AU
		endif
		W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))
		do j=1,D%nTheta-1
			do ii=1,ngrains
				Tevap=1500d0
				minT=0d0
				maxT=real(TMAX)*dT
				do iter=1,100
					if(determinegasfrac(Tevap,i,j,ii).gt.1d0) then
						maxT=Tevap
						Tevap=(Tevap+minT)/2d0
					else
						minT=Tevap
						Tevap=(Tevap+maxT)/2d0
					endif
				enddo
				iT=int(Tevap/dT)
				kp=0d0
				Kpabsstar=0d0
				do iopac=1,Grain(ii)%nopac
					kp=kp+Grain(ii)%Kp(iopac,iT)*C(i,j)%wopac(ii,iopac)
					Kpabsstar=Kpabsstar+Grain(ii)%Kpabsstar(iopac)*C(i,j)%wopac(ii,iopac)
				enddo
				eps=kp/(BBint(iT)*Kpabsstar)
				f1=eps*(muRad(j)*(2d0+3d0*muRad(j)*eps)+(1d0/eps-3d0*eps*muRad(j)**2)*exp(-tau(j)/(muRad(j)*eps)))
				fBW(ii,j)=f1
			enddo
				
			do iter=1,3
			phot%i=i
			phot%j=j
			spec(1:nlam)=0d0
			BBGrains=.false.
			do ii=1,ngrains
				if(.not.Grain(ii)%qhp.and.C(i,j)%w(ii).gt.0d0) then
					do iopac=1,Grain(ii)%nopac
						if(BW) then
							spec(1:nlam)=spec(1:nlam)+fBW(ii,j)*D%Fstar(1:nlam)
     &				*Grain(ii)%Kabs(iopac,1:nlam)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)/(pi*D%Rstar**2)
						else
							spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)
     &				*Grain(ii)%Kabs(iopac,1:nlam)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)/(pi*D%Rstar**2)
						endif
						BBGrains=.true.
					enddo
				endif
			enddo
			if(BBGrains) then
				call integrate(spec,phot%E)
				phot%E=W*phot%E
				C(i,j)%EJv=phot%E
				C(i,j)%T=determineT(phot)
				C(i,j)%TMC=C(i,j)%T
			else
				C(i,j)%T=D%Tstar*W**(0.25)
				C(i,j)%TMC=C(i,j)%T
			endif
			if(.not.tcontact.or.tdes_iter) then
			do ii=1,ngrains
				spec=0d0
				do iopac=1,Grain(ii)%nopac
					spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)*C(i,j)%wopac(ii,iopac)/(pi*D%Rstar**2)
				enddo
				call integrate(spec,phot%E)
				phot%E=W*phot%E
				if(BW) phot%E=phot%E*fBW(ii,j)
				C(i,j)%EJvP(ii)=phot%E
				C(i,j)%TP(ii)=determineTP(phot,ii)
			enddo
			endif
			C(i,j)%dT=dT

			do ii=1,ngrains
				if(.not.tcontact.or.tdes_iter) then
					iT=int(C(i,j)%TP(ii)/dT)
				else
					iT=int(C(i,j)%T/dT)
				endif
				kp=0d0
				Kpabsstar=0d0
				do iopac=1,Grain(ii)%nopac
					kp=kp+Grain(ii)%Kp(iopac,iT)*C(i,j)%wopac(ii,iopac)
					Kpabsstar=Kpabsstar+Grain(ii)%Kpabsstar(iopac)*C(i,j)%wopac(ii,iopac)
				enddo
				eps=kp/(BBint(iT)*Kpabsstar)
				f1=eps*(muRad(j)*(2d0+3d0*muRad(j)*eps)+(1d0/eps-3d0*eps*muRad(j)**2)*exp(-tau(j)/(muRad(j)*eps)))
				fBW(ii,j)=f1
			enddo
			enddo
			do ii=1,ngrains
				do iopac=1,Grain(ii)%nopac
					tau(j)=tau(j)+C(i,j)%dens*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)*Grain(ii)%Kpstar(iopac)*(D%R(i+1)-D%R(i))*AU
				enddo
			enddo
		enddo
	enddo

	return
	end



	subroutine EmitViscous(phot)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 Rad,G,T
	integer iT
	parameter(G=6.67300d-8) ! in cm^3/g/s

	Rad=sqrt(phot%x**2+phot%y**2+phot%z**2)

	T=(3d0*G*D%Mstar*D%Mdot*(1d0-sqrt(D%Rstar/(Rad*AU)))/(8d0*pi*(Rad*AU)**3*sigma))**0.25
	iT=(T/dT+0.5)
	if(iT.lt.1) iT=1
	if(iT.gt.TMAX) iT=TMAX

	call emit(phot,BB(1:nlam,iT),BBint(iT))

	return
	end




	subroutine EmitViscousOld(phot)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 spec(nlam),ran2,T0,T1,increaseT,increaseTP
	real*8 x,y,z,phi,theta,r,Ksca,Kext,kp,Kabs,A(ngrains)
	real*8 epsT1,epsT0,KabsQHP,KabsTot,Ks(nlam),Egas,Efrac
	integer i,iT,iT0,iT1,l,iopac

	spec(1:nlam)=0d0
	kp=0d0
	if(C(phot%i,phot%j)%useFE.and.g2d_heat) then

c absorption
	phot%scatt=.false.
	phot%pol=.false.
	phot%Q=0d0
	phot%U=0d0
	phot%V=0d0

c normal absorption and reemission
	if(.not.tcontact) then
c No thermal contact

	do i=1,ngrains

	if(.not.Grain(i)%qhp) then

	Efrac=phot%E*C(phot%i,phot%j)%FE(i)

	C(phot%i,phot%j)%EabsP(i)=C(phot%i,phot%j)%EabsP(i)+Efrac
	C(phot%i,phot%j)%Eabs=C(phot%i,phot%j)%Eabs+Efrac

	C(phot%i,phot%j)%EJvP(i)=C(phot%i,phot%j)%EJvP(i)+Efrac/(C(phot%i,phot%j)%mass*C(phot%i,phot%j)%w(i))
	C(phot%i,phot%j)%EJv=C(phot%i,phot%j)%EJv+Efrac/(C(phot%i,phot%j)%mass)
	
	T0=C(phot%i,phot%j)%TP(i)
	T1=increaseTP(phot,i)

	iT0=int(T0/dT)
	iT1=int(T1/dT)
	epsT0=T0-real(iT0)*dT
	epsT1=T1-real(iT1)*dT
	if(iT0.ge.TMAX-1) iT0=TMAX-2
	if(iT1.ge.TMAX-1) iT1=TMAX-1
	if(iT0.eq.iT1) then
		do iopac=1,Grain(i)%nopac
			do l=1,nlam
				spec(l)=spec(l)+(BB(l,iT0+1)-BB(l,iT0))*Grain(i)%Kabs(iopac,l)
     &					*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
			enddo
			kp=kp+(Grain(i)%Kp(iopac,iT0+1)-Grain(i)%Kp(iopac,iT0))*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
		enddo
	else
		do iopac=1,Grain(i)%nopac
			do l=1,nlam
				spec(l)=spec(l)+(epsT1*BB(l,iT1+1)+(1d0-epsT1)*BB(l,iT1)-epsT0*BB(l,iT0+1)-(1d0-epsT0)*BB(l,iT0))*
     &				Grain(i)%Kabs(iopac,l)*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
			enddo
			kp=kp+(epsT1*Grain(i)%Kp(iopac,iT1+1)+(1d0-epsT1)*Grain(i)%Kp(iopac,iT1)
     &			-epsT0*Grain(i)%Kp(iopac,iT0+1)-(1d0-epsT0)*Grain(i)%Kp(iopac,iT0))
     &			*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
		enddo
	endif

	C(phot%i,phot%j)%TP(i)=T1
	endif
	enddo

	C(phot%i,phot%j)%T=increaseT(phot)
	
	else
c Thermal contact (single temperature)
	do i=1,ngrains
		Efrac=phot%E*C(phot%i,phot%j)%FE(i)
		C(phot%i,phot%j)%Eabs=C(phot%i,phot%j)%Eabs+Efrac
		C(phot%i,phot%j)%EJv=C(phot%i,phot%j)%EJv+Efrac/(C(phot%i,phot%j)%mass)
	enddo

	T0=C(phot%i,phot%j)%T
	T1=increaseT(phot)

	iT0=int(T0/dT)
	iT1=int(T1/dT)

	epsT0=T0-real(iT0)*dT
	epsT1=T1-real(iT1)*dT
	if(iT0.ge.TMAX-1) iT0=TMAX-2
	if(iT1.ge.TMAX-1) iT1=TMAX-1

	if(iT0.eq.iT1) then
		do l=1,nlam
			Kabs=0d0
			do i=1,ngrains
				if(.not.Grain(i)%qhp) then
					do iopac=1,Grain(i)%nopac
						Kabs=Kabs+Grain(i)%Kabs(iopac,l)*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
					enddo
				endif
			enddo
			spec(l)=(BB(l,iT0+1)-BB(l,iT0))*Kabs
		enddo
		kp=0d0
		do i=1,ngrains
			if(.not.Grain(i)%qhp) then
				do iopac=1,Grain(i)%nopac
					kp=kp+(Grain(i)%Kp(iopac,iT0+1)-Grain(i)%Kp(iopac,iT0))*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
				enddo
			endif
		enddo
	else
		do l=1,nlam
			Kabs=0d0
			do i=1,ngrains
				if(.not.Grain(i)%qhp) then
					do iopac=1,Grain(i)%nopac
						Kabs=Kabs+Grain(i)%Kabs(iopac,l)*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
					enddo
				endif
			enddo
			spec(l)=(epsT1*BB(l,iT1+1)+(1d0-epsT1)*BB(l,iT1)-epsT0*BB(l,iT0+1)-(1d0-epsT0)*BB(l,iT0))*Kabs
		enddo
		kp=0d0
		do i=1,ngrains
			if(.not.Grain(i)%qhp) then
				do iopac=1,Grain(i)%nopac
					kp=kp+(epsT1*Grain(i)%Kp(iopac,iT1+1)+(1d0-epsT1)*Grain(i)%Kp(iopac,iT1)
     &			-epsT0*Grain(i)%Kp(iopac,iT0+1)-(1d0-epsT0)*Grain(i)%Kp(iopac,iT0))
     &			*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
    			enddo
			endif
		enddo
	endif

	if(kp.ne.0d0.and.use_qhp) then
		spec(1:nlam)=spec(1:nlam)*(KabsTot-KabsQHP)/(kp*KabsTot)
		kp=1d0
	endif
	do i=1,ngrains
		if(Grain(i)%qhp) then
			do iopac=1,Grain(i)%nopac
				KabsQHP=(Grain(i)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(i)%Kabs(iopac,phot%ilam2)*phot%wl2)*C(phot%i,phot%j)%w(i)*C(phot%i,phot%j)%wopac(i,iopac)
			enddo
			if(use_qhp) then
				spec(1:nlam)=spec(1:nlam)
     &	+KabsQHP*C(phot%i,phot%j)%QHP(Grain(i)%qhpnr,1:nlam)/KabsTot
			endif
		endif
	enddo

	C(phot%i,phot%j)%T=T1

	endif

	spec=spec/kp
	kp=1d0

	Efrac=phot%E*C(phot%i,phot%j)%FE(0)

	T0=C(phot%i,phot%j)%Tgas
	iT0=int(T0/dT)
	if(iT0.gt.TMAX-1) iT0=TMAX-1
	spec(1:nlam)=spec(1:nlam)*(phot%E-Efrac)+Efrac*BB(1:nlam,iT0)/BBint(iT0)
	kp=phot%E

	else

	T0=C(phot%i,phot%j)%Tgas
	iT0=int(T0/dT)
	if(iT0.lt.1) iT0=1
	if(iT0.gt.TMAX-1) iT0=TMAX-1
	spec(1:nlam)=BB(1:nlam,iT0)/BBint(iT0)
	kp=1d0
	
	endif

	call emit(phot,spec,kp)

	return
	end


	real*8 function KappaGas(dens,T)
	IMPLICIT NONE
	real*8 dens,T,aK
	call gop(dens,T,aK)
	KappaGas=aK
	return
	end
	
	
	subroutine determineTgas(i,j)
	use Parameters
	IMPLICIT NONE
	real*8 H,Tgas,Mdot,nH,nd,mH,T,abar,KappaGas,A,G,mu,diff,mindiff,Sig
	real*8 Evisc,Egas
	parameter(mH=1.67d-24)
	parameter(G=6.67300d-8) ! in cm^3/g/s
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
	integer ii,iopac,i,j,iT,jj,iter,i2
	real*8 Tmaxi,Tmini,W,r,Eabs,Kappa,tau0

	if(D%R_av(i)/AU.gt.HotGasMinRad.and.D%R_av(i)/AU.lt.HotGasMaxRad) then
		Tgas=HotGasT
		goto 2
	endif

	if(deadzone.or.gravstable.or.D%Rexp.lt.1d10) then
	   Mdot=D%MdotR(i)
	else
	   Mdot=D%Mdot
	endif

	Evisc=(2d0*sqrt(D%Rstar/(AU*D%R(i+1)))-3d0)/(3d0*D%R(i+1)*AU)
	Evisc=Evisc-(2d0*sqrt(D%Rstar/(AU*D%R(i)))-3d0)/(3d0*D%R(i)*AU)
	Evisc=2d0*pi*Evisc*3d0*G*D%Mstar*Mdot/(4d0*pi)
	Evisc=Evisc/(4d0*pi)

	Sig=0d0
	do jj=1,D%nTheta-1
		Sig=Sig+C(i,jj)%gasdens*C(i,jj)%V
	enddo
	Evisc=Evisc*C(i,j)%gasdens/Sig

	tau0=0d0
	do i2=1,i-1
		do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				tau0=tau0+(D%R(i+1)-D%R(i))*AU*C(i2,j)%dens*C(i2,j)%w(ii)*Grain(ii)%Kpstar(iopac)*C(i2,j)%wopac(ii,iopac)
			enddo
		enddo
	enddo

	if(i.ne.0) then
		r=D%R(i)*AU
	else
		r=D%R(1)*AU
	endif
	W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))*exp(-tau0)

	nH=C(i,j)%gasdens*gas2dust/mH
	Tgas=C(i,j)%T
	Tmaxi=TMAX*dT
	Tmini=1d0
c=======================================================================
c set the minimum gas temperature to the dust temperature (02-09-2011)
c=======================================================================
	Tmini=C(i,j)%T
c=======================================================================
c
c=======================================================================
	if(Tgas.le.Tmini.or.Tgas.ge.Tmaxi) Tgas=(Tmaxi+Tmini)/2d0
1	continue
		Kappa=KappaGas(C(i,j)%gasdens*gas2dust,D%Tstar)
		Eabs=Kappa*W*D%Lstar/(pi*D%Rstar**2)
		Eabs=Eabs*C(i,j)%gasdens*gas2dust

		Egas=Evisc+C(i,j)%Egas*C(i,j)%dens/C(i,j)%V!+Eabs
		H=0d0
		iT=int(Tgas/dT)
		do ii=1,ngrains
			nd=C(i,j)%dens*C(i,j)%w(ii)/(4d0*pi*Grain(ii)%rv**3*Grain(ii)%rho/3d0)
			if(.not.Grain(ii)%qhp) then
			if(tcontact) then
				T=C(i,j)%T
			else
				T=C(i,j)%TP(ii)
			endif
			abar=0.35d0*exp(-sqrt((T+Tgas)/500d0))+0.1d0
			H=H+nH*nd*pi*Grain(ii)%rv**2*sqrt(8d0*kb*Tgas/(pi*mH))*abar*2d0*kb*(T-Tgas)
			endif
		enddo
		H=H/(4d0*pi)
		Kappa=KappaGas(C(i,j)%gasdens*gas2dust,Tgas)
		diff=Egas+H-BBint(iT)*Kappa*C(i,j)%gasdens*gas2dust
		if(diff.lt.0d0) then
			Tmaxi=Tgas
			Tgas=(Tgas+Tmini)/2d0
			if(abs(Tmaxi-Tmini).gt.1d0) goto 1
		else
			Tmini=Tgas
			Tgas=(Tgas+Tmaxi)/2d0
			if(abs(Tmaxi-Tmini).gt.1d0) goto 1
		endif
2	continue
	
	C(i,j)%Tgas=Tgas
	iT=int(Tgas/dT)

	C(i,j)%KappaGas=KappaGas(C(i,j)%gasdens*gas2dust,C(i,j)%Tgas)
	if(.not.allocated(C(i,j)%FE)) allocate(C(i,j)%FE(0:ngrains))
	C(i,j)%useFE=.true.
	if(iT.lt.TMAX) then
		C(i,j)%FE(0)=BBint(iT)*C(i,j)%KappaGas/(BBint(iT)*C(i,j)%KappaGas+C(i,j)%EJv)
	else
		C(i,j)%FE(0)=BBint(TMAX)*C(i,j)%KappaGas/(BBint(TMAX)*C(i,j)%KappaGas+C(i,j)%EJv)
	endif

	A=0d0
	do ii=1,ngrains
		A=A+C(i,j)%w(ii)/Grain(ii)%rv
	enddo
	do ii=1,ngrains
		C(i,j)%FE(ii)=((1d0-C(i,j)%FE(0))*C(i,j)%w(ii)/Grain(ii)%rv)/A
		if(viscous) C(i,j)%EviscDirect(ii)=Evisc*C(i,j)%FE(ii)*C(i,j)%V
	enddo
	
	return
	end

	subroutine MakeAdiabatic(nabla_ad)
	use Parameters
	IMPLICIT NONE
	real*8 nabla_ad,nabla
	integer i,j
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Correcting the temperature structure for convection")')
	write(9,'("Correcting the temperature structure for convection")')

	do i=1,D%nR-1
	do j=2,D%nTheta-1
		nabla=log(C(i,j-1)%T/C(i,j)%T)/log(C(i,j-1)%gasdens*C(i,j-1)%T/(C(i,j)%gasdens*C(i,j)%T))
		if(nabla.gt.nabla_ad) then
			C(i,j)%T=(C(i,j-1)%T*(C(i,j)%gasdens/(C(i,j-1)%gasdens*C(i,j-1)%T))**nabla_ad)**(1d0/(1d0-nabla_ad))
		endif
	enddo
	enddo
	
	return
	end
	

	module InnerGasDiskEmission
		integer NRAD
		parameter(NRAD=250)
		real*8 R(NRAD),Ftot(NRAD)
	end module InnerGasDiskEmission


	subroutine InitInnerGasDisk(lam0,Etot)
	use InnerGasDiskEmission
	use Parameters
	IMPLICIT NONE
	integer i
	real*8 lam0,T,G,Planck,F(NRAD),Etot,Rad
	parameter(G=6.67300d-8) ! in cm^3/g/s

	do i=1,NRAD
		R(i)=D%R(0)*Rinner_gas+(D%R(1)-D%R(0)*Rinner_gas)*real(i-1)/real(NRAD-1)
	enddo
	Ftot(1)=0d0
	do i=1,NRAD-1
		Rad=sqrt(R(i)*R(i+1))*AU
		T=(3d0*G*D%Mstar*D%Mdot*(1d0-sqrt(D%Rstar/(Rad)))/(8d0*pi*(Rad)**3*sigma))**0.25
		F(i)=Planck(T,lam0)/4d0
		Ftot(i+1)=Ftot(i)+AU**2*pi*(R(i+1)**2-R(i)**2)*F(i)
	enddo
	Etot=Ftot(NRAD)

	return
	end

	real*8 function RadInnerGasDisk(E)
	use InnerGasDiskEmission
	use Parameters
	IMPLICIT NONE
	integer i
	real*8 E,eps
	
	call hunt(Ftot,NRAD,E,i)

	eps=(E-Ftot(i))/(Ftot(i+1)-Ftot(i))
	RadInnerGasDisk=R(i)+eps*(R(i+1)-R(i))

	return
	end



