c--------------------------------------------------------------------
c This subrotuine regrids the radial grid around R=R1 and sets the
c edge of the disk to R=R2
c 2007-10-09: Improved interpolation for the regridding of the density
c 2007-10-11: Regrid using the optical depth. 
c             Finest grid is from tau=0.1-10
c 2008-01-22: Grid refinement improved for stability
c 2008-02-08: Grid refinement now included for TdiskStruct
c 2008-02-08: Changed R2 to be the limitting radius of refinement
c--------------------------------------------------------------------
	subroutine RegridR(R1,R2)
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	real*8 Rnew(0:D%nR),R1,R2,Rav(0:D%nR),tot,tot0
	real*8 MassTot,MassTot0
	real*8 V,M,M0,Ni,determineT,determineTP
	type(cell) Cold(0:D%nR)
	type(photon) phot
	integer k,i1,i2,ii,itot,jj,l,inext,nrefine,nused,iopac
	character*100 surfdensfile
	real*8 tau,tau0,lam0,Kext,dens,w(ngrains),taumax,R
	real*8 taunext,dltau,RR1,lR1,lR2,totM0,Mg
	real*8 taustart,dtaumax,tauend,ct,Rmultmax,dtaumaxabs
	logical escape,hitstar,locfield,RNDWBackup
	integer iT,irefine,inormal
	logical refineR,refining,error
	real*8 taustar(0:D%nTheta),taulocal(0:D%nTheta),shtemp(0:D%nR)
	real*8 taustarabs(0:D%nTheta),taulocalabs(0:D%nTheta),MaxGasMass(ngrains)
	real*8 Kpext(ngrains,0:TMAX),Kpabs(ngrains,0:TMAX),Kpstarabs(ngrains)
	real*8 Rnext(D%nTheta),Rprev,R1refine,TdiskStruct(0:D%nR),spec(nlam)

	if(.not.allocated(Temp)) then
		allocate(Temp(0:D%nR,0:D%nTheta))
		Temp(0:D%nR,0:D%nTheta)=C(0:D%nR,0:D%nTheta)%T
	endif

	MassTot0=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		MassTot0=MassTot0+C(i,j)%dens0*C(i,j)%V
	enddo
	enddo

	if(R1.lt.D%Rin) R1=D%Rin

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Refining radial grid around:  ",f10.3," AU")') R1
	write(9,'("Refining radial grid around:  ",f10.3," AU")') R1

	Rnew(0)=D%R(0)
	Rnew(1)=D%Rin

c regridding with the optical depth

	do ii=1,ngrains
		Kpext(ii,0:TMAX)=0d0
		Kpabs(ii,0:TMAX)=0d0
		do j=1,TMAX
			call integrate(BB(1:nlam,j),tot)
			spec(1:nlam)=BB(1:nlam,j)*Grain(ii)%Kext(1,1:nlam)
			!! assumed iopac=1 (no loop)
			call integrate(spec,Kpext(ii,j))
			Kpext(ii,j)=Kpext(ii,j)/tot
			spec(1:nlam)=BB(1:nlam,j)*Grain(ii)%Kabs(1,1:nlam)
			!! assumed iopac=1 (no loop)
			call integrate(spec,Kpabs(ii,j))
			Kpabs(ii,j)=Kpabs(ii,j)/tot
		enddo
		spec(1:nlam)=D%Fstar(1:nlam)*Grain(ii)%Kabs(1,1:nlam)
		!! assumed iopac=1 (no loop)
		call integrate(spec,Kpstarabs(ii))
		Kpstarabs(ii)=Kpstarabs(ii)/D%Lstar
	enddo

	RNDWBackup=RNDW
	RNDW=.false.

	locfield=.false.
	taustar(1:D%nTheta-1)=0d0
	taulocal(1:D%nTheta-1)=0d0
	taustar(0)=1d0
	taulocal(0)=1d0
	taustar(D%nTheta)=1d0
	taulocal(D%nTheta)=1d0

	taustarabs(1:D%nTheta-1)=0d0
	taulocalabs(1:D%nTheta-1)=0d0
	taustarabs(0)=1d0
	taulocalabs(0)=1d0
	taustarabs(D%nTheta)=1d0
	taulocalabs(D%nTheta)=1d0

	taustart=0.1d0
	dtaumax=0.5d0
	tauend=10d0
	dtaumaxabs=200d0

	nrefine=nlev*nspan-nspan
	irefine=0
	inormal=0
	R1refine=Rnew(1)
	nused=1
	
	do i=1,D%nR-D%nRfix-nspan-3
		do jj=1,D%nR-1
			if(Rnew(i).ge.D%R(jj).and.Rnew(i).lt.D%R(jj+1)) goto 1
		enddo
		jj=1
1		continue
		Rmultmax=(D%Rout/Rnew(i))**(1d0/real(D%nR-D%nRfix-nspan-irefine-inormal-3))
		if(Rmultmax.le.1d0) Rmultmax=1.001d0
		Rnew(i+1)=Rnew(i)*Rmultmax
		if(i.eq.1.and.tdes_iter) Rnew(i+1)=Rbwdes(D%nTheta-1)
		if(Rnew(i+1).le.Rnew(i)) Rnew(i+1)=Rnew(i)*Rmultmax
		if(Rnew(i+1).ge.D%R(D%nR-1)) Rnew(i+1)=(D%R(D%nR-1)+Rnew(i))/2d0
		refineR=.false.
		refining=.false.
		do j=1,D%nTheta-1
			if(taustar(j).lt.tauend) then
c     &     .or.taustar(j-1).lt.taustart.or.taustar(j+1).lt.taustart) then
				if(taustar(j).gt.taustart) then
					tau0=taustar(j)*dtaumax
c					if(tau0.gt.dtaumaxabs.and.
c     &	(taustar(j-1).lt.taustart.or.taustar(j+1).lt.taustart)) tau0=dtaumaxabs
				else
					tau0=taustart
				endif
				k=jj
				Rprev=Rnew(i)
10				continue
				Kext=0d0
				do ii=1,ngrains
				   do iopac=1,Grain(ii)%nopac
				      Kext=Kext+C(k,j)%w(ii)*C(k,j)%wopac(ii,iopac)*Grain(ii)%Kpstar(iopac)
				   enddo
				enddo
				tau=Kext*AU*C(k,j)%dens*(D%R(k+1)-Rprev)
				if(tau.ge.tau0.or.k.ge.D%nR-1) then
					Rnext(j)=Rprev+(D%R(k+1)-Rprev)*tau0/tau
					goto 20
				endif
				tau0=tau0-tau
				Rprev=D%R(k+1)
				k=k+1
				goto 10
20				continue
				if(Rnext(j).lt.D%R(D%nR-1)) refining=.true.
				if(Rnext(j).lt.Rnew(i+1).and.Rnext(j).gt.Rnew(i)) then
					Rnew(i+1)=Rnext(j)
					inext=k
					refineR=.true.
				endif
			endif
			if(taulocal(j).lt.tauend) then
c     &     .or.taulocal(j-1).lt.taustart.or.taulocal(j+1).lt.taustart) then
				if(taulocal(j).gt.taustart) then
					tau0=taulocal(j)*dtaumax
c					if(tau0.gt.dtaumaxabs.and.
c     &	(taulocal(j-1).lt.taustart.or.taulocal(j+1).lt.taustart)) tau0=dtaumaxabs
				else
					tau0=taustart
				endif
				k=jj
				Rprev=Rnew(i)
30				continue
				Kext=0d0
				iT=int(C(k,j)%T/dT)
				if(iT.gt.TMAX) iT=TMAX
				do ii=1,ngrains
					Kext=Kext+C(k,j)%w(ii)*Kpext(ii,iT)
				enddo
				tau=Kext*AU*C(k,j)%dens*(D%R(k+1)-Rprev)
				if(tau.ge.tau0.or.k.ge.D%nR-1) then
					Rnext(j)=Rprev+(D%R(k+1)-Rprev)*tau0/tau
					goto 40
				endif
				tau0=tau0-tau
				Rprev=D%R(k+1)
				k=k+1
				goto 30
40				continue
				if(Rnext(j).lt.D%R(D%nR-1)) refining=.true.
				if(Rnext(j).lt.Rnew(i+1).and.Rnext(j).gt.Rnew(i)) then
					Rnew(i+1)=Rnext(j)
					inext=k
					refineR=.true.
				endif
			endif
		enddo
		if(refineR) then
			if(irefine.eq.0) R1refine=Rnew(i+1)
			irefine=irefine+1
			nused=i+1
			if(irefine.ge.nrefine) goto 2
		else
			inormal=inormal+1
		endif
		if(Rnew(i+1).ge.R2) then
			Rnew(i+1)=R2
			goto 2
		endif
		do inext=1,D%nR-1
			if(Rnew(i+1).ge.D%R(inext).and.Rnew(i+1).lt.D%R(inext+1)) goto 3
		enddo
		inext=D%nR-1
3		continue
		if(.not.refining) goto 2
		do j=1,D%nTheta-1
			if(taustar(j).lt.tauend) then
				if(jj.eq.inext) then
					Kext=0d0
					do ii=1,ngrains
					   do iopac=1,Grain(ii)%nopac
					      Kext=Kext+C(jj,j)%w(ii)*C(jj,j)%wopac(ii,iopac)*Grain(ii)%Kpstar(iopac)
					   enddo
					enddo
					taustar(j)=taustar(j)+(Rnew(i+1)-Rnew(i))*Kext*C(jj,j)%dens*AU
				else
					Kext=0d0
					do ii=1,ngrains
					   do iopac=1,Grain(ii)%nopac
					      Kext=Kext+C(jj,j)%w(ii)*C(jj,j)%wopac(ii,iopac)*Grain(ii)%Kpstar(iopac)
					   enddo
					enddo
					taustar(j)=taustar(j)+(D%R(jj+1)-Rnew(i))*Kext*C(jj,j)%dens*AU
					do l=jj+1,inext-1
						Kext=0d0
						do ii=1,ngrains
						   do iopac=1,Grain(ii)%nopac
						      Kext=Kext+C(l,j)%w(ii)*C(l,j)%wopac(ii,iopac)*Grain(ii)%Kpstar(iopac)
						   enddo
						enddo
						taustar(j)=taustar(j)+(D%R(l+1)-D%R(l))*Kext*C(l,j)%dens*AU
					enddo
					Kext=0d0
					do ii=1,ngrains
					   do iopac=1,Grain(ii)%nopac
					      Kext=Kext+C(inext,j)%w(ii)*C(inext,j)%wopac(ii,iopac)*Grain(ii)%Kpstar(iopac)
					   enddo
					enddo
					taustar(j)=taustar(j)+(Rnew(i+1)-D%R(inext))*Kext*C(inext,j)%dens*AU
				endif
			endif
			if(taulocal(j).lt.tauend) then
				if(jj.eq.inext) then
					Kext=0d0
					iT=int(C(jj,j)%T/dT)
					if(iT.gt.TMAX) iT=TMAX
					do ii=1,ngrains
						Kext=Kext+C(jj,j)%w(ii)*Kpext(ii,iT)
					enddo
					taulocal(j)=taulocal(j)+(Rnew(i+1)-Rnew(i))*Kext*C(jj,j)%dens*AU
				else
					Kext=0d0
					iT=int(C(jj,j)%T/dT)
					if(iT.gt.TMAX) iT=TMAX
					do ii=1,ngrains
						Kext=Kext+C(jj,j)%w(ii)*Kpext(ii,iT)
					enddo
					taulocal(j)=taulocal(j)+(D%R(jj+1)-Rnew(i))*Kext*C(jj,j)%dens*AU
					do l=jj+1,inext-1
						Kext=0d0
						iT=int(C(l,j)%T/dT)
						if(iT.gt.TMAX) iT=TMAX
						do ii=1,ngrains
							Kext=Kext+C(l,j)%w(ii)*Kpext(ii,iT)
						enddo
						taulocal(j)=taulocal(j)+(D%R(l+1)-D%R(l))*Kext*C(l,j)%dens*AU
					enddo
					Kext=0d0
					iT=int(C(inext,j)%T/dT)
					if(iT.gt.TMAX) iT=TMAX
					do ii=1,ngrains
						Kext=Kext+C(inext,j)%w(ii)*Kpext(ii,iT)
					enddo
					taulocal(j)=taulocal(j)+(Rnew(i+1)-D%R(inext))*Kext*C(inext,j)%dens*AU
				endif
			endif
		enddo
	enddo

2	continue

	call sort(Rnew(1:nused),nused)

	R2=Rnew(nused)

	write(*,'("Refined region:   ",f10.3," -",f10.3," AU")') R1refine,Rnew(nused)
	write(9,'("Refined region:   ",f10.3," -",f10.3," AU")') R1refine,Rnew(nused)

	lR1=log10(Rnew(nused))
	lR2=log10(D%Rout)

	if(tdes_iter.or.Rnew(2).gt.(2d0*Rnew(1))) then
		do i=nused+1,nused+nspan
			Rnew(i)=Rnew(1)+(Rnew(2)-Rnew(1))*(real(i-nused)/real(nspan+1))**0.25
		enddo
		nused=nused+nspan
	endif

	do i=1,D%nRfix
		if(nused.lt.D%nR) then
			nused=nused+1
			Rnew(nused)=D%Rfix(i)
		endif
	enddo

C	do i=1,D%nRfix
C		phot%i=-1
C		phot%j=D%nTheta-1
C		r=sqrt(phot%x**2+phot%y**2+phot%z**2)
C		do j=0,D%nR-1
C			if(D%Rfix(i).gt.D%R(j).and.D%Rfix(i).le.D%R(j+1)) then
C				phot%i=j
C				r=D%R_av(j)/AU
C				phot%x=r*sin(D%theta_av(phot%j))*cos(0.1)
C				phot%y=r*sin(D%theta_av(phot%j))*sin(0.1)
C				phot%z=r*cos(D%theta_av(phot%j))
C			endif
C		enddo
C		if(phot%i.eq.-1) then
C			phot%i=0
C			r=D%R_av(phot%i)/AU
C			phot%x=r*sin(D%theta_av(phot%j))*cos(0.1)
C			phot%y=r*sin(D%theta_av(phot%j))*sin(0.1)
C			phot%z=r*cos(D%theta_av(phot%j))
C		endif
C		phot%vx=phot%x/r
C		phot%vy=phot%y/r
C		phot%vz=phot%z/r
C		r=sqrt(phot%vx**2+phot%vy**2+phot%vz**2)
C		phot%vx=phot%vx/r
C		phot%vy=phot%vy/r
C		phot%vz=phot%vz/r
C
C		tau=0.25
C		phot%lam=0.55
C		phot%onEdge=.false.
C		do j=1,nlam-1
C			if(phot%lam.ge.lam(j).and.phot%lam.le.lam(j+1)) then
C				phot%ilam1=j
C				phot%ilam2=j+1
C				phot%wl1=(lam(j+1)-phot%lam)/(lam(j+1)-lam(j))
C				phot%wl2=(phot%lam-lam(j))/(lam(j+1)-lam(j))
C			endif
C		enddo
C		phot%nu=nu(phot%ilam1)*phot%wl1+nu(phot%ilam2)*phot%wl2
C		phot%E=1d0
C		do ii=1,ngrains
C			do iopac=1,Grain(ii)%nopac
C				Grain(ii)%KabsL(iopac)=Grain(ii)%Kabs(iopac,phot%ilam1)*phot%wl1+
C     &     Grain(ii)%Kabs(iopac,phot%ilam2)*phot%wl2
C				Grain(ii)%KextL(iopac)=Grain(ii)%Kext(iopac,phot%ilam1)*phot%wl1+
C     &     Grain(ii)%Kext(iopac,phot%ilam2)*phot%wl2
C			enddo
C		enddo
C		do j=1,nspan*nlev/2
C			call trace2d(phot,tau,escape,hitstar,.false.)
C			if(escape.or.hitstar) goto 91
C			tau=tau*2d0
C			if(nused.lt.D%nR) then
C				nused=nused+1
C				Rnew(nused)=sqrt(phot%x**2+phot%y**2+phot%z**2)
C				if(Rnew(nused).le.D%Rin) nused=nused-1
C			endif
C		enddo
C91		continue
C
C		phot%i=-1
C		phot%j=D%nTheta-1
C		r=sqrt(phot%x**2+phot%y**2+phot%z**2)
C		do j=0,D%nR-1
C			if(D%Rfix(i).gt.D%R(j).and.D%Rfix(i).le.D%R(j+1)) then
C				phot%i=j
C				r=D%R_av(j)/AU
C				phot%x=r*sin(D%theta_av(phot%j))*cos(0.1)
C				phot%y=r*sin(D%theta_av(phot%j))*sin(0.1)
C				phot%z=r*cos(D%theta_av(phot%j))
C			endif
C		enddo
C		if(phot%i.le.0) then
C			phot%i=1
C			r=D%R_av(phot%i)/AU
C			phot%x=r*sin(D%theta_av(phot%j))*cos(0.1)
C			phot%y=r*sin(D%theta_av(phot%j))*sin(0.1)
C			phot%z=r*cos(D%theta_av(phot%j))
C		endif
C		phot%vx=-phot%x/r
C		phot%vy=-phot%y/r
C		phot%vz=phot%z/r
C		r=sqrt(phot%vx**2+phot%vy**2+phot%vz**2)
C		phot%vx=phot%vx/r
C		phot%vy=phot%vy/r
C		phot%vz=phot%vz/r
C		tau=0.25
C		phot%lam=0.55
C		phot%onEdge=.false.
C		do j=1,nlam-1
C			if(phot%lam.ge.lam(j).and.phot%lam.le.lam(j+1)) then
C				phot%ilam1=j
C				phot%ilam2=j+1
C				phot%wl1=(lam(j+1)-phot%lam)/(lam(j+1)-lam(j))
C				phot%wl2=(phot%lam-lam(j))/(lam(j+1)-lam(j))
C			endif
C		enddo
C		phot%nu=nu(phot%ilam1)*phot%wl1+nu(phot%ilam2)*phot%wl2
C		phot%E=1d0
C		do ii=1,ngrains
C			do iopac=1,Grain(ii)%nopac
C				Grain(ii)%KabsL(iopac)=Grain(ii)%Kabs(iopac,phot%ilam1)*phot%wl1+
C     &     Grain(ii)%Kabs(iopac,phot%ilam2)*phot%wl2
C				Grain(ii)%KextL(iopac)=Grain(ii)%Kext(iopac,phot%ilam1)*phot%wl1+
C     &     Grain(ii)%Kext(iopac,phot%ilam2)*phot%wl2
C			enddo
C		enddo
C		do j=1,nspan*nlev/2
C			call trace2d(phot,tau,escape,hitstar,.false.)
C			if(escape.or.hitstar.or.phot%x.lt.0d0) goto 92
C			tau=tau*2d0
C			if(nused.lt.D%nR) then
C				nused=nused+1
C				Rnew(nused)=sqrt(phot%x**2+phot%y**2+phot%z**2)
C				if(Rnew(nused).le.D%Rin) nused=nused-1
C			endif
C		enddo
C92		continue
C	enddo

	if(nused.lt.D%nR) then
		do i=nused+1,D%nR
			Rnew(i)=10d0**(lR1+(lR2-lR1)*real(i-nused)/real(D%nR-nused))
		enddo
	endif
	
	RNDW=RNDWBackup
70	continue
	error=.false.
	do i=1,D%nR-1
		if(Rnew(i-1).eq.Rnew(i)) then
			Rnew(i)=(Rnew(i+1)+Rnew(i))/2d0
			error=.true.
		endif
	enddo
	if(error) goto 70

	call sort(Rnew(0:D%nR),D%nR+1)
	Rnew(D%nR)=D%Rout

71	continue
	if(Rnew(1).eq.Rnew(2)) then
		Rnew(2)=(Rnew(2)+Rnew(3))/2d0
		goto 71
	endif
	do i=1,D%nR-1
		if(Rnew(i).eq.Rnew(i+1)) then
			Rnew(i)=(Rnew(i-1)+Rnew(i))/2d0
			goto 71
		endif
	enddo
	
	do i=0,D%nR-1
		Rav(i)=AU*10d0**((log10(Rnew(i))+log10(Rnew(i+1)))/2d0)
	enddo

	do i=0,D%nR
		allocate(Cold(i)%w(ngrains))
		allocate(Cold(i)%w0(ngrains))
		if(.not.tcontact.or.tdes_iter) then
			allocate(Cold(i)%TP(ngrains))
			allocate(Cold(i)%EJvP(ngrains))
		endif
		if(use_qhp) then
			allocate(Cold(i)%LRF(nlam))
			allocate(Cold(i)%QHP(nqhp,nlam))
			allocate(Cold(i)%tdistr(nqhp,NTQHP))
			allocate(Cold(i)%Tqhp(nqhp))
			allocate(Cold(i)%EJvQHP(nqhp))
		endif
	enddo
	
	call UnMakeGaps()

	MassTot=0d0
	do j=1,D%nTheta-1
	call tellertje(j,D%nTheta-1)
	Cold(0:D%nR)=C(0:D%nR,j)
	Cold(D%nR)%gasfrac=C(D%nR-1,j)%gasfrac
	Cold(D%nR)%dens=C(D%nR-1,j)%dens
	Cold(D%nR)%dens0=C(D%nR-1,j)%dens0
	Cold(D%nR)%gasdens=C(D%nR-1,j)%gasdens
	Cold(D%nR)%w(1:ngrains)=C(D%nR-1,j)%w(1:ngrains)
	Cold(D%nR)%w0(1:ngrains)=C(D%nR-1,j)%w0(1:ngrains)
	if(use_qhp) Cold(D%nR)%Tqhp=C(D%nR-1,j)%Tqhp
	if(use_qhp) Cold(D%nR)%LRF=C(D%nR-1,j)%LRF
	if(use_qhp) Cold(D%nR)%tdistr=C(D%nR-1,j)%tdistr
	if(use_qhp) Cold(D%nR)%EJvQHP=C(D%nR-1,j)%EJvQHP
	if(.not.tcontact.or.tdes_iter) then
		Cold(D%nR)%EJvP(1:ngrains)=C(D%nR-1,j)%EJvP(1:ngrains)
	endif
	TdiskStruct(0:D%nR)=Temp(0:D%nR,j)
	TdiskStruct(D%nR)=Temp(D%nR-1,j)
	do i=1,D%nR-1
		C(i,j)%xedge(1)=Rnew(i)
		C(i,j)%xedge(2)=Rnew(i+1)
		C(i,j)%xedge2(1)=Rnew(i)**2
		C(i,j)%xedge2(2)=Rnew(i+1)**2
		i1=1
		i2=D%nR-1
		do k=1,D%nR-1
			if(Rnew(i).ge.D%R(k).and.Rnew(i).lt.D%R(k+1)) then
				i1=k
				goto 4
			endif
		enddo
4		continue
		do k=1,D%nR-1
			if(Rnew(i+1).ge.D%R(k).and.Rnew(i+1).lt.D%R(k+1)) then
				i2=k
				goto 5
			endif
		enddo
5		continue


		C(i,j)%V=(4d0*pi/3d0)*(Rnew(i+1)**3-Rnew(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3

		Temp(i,j)=0d0
		if(i1.ne.i2) then
			V=(4d0*pi/3d0)*(D%R(i1+1)**3-Rnew(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
			R=10d0**(log10(Rnew(i)*D%R(i1+1))/2d0)
			dens=log10(Cold(i1)%gasdens)+
     &				log10(Cold(i1+1)%gasdens/Cold(i1)%gasdens)*
     &				log10(R*AU/D%R_av(i1))/log10(D%R_av(i1+1)/D%R_av(i1))
			Mg=V*10d0**dens
			dens=log10(Cold(i1)%dens0)+
     &				log10(Cold(i1+1)%dens0/Cold(i1)%dens0)*
     &				log10(R*AU/D%R_av(i1))/log10(D%R_av(i1+1)/D%R_av(i1))
			M0=V*10d0**dens
			dens=log10(Cold(i1)%dens)+
     &				log10(Cold(i1+1)%dens/Cold(i1)%dens)*
     &				log10(R*AU/D%R_av(i1))/log10(D%R_av(i1+1)/D%R_av(i1))
			M=V*10d0**dens
			Temp(i,j)=Temp(i,j)+TdiskStruct(i1)**4*M0
			Ni=real(Cold(i1)%Ni)*M
			C(i,j)%EJv=Cold(i1)%EJv*M
			if(use_qhp) C(i,j)%EJvQHP=Cold(i1)%EJvQHP*M
			if(use_qhp) C(i,j)%LRF(1:nlam)=Cold(i1)%LRF(1:nlam)*M
			if(use_qhp) C(i,j)%tdistr(1:nqhp,1:NTQHP)=Cold(i1)%tdistr(1:nqhp,1:NTQHP)*M
			C(i,j)%EJv2=Cold(i1)%EJv2*M**2
			if(use_qhp) C(i,j)%Tqhp(1:nqhp)=Cold(i1)%Tqhp(1:nqhp)**4*M0
			do ii=1,ngrains
				C(i,j)%w(ii)=Cold(i1)%w(ii)*M
				C(i,j)%w0(ii)=Cold(i1)%w0(ii)*M
				if(.not.tcontact.or.tdes_iter) then
					C(i,j)%EJvP(ii)=Cold(i1)%EJvP(ii)*M
				endif
			enddo
			C(i,j)%mass=M
			totM0=M0
			V=(4d0*pi/3d0)*(Rnew(i+1)**3-D%R(i2)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
			R=10d0**(log10(Rnew(i+1)*D%R(i2))/2d0)
			dens=log10(Cold(i2-1)%gasdens)+
     &				log10(Cold(i2)%gasdens/Cold(i2-1)%gasdens)*
     &				log10(R*AU/D%R_av(i2-1))/log10(D%R_av(i2)/D%R_av(i2-1))
			Mg=Mg+V*10d0**dens
			dens=log10(Cold(i2-1)%dens0)+
     &				log10(Cold(i2)%dens0/Cold(i2-1)%dens0)*
     &				log10(R*AU/D%R_av(i2-1))/log10(D%R_av(i2)/D%R_av(i2-1))
			M0=V*10d0**dens
			dens=log10(Cold(i2-1)%dens)+
     &				log10(Cold(i2)%dens/Cold(i2-1)%dens)*
     &				log10(R*AU/D%R_av(i2-1))/log10(D%R_av(i2)/D%R_av(i2-1))
			M=V*10d0**dens
			Temp(i,j)=Temp(i,j)+TdiskStruct(i2)**4*M0
			Ni=Ni+real(Cold(i2)%Ni)*M
			C(i,j)%EJv=C(i,j)%EJv+Cold(i2)%EJv*M
			if(use_qhp) C(i,j)%EJvQHP=C(i,j)%EJvQHP+Cold(i2)%EJvQHP*M
			if(use_qhp) C(i,j)%LRF(1:nlam)=C(i,j)%LRF(1:nlam)+Cold(i2)%LRF(1:nlam)*M
			if(use_qhp) C(i,j)%tdistr(1:nqhp,1:NTQHP)=C(i,j)%tdistr(1:nqhp,1:NTQHP)
     &					+Cold(i2)%tdistr(1:nqhp,1:NTQHP)*M
			C(i,j)%EJv2=C(i,j)%EJv2+Cold(i2)%EJv2*M**2
			if(use_qhp) C(i,j)%Tqhp(1:nqhp)=C(i,j)%Tqhp(1:nqhp)+Cold(i2)%Tqhp(1:nqhp)**4*M0
			do ii=1,ngrains
				C(i,j)%w(ii)=C(i,j)%w(ii)+Cold(i2)%w(ii)*M
				C(i,j)%w0(ii)=C(i,j)%w0(ii)+Cold(i2)%w0(ii)*M
				if(.not.tcontact.or.tdes_iter) then
					C(i,j)%EJvP(ii)=C(i,j)%EJvP(ii)+Cold(i2)%EJvP(ii)*M
				endif
			enddo
			C(i,j)%mass=C(i,j)%mass+M
			totM0=totM0+M0
			do k=i1+1,i2-1
				Mg=Mg+Cold(k)%V*Cold(k)%gasdens
				M=Cold(k)%mass
				M0=Cold(k)%V*Cold(k)%dens0
				Temp(i,j)=Temp(i,j)+TdiskStruct(k)**4*M0
				Ni=Ni+real(Cold(k)%Ni)*M
				C(i,j)%EJv=C(i,j)%EJv+Cold(k)%EJv*M
				if(use_qhp) C(i,j)%EJvQHP=C(i,j)%EJvQHP+Cold(k)%EJvQHP*M
				if(use_qhp) C(i,j)%LRF(1:nlam)=C(i,j)%LRF(1:nlam)+Cold(k)%LRF(1:nlam)*M
				if(use_qhp) C(i,j)%tdistr(1:nqhp,1:NTQHP)=C(i,j)%tdistr(1:nqhp,1:NTQHP)
     &						+Cold(k)%tdistr(1:nqhp,1:NTQHP)*M
				C(i,j)%EJv2=C(i,j)%EJv2+Cold(k)%EJv2*M**2
				if(use_qhp) C(i,j)%Tqhp(1:nqhp)=C(i,j)%Tqhp(1:nqhp)+Cold(k)%Tqhp(1:nqhp)**4*M0
				do ii=1,ngrains
					C(i,j)%w(ii)=C(i,j)%w(ii)+Cold(k)%w(ii)*M
					C(i,j)%w0(ii)=C(i,j)%w0(ii)+Cold(k)%w0(ii)*M
					if(.not.tcontact.or.tdes_iter) then
						C(i,j)%EJvP(ii)=C(i,j)%EJvP(ii)+Cold(k)%EJvP(ii)*M
					endif
				enddo
				C(i,j)%mass=C(i,j)%mass+M
				totM0=totM0+M0
			enddo
			M0=totM0
			C(i,j)%gasdens=Mg/C(i,j)%V
			C(i,j)%dens=C(i,j)%mass/C(i,j)%V
			C(i,j)%dens0=M0/C(i,j)%V
			Temp(i,j)=(Temp(i,j)/M0)**0.25d0
			C(i,j)%gasfrac=1d0-C(i,j)%dens/C(i,j)%dens0
			C(i,j)%Ni=Ni/C(i,j)%mass
			C(i,j)%mass=C(i,j)%dens*C(i,j)%V
			C(i,j)%EJv=C(i,j)%EJv/C(i,j)%mass
			if(use_qhp) C(i,j)%EJvQHP=C(i,j)%EJvQHP/C(i,j)%mass
			if(use_qhp) C(i,j)%LRF(1:nlam)=C(i,j)%LRF(1:nlam)/C(i,j)%mass
			if(use_qhp) C(i,j)%tdistr(1:nqhp,1:NTQHP)=
     &					C(i,j)%tdistr(1:nqhp,1:NTQHP)/C(i,j)%mass
			C(i,j)%EJv2=C(i,j)%EJv2/C(i,j)%mass**2
			if(use_qhp) C(i,j)%Tqhp(1:nqhp)=(C(i,j)%Tqhp(1:nqhp)/M0)**0.25d0
			do ii=1,ngrains
				C(i,j)%w(ii)=C(i,j)%w(ii)/C(i,j)%mass
				C(i,j)%w0(ii)=C(i,j)%w0(ii)/C(i,j)%mass
				if(.not.tcontact.or.tdes_iter) then
					C(i,j)%EJvP(ii)=C(i,j)%EJvP(ii)/C(i,j)%mass
				endif
			enddo
			call CheckMinimumDensity(i,j)
		else
			if(Rav(i).gt.D%R_av(i1).or.i1.eq.1) then
				dens=log10(Cold(i1)%gasdens)+
     &					log10(Cold(i1+1)%gasdens/Cold(i1)%gasdens)*
     &					log10(Rav(i)/D%R_av(i1))/log10(D%R_av(i1+1)/D%R_av(i1))
				C(i,j)%gasdens=10d0**dens
				dens=log10(Cold(i1)%dens0)+
     &					log10(Cold(i1+1)%dens0/Cold(i1)%dens0)*
     &					log10(Rav(i)/D%R_av(i1))/log10(D%R_av(i1+1)/D%R_av(i1))
				C(i,j)%dens0=10d0**dens
				dens=log10(Cold(i1)%dens)+
     &					log10(Cold(i1+1)%dens/Cold(i1)%dens)*
     &					log10(Rav(i)/D%R_av(i1))/log10(D%R_av(i1+1)/D%R_av(i1))
				C(i,j)%dens=10d0**dens
			else
				dens=log10(Cold(i1-1)%gasdens)+
     &					log10(Cold(i1)%gasdens/Cold(i1-1)%gasdens)*
     &					log10(Rav(i)/D%R_av(i1-1))/log10(D%R_av(i1)/D%R_av(i1-1))
				C(i,j)%gasdens=10d0**dens
				dens=log10(Cold(i1-1)%dens0)+
     &					log10(Cold(i1)%dens0/Cold(i1-1)%dens0)*
     &					log10(Rav(i)/D%R_av(i1-1))/log10(D%R_av(i1)/D%R_av(i1-1))
				C(i,j)%dens0=10d0**dens
				dens=log10(Cold(i1-1)%dens)+
     &					log10(Cold(i1)%dens/Cold(i1-1)%dens)*
     &					log10(Rav(i)/D%R_av(i1-1))/log10(D%R_av(i1)/D%R_av(i1-1))
				C(i,j)%dens=10d0**dens
			endif
			C(i,j)%mass=C(i,j)%dens*C(i,j)%V
			C(i,j)%gasfrac=Cold(i1)%gasfrac
			Temp(i,j)=TdiskStruct(i1)
			C(i,j)%EJv=Cold(i1)%EJv
			if(use_qhp) C(i,j)%EJvQHP=Cold(i1)%EJvQHP
			if(use_qhp) C(i,j)%LRF(1:nlam)=Cold(i1)%LRF(1:nlam)
			if(use_qhp) C(i,j)%tdistr(1:nqhp,1:NTQHP)=Cold(i1)%tdistr(1:nqhp,1:NTQHP)
			C(i,j)%EJv2=Cold(i1)%EJv2
			if(use_qhp) C(i,j)%Tqhp(1:nqhp)=Cold(i1)%Tqhp(1:nqhp)
			C(i,j)%Ni=Cold(i1)%Ni*C(i,j)%mass/Cold(i1)%mass
			do ii=1,ngrains
				C(i,j)%w(ii)=Cold(i1)%w(ii)
				C(i,j)%w0(ii)=Cold(i1)%w0(ii)
				if(.not.tcontact.or.tdes_iter) then
					C(i,j)%EJvP(ii)=Cold(i1)%EJvP(ii)
				endif
			enddo
		endif
		call CheckMinimumDensity(i,j)
		tot=0d0
		do ii=1,ngrains
			tot=tot+C(i,j)%w(ii)
		enddo
		if(tot.eq.0d0) C(i,j)%w(1:ngrains)=C(i,j)%w0(1:ngrains)
		phot%i=i
		phot%j=j
		phot%E=C(i,j)%EJv
		C(i,j)%T=determineT(phot)
		if(C(i,j)%Ni.ne.0) then
			C(i,j)%dEJv=(C(i,j)%EJv
     &				+sqrt(abs(C(i,j)%EJv2-C(i,j)%EJv**2/real(C(i,j)%Ni))))
     &				/sqrt(real(C(i,j)%Ni))
			C(i,j)%dT=0.25d0*C(i,j)%T*C(i,j)%dEJv/C(i,j)%EJv
		else
			C(i,j)%dT=real(TMAX)/dT
			C(i,j)%dEJv=0d0
			do ii=1,ngrains
			   do iopac=1,Grain(ii)%nopac
			      C(i,j)%dEJv=C(i,j)%dEJv+C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)*Grain(ii)%Kp(iopac,TMAX)
			   enddo
			enddo
		endif
		if(.not.tcontact.or.tdes_iter) then
		do ii=1,ngrains
			phot%i=i
			phot%j=j
c			C(i,j)%EJvP(ii)=C(i,j)%EJvP(ii)/C(i,j)%mass
			phot%E=C(i,j)%EJvP(ii)
			C(i,j)%TP(ii)=determineTP(phot,ii)
		enddo
		endif
		MassTot=MassTot+C(i,j)%mass
	enddo
	enddo
	
	do ii=1,ngrains
	if(Grain(ii)%settle) then
		do i=1,D%nR-1
			if(Rav(i).le.D%R_av(1)) then
				j=1
				goto 50
			endif
			do j=1,D%nR-2
				if(Rav(i).ge.D%R_av(j).and.Rav(i).le.D%R_av(j+1)) goto 50
			enddo
			j=D%nR-2
50			shtemp(i)=Grain(ii)%shscale(j)+
     &	(Grain(ii)%shscale(j+1)-Grain(ii)%shscale(j))*(Rav(i)-D%R_av(j))/(D%R_av(j+1)-D%R_av(j))
		enddo
		Grain(ii)%shscale(1:D%nR-1)=shtemp(1:D%nR-1)
	endif
	enddo
	do i=1,D%nR-1
		if(Rav(i).le.D%R_av(1)) then
			j=1
			goto 60
		endif
		do j=1,D%nR-2
			if(Rav(i).ge.D%R_av(j).and.Rav(i).le.D%R_av(j+1)) goto 60
		enddo
		j=D%nR-2
60		shtemp(i)=shscale(j)+(shscale(j+1)-shscale(j))*(Rav(i)-D%R_av(j))/(D%R_av(j+1)-D%R_av(j))
	enddo
	shscale(1:D%nR-1)=shtemp(1:D%nR-1)
	
	do i=1,D%nR
		D%R(i)=Rnew(i)
	enddo
	do i=0,D%nR-1
		D%R_av(i)=Rav(i)
	enddo

	do i=0,D%nR
		deallocate(Cold(i)%w)
		deallocate(Cold(i)%w0)
		if(.not.tcontact.or.tdes_iter) then
			deallocate(Cold(i)%TP)
			deallocate(Cold(i)%EJvP)
		endif
		if(use_qhp) then
			deallocate(Cold(i)%LRF)
			deallocate(Cold(i)%QHP)
			deallocate(Cold(i)%tdistr)
			deallocate(Cold(i)%Tqhp)
			deallocate(Cold(i)%EJvQHP)
		endif
	enddo

	call MakeGaps()

	MassTot=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		MassTot=MassTot+C(i,j)%dens0*C(i,j)%V
	enddo
	enddo

	if(dabs(1d0-MassTot0/MassTot).gt.1d-2) then
		write(*,'("Large error in regridding: ",f10.2,"%")') 100d0*(1d0-MassTot0/MassTot)
		write(9,'("Large error in regridding: ",f10.2,"%")') 100d0*(1d0-MassTot0/MassTot)
	endif
	do i=1,D%nR-1
	do j=1,D%nTheta
		C(i,j)%mass=C(i,j)%mass*MassTot0/MassTot
		C(i,j)%dens=C(i,j)%dens*MassTot0/MassTot
		C(i,j)%dens0=C(i,j)%dens0*MassTot0/MassTot
		C(i,j)%gasdens=C(i,j)%gasdens*MassTot0/MassTot
		C(i,j)%Tav=Temp(i,j)
	enddo
	enddo

	write(surfdensfile,'(a,"surfacedens.dat")') outdir(1:len_trim(outdir))
	open(unit=90,file=surfdensfile,RECL=1000)
	do i=1,D%nR-1
	MassTot=0d0
	MassTot0=0d0
	do j=1,D%nTheta
		MassTot=MassTot+C(i,j)%mass
		MassTot0=MassTot0+C(i,j)%dens0*C(i,j)%V
	enddo
	write(90,*) D%R_av(i)/AU,MassTot/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2),MassTot0/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)
	enddo
	close(unit=90)

	return
	end
	

	real*8 function powinterpol(x,x1,x2,y1,y2)
	IMPLICIT NONE
	real*8 x1,x2,y1,y2,x,pow
	pow=log10(y2/y1)/log10(x2/x1)
	powinterpol=y1*(x/x1)**pow
	return
	end


