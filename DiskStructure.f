c History:	20071003 MK: implemented the improved dens-Temp fit suggested by Carsten&Michiel
c		20071008 MK: if Monte Carlo temperature less than 1K, gas fraction set to zero.
c		20071009 MM: Improved interpolation for the regridding of the density
c		20071011 MM: Regrid using the optical depth. 
c					 Finest grid is from tau=0.1-10
c		20071015 MM: Fixed the smoothing of the T-structure before
c					 determining the scaleheight. There was a bug
c					 causing it to be off.
c		20071022 MM: Turned the regridding with the optical depth
c					 on for all runs (also the ones without destruction)
c		20080410 MM: Fixed an issue with the w0 and w in the DiskStructure
c					 subroutine. This caused strange effects when using
c					 thermal contact off and dust destruction on simultaneously.

	Module Diskstruct
		use Parameters
		integer i,j,j1,j2
		integer,allocatable :: iRfirst(:,:)
		real*8,allocatable :: z(:),rho(:,:),Temp(:,:),Temp2(:)
		real*8 scale
		real*8,allocatable :: fBW(:,:),Rthindes(:),Rbwdes(:)
	end module Diskstruct


	subroutine TempAverage(fw)
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	integer ii
	real*8 T,tot,fw
	
c ------ for the IN05 interior T -------------
c	call InteriorIN05
c ------ for the IN05 interior T -------------

	if(.not.allocated(Temp)) then
		allocate(Temp(0:D%nR,0:D%nTheta))
		if(useTgas.or.viscous) then
			Temp(0:D%nR,0:D%nTheta)=C(0:D%nR,0:D%nTheta)%Tgas
		else
			Temp(0:D%nR,0:D%nTheta)=C(0:D%nR,0:D%nTheta)%T
		endif
		return
	endif

	do i=0,D%nR-1
		do j=1,D%nTheta-1
			T=0d0
			tot=0d0
			if(useTgas.or.viscous) then
				T=C(i,j)%Tgas
			else
c				do ii=1,ngrains
c					if(.not.Grain(ii)%qhp) then
c						if(tcontact.or.C(i,j)%diff) then
c							T=T+C(i,j)%w(ii)*C(i,j)%T
c							tot=tot+C(i,j)%w(ii)
c						else
c							T=T+C(i,j)%w(ii)*C(i,j)%TP(ii)
c							tot=tot+C(i,j)%w(ii)
c						endif
c					endif	
c				enddo
c				if(tot.gt.1d-6) then
c					T=T/tot
c				else
					T=C(i,j)%T
c				endif
			endif
			if(T.lt.2.7) T=2.7
			Temp(i,j)=(1d0-fw)*Temp(i,j)+fw*T
			C(i,j)%Tav=Temp(i,j)
		enddo
	enddo

	return
	end

c-----------------------------------------------------------------------
c This subroutine determines the vertical structure of the disk
c using the temperature and density structure.
c Note that the disk will collapse if at any point the temperature
c is set to 0.
c-----------------------------------------------------------------------
	subroutine DiskStructure()
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	real*8 M0(0:ngrains),M1(0:ngrains),tot,gasf(ngrains),T,M(0:ngrains)
	real*8 ystart(2),eps,h1,hmin,gasdev,rhosig(D%nTheta),Rcyl,x,y
	doubleprecision errtot
	external derivs
	integer ii,fix,iz
	real*8 rr,zz,hr,f1,f2
	
	real*8 totG,xxx,yyy,zzz,xx0,yy0,zz0,sint
	integer iii,jjj,kkk
	
	write(*,'("Solving vertical structure")')
	write(9,'("Solving vertical structure")')

	allocate(z(D%nTheta+1))
	allocate(Temp2(D%nTheta+1))
	allocate(rho(D%nTheta+1,0:ngrains))

	call UnMakeGaps()

	if (fixmpset) call FixedMPSettling()
	if (mpset) call MPSettling() !Gijsexp

	eps=1d-6
	do i=1,D%nR-1

c	do j=1,D%nTheta-1
c		xx0=D%R_av(i)*sin(D%theta_av(j))
c		yy0=0d0
c		zz0=D%R_av(i)*cos(D%theta_av(j))
c		totG=0d0
c		do iii=1,D%nR-1
c		do jjj=1,D%nTheta-1
c			sint=D%R_av(iii)*sin(D%theta_av(jjj))
c			zzz=D%R_av(iii)*cos(D%theta_av(jjj))
c			do kkk=1,4
c				xxx=cos(pi*(real(kkk)-0.5d0)/9d0)*sint
c				yyy=sin(pi*(real(kkk)-0.5d0)/9d0)*sint
c				totG=totG+(C(iii,jjj)%mass/8d0)*(zz0-zzz)/(sqrt((xxx-xx0)**2+(yyy-yy0)**2+(zzz-zz0)**2)**3)
c				zzz=-zzz
c				totG=totG+(C(iii,jjj)%mass/8d0)*(zz0-zzz)/(sqrt((xxx-xx0)**2+(yyy-yy0)**2+(zzz-zz0)**2)**3)
c			enddo
c		enddo
c		enddo
c		C(i,j)%FradZ=totG*gas2dust/(D%Mstar*zz0/(sqrt(xx0**2+yy0**2+zz0**2)**3))
c	enddo

		call tellertje(i,D%nR-1)
		do ii=0,ngrains
		rhosig(1:D%nTheta-1)=0d0
		do j=1,D%nTheta-1
			Temp2(j)=Temp(i,j)
		enddo

		! Use T=Tmidplane for vertical structure
		if (mpstr) then
		   do j=1,D%nTheta-1
		      Temp2(j)=Temp(i,D%nTheta-1)
		   enddo
		endif

		! Scale factor for grain ii
		if(ii.ne.0) then
			scale=shscale(i)
			if(Grain(ii)%settle) scale=scale*Grain(ii)%shscale(i)
			if(Grain(ii)%gascoupled.and.scale.lt.1d0) scale=1d0
		else
			scale=1d0*shscale(i)
		endif

		! Total dust mass for grain ii in cells i
		M0(ii)=0d0
		M(ii)=0d0
		do j=1,D%nTheta-1
			if(ii.ne.0) then
				M0(ii)=M0(ii)+C(i,j)%dens0*C(i,j)%V*C(i,j)%w0(ii)
				M(ii)=M(ii)+C(i,j)%dens*C(i,j)%V*C(i,j)%w(ii)
			else
				M0(ii)=M0(ii)+C(i,j)%gasdens*C(i,j)%V
			endif
			z(j)=D%R_av(i)*cos(D%theta_av(j))
		enddo

		! Use zones?
		if(ii.ne.0) then
			fix=0
			do iz=1,nzones
				if(D%R_av(i).gt.(Zone(iz)%Rin*AU).and.D%R_av(i).lt.(Zone(iz)%Rout*AU)) then
					if(Zone(iz)%fix_struct.and.Zone(iz)%inc_grain(ii)) then
						fix=iz
					endif
				endif
			enddo
		else
			fix=0
			do iz=1,nzones
				if(D%R_av(i).gt.(Zone(iz)%Rin*AU).and.D%R_av(i).lt.(Zone(iz)%Rout*AU)) then
					if(Zone(iz)%fix_struct) then
						fix=iz
					endif
				endif
			enddo
			do iz=1,nzones
				if(D%R_av(i).gt.(Zone(iz)%Rin*AU).and.D%R_av(i).lt.(Zone(iz)%Rout*AU)) then
					if(.not.Zone(iz)%fix_struct) then
						fix=0
					endif
				endif
			enddo
		endif
		if(fix.eq.0) then
c	go ahead, no zones with fixed structure

		! Solve the vertical structure of the gas 
		if(ii.eq.0) then
			if(scale.gt.0d0) then
				rho(D%nTheta-1,ii)=0d0
				do j=D%nTheta-2,1,-1
					j1=j+1
					j2=j
					ystart(1)=rho(j1,ii)
					h1=abs(z(j2)-z(j1))*1d-6
					hmin=0d0
					call odeint(ystart,1,z(j1),z(j2),eps,h1,hmin,derivs)
					rho(j2,ii)=ystart(1)
				enddo
			else ! wedge
				rho(D%nTheta-1,ii)=0d0
				do j=D%nTheta-2,1,-1
					if((1d0/tan(D%theta_av(j))).lt.(-scale)) then
						rho(j,ii)=0d0
					else
						rho(j,ii)=-60d0
					endif
				enddo
			endif

		! Solve the vertical structure of the dust
		else if(Grain(ii)%shtype.eq.'DISK') then
			if(scale.gt.0d0) then
				rho(D%nTheta-1,ii)=0d0
				do j=D%nTheta-2,1,-1
					j1=j+1
					j2=j
					ystart(1)=rho(j1,ii)
					h1=abs(z(j2)-z(j1))*1d-6
					hmin=0d0
					call odeint(ystart,1,z(j1),z(j2),eps,h1,hmin,derivs)
					rho(j2,ii)=ystart(1)
				enddo
			else ! wedge
				rho(D%nTheta-1,ii)=0d0
				do j=D%nTheta-2,1,-1
					if((1d0/tan(D%theta_av(j))).lt.(-scale)) then
						rho(j,ii)=0d0
					else
						rho(j,ii)=-60d0
					endif
				enddo
			endif
		else if(Grain(ii)%shtype.eq.'ZODI') then
			do j=1,D%nTheta-1
c models for zodiacal dust from Giese et al. (1986)
c cosine model from Rittich (1986)
				rho(j,ii)=log(0.15d0+0.85*(cos(pi/2d0-D%theta_av(j)))**28d0)
c model from Leinert et al. 2002
c				rho(j,ii)=-4.8d0*abs(sin(pi/2d0-D%theta_av(j)))**1.3d0
			enddo
		else if(Grain(ii)%shtype.eq.'HALO') then
		   do j=1,D%nTheta-1
		      rho(j,ii)=0d0

		      !  a wedge with constant density for z/r < -settle#
		      if(Grain(ii)%shscale(i).lt.0d0.and.
     1                     pi/2d0-D%theta_av(j).gt.-(Grain(ii)%shscale(i))) then
			 rho(j,ii)=-60d0
		      endif
		   enddo
		else if(Grain(ii)%shtype.eq.'LOBE') then
			do j=1,D%nTheta-1
				rho(j,ii)=log10(exp(-((D%theta_av(j)*180d0/pi)/scale)**2))
			enddo
		else if(Grain(ii)%shtype.eq.'JET') then
			do j=1,D%nTheta-1
				rho(j,ii)=-((D%thet(j)*180d0/pi)/(scale*(D%R_av(i)/D%R_av(D%nR-1))))**2
c				if(D%theta_av(j)*180d0/pi.lt.scale*(D%R_av(i)/D%R_av(D%nR-1)).or.j.eq.1) then
c					rho(j,ii)=0d0	!(D%theta_av(j)*180d0/pi)**2-scale**2
c				else
c					rho(j,ii)=-60d0
c				endif
			enddo
		else if(Grain(ii)%shtype.eq.'SHELL') then
			do j=1,D%nTheta-1
				rho(j,ii)=log(1d0+Grain(ii)%shscale(i)*(sin(D%theta_av(j)))**3d0)
c				if(Grain(ii)%shscale(i).gt.0d0) then
c					rho(j,ii)=log((cos(D%theta_av(j)))**Grain(ii)%shscale(i))
c				else
c					rho(j,ii)=log((sin(D%theta_av(j)))**(-Grain(ii)%shscale(i)))
c				endif
			enddo
		else if(Grain(ii)%shtype.eq.'CYLI') then
			rho(D%nTheta-1,ii)=0d0
			do j=D%nTheta-2,1,-1
				if((z(j)/AU).lt.(scale)) then
					rho(j,ii)=0d0
				else
					rho(j,ii)=-60d0
				endif
			enddo
		else if(Grain(ii)%shtype.eq.'CAVITY') then
			rho(D%nTheta-1,ii)=0d0
			rho(1,ii)=-60d0
			do j=D%nTheta-2,2,-1
				Rcyl=D%R_av(i)*sin(D%theta_av(j))
				if((z(j)/AU).lt.((Rcyl/AU)**2*scale)) then
					rho(j,ii)=0d0
				else
					rho(j,ii)=-60d0
				endif
			enddo
		else if(Grain(ii)%shtype.eq.'WEDGE') then
			if(scale.gt.0d0) then
				do j=1,D%nTheta-1
					if((180d0*D%theta_av(j)/pi).gt.scale) then
						rho(j,ii)=0d0
					else
						rho(j,ii)=-60d0
					endif
				enddo
			else
				do j=1,D%nTheta-1
					if((180d0*D%theta_av(j)/pi).lt.(-scale)) then
						rho(j,ii)=0d0
					else
						rho(j,ii)=-60d0
					endif
				enddo
			endif				
c			do j=D%nTheta-1,1,-1
c				rho(j,ii)=log(1d0+Grain(ii)%shscale(i)*exp(-((90d0-180d0*D%theta_av(j)/pi)/shscale(i))**2d0))
c			enddo
		else
			if(i.eq.1) then
c 			GFORTRAN changed i to i5 otherwise get \
c                Nonnegative width required in format string at with gfortran
				write(*,'("Scale height type unknown!!",i2)') ii
				write(9,'("Scale height type unknown!!",i2)') ii
				write(*,'("Leaving it as it is")')
				write(9,'("Leaving it as it is")')
			endif
			do j=1,D%nTheta-1
				if(C(i,j)%dens*C(i,j)%w(ii).gt.1d-60) then
					rho(j,ii)=log(C(i,j)%dens*C(i,j)%w(ii))
				else
					rho(j,ii)=-60d0
				endif
			enddo
c			stop
		endif ! end setting dust scale height

		! continue with zones
		else
c	fix scaleheight to the one defined in Zone(fix)
			if(Zone(fix)%sh.lt.0d0) then
c	this is a wedge zone!
				rho(D%nTheta-1,ii)=0d0
				do j=1,D%nTheta-2
					if((180d0*D%theta_av(j)/pi).gt.(-Zone(fix)%sh)) then
						rho(j,ii)=0d0
					else
						rho(j,ii)=-60d0
					endif
				enddo
			else
c	this is not a wedge zone!
				do j=1,D%nTheta-1
					rr=D%R_av(i)*sin(D%theta_av(j))/AU
					zz=D%R_av(i)*cos(D%theta_av(j))/AU
					hr=scale*Zone(fix)%sh*(rr/Zone(fix)%Rsh)**Zone(fix)%shpow
					f1=rr**(-Zone(fix)%denspow)
c					f2=exp(-(zz/hr)**2)
c					if(j.eq.D%nTheta-1) f2=1d0
					f2=-(zz/hr)**2
					if(j.eq.D%nTheta-1) f2=0d0
c					rho(j,ii)=log(f1*f2/hr)
					rho(j,ii)=log(f1/hr)+f2
				enddo
			endif
		endif
		
		! Recalculate mass of grain ii in cells i
		M1(ii)=0d0
		do j=1,D%nTheta-1
			M1(ii)=M1(ii)+exp(rho(j,ii))*C(i,j)%V
		enddo

		! Renormalize mass of grain ii in cells i
		do j=1,D%nTheta-1
			rho(j,ii)=exp(rho(j,ii))*M0(ii)/M1(ii)
c			if(IsNaN(rho(j,ii))) rho(j,ii)=1d-50
			if(rho(j,ii).lt.1d-50) rho(j,ii)=1d-50
		enddo

		enddo ! ii

		! Settling in a vertical slab, for all grains
		if (scset) then
		   call settling(rho(1:D%nTheta-1,0:ngrains),
     &		                 D%nTheta-1,ngrains,i)
		endif

		! Set density, abundance etc.
		do j=1,D%nTheta-1
			do ii=1,ngrains
				if(C(i,j)%w0(ii).ge.1d-100) then
					gasf(ii)=1d0-C(i,j)%dens*C(i,j)%w(ii)/(C(i,j)%dens0*C(i,j)%w0(ii))
				else
					gasf(ii)=0d0
				endif
				if(gasf(ii).lt.0d0) gasf(ii)=0d0
				if(gasf(ii).gt.1d0) gasf(ii)=1d0
				if(IsNaN(gasf(ii))) gasf(ii)=0d0
			enddo
			C(i,j)%gasdens=rho(j,0)
			C(i,j)%dens0=0d0
			do ii=1,ngrains
				C(i,j)%dens0=C(i,j)%dens0+rho(j,ii)
			enddo
c			if(.not.tdes_iter) then
			C(i,j)%dens=0d0
			do ii=1,ngrains
				if(M0(ii).ne.0d0) C(i,j)%dens=C(i,j)%dens+rho(j,ii)*M(ii)/M0(ii)
			enddo
c			endif
			if(C(i,j)%dens0.gt.1d-50) then
				do ii=1,ngrains
					C(i,j)%w0(ii)=rho(j,ii)/C(i,j)%dens0
				enddo
			endif
			tot=0d0
			do ii=1,ngrains
				tot=tot+C(i,j)%w0(ii)
			enddo
			do ii=1,ngrains
				C(i,j)%w0(ii)=C(i,j)%w0(ii)/tot
				C(i,j)%w(ii)=C(i,j)%w0(ii)*(1d0-gasf(ii))
			enddo
			tot=0d0
			do ii=1,ngrains
				tot=tot+C(i,j)%w(ii)
			enddo
			if(tot.gt.1d-50) then
				do ii=1,ngrains
					C(i,j)%w(ii)=C(i,j)%w(ii)/tot
				enddo
			else
				do ii=1,ngrains
					C(i,j)%w(ii)=C(i,j)%w0(ii)
				enddo
			endif
			C(i,j)%dens=C(i,j)%dens0*tot
			call CheckMinimumDensity(i,j)
		enddo ! j
	enddo ! i

	if(haloswitch) call MakeHalo()

	call MakeGaps()

	deallocate(z)
	deallocate(Temp2)
	deallocate(rho)

	call CheckCells()

	if(raditer.or.getalpha) call RadialStruct()

	return
	end


	subroutine derivs(x,y,dydx)
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	real*8 x,y,dydx
	real*8 mu,G
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
c	parameter(mu=1.4*1.67262158d-24) !1.4 times the proton mass in gram
	parameter(G=6.67300d-8) ! in cm^3/g/s
	real*8 T,dlnT,beta
	integer ii
	
	beta=C(i,j1)%FradZ+(C(i,j2)%FradZ-C(i,j1)%FradZ)*(x-z(j1))/(z(j2)-z(j1))

c	T=Temp(i,j1)+(Temp(i,j2)-Temp(i,j1))*(x-z(j1))/(z(j2)-z(j1))
c	dlnT=(log(Temp(i,j1))-log(Temp(i,j2)))/(z(j1)-z(j2))
c	dydx=-mu*G*D%Mstar*x/(scale**2*kb*T*D%R_av(i)**3)-dlnT

	T=Temp2(j1)+(Temp2(j2)-Temp2(j1))*(x-z(j1))/(z(j2)-z(j1))
	dlnT=(log(Temp2(j1))-log(Temp2(j2)))/(z(j1)-z(j2))
c	dydx=-mu*G*D%Mstar*x/(scale**2*kb*T*D%R_av(i)**3)-dlnT
	dydx=-mu*G*D%Mstar*(1d0+beta)*x/(scale**2*kb*T*D%R_av(i)**3)-dlnT

	return
	end

	subroutine MakeHalo()
	use Parameters
	IMPLICIT NONE
	real*8 mindens,mass,massgrains,Agrains
	real*8 masstot,masstot0
	integer i,j,ii,jj,iopac
	
	
	do i=1,D%nR-1
	mindens=1d-50
	masstot0=0d0
	do j=1,D%nTheta-1
		masstot0=masstot0+C(i,j)%dens*C(i,j)%V
	enddo
	do j=D%nTheta-2,1,-1
		mass=(C(i,j)%gasdens*gas2dust+C(i,j)%dens)*D%R_av(i)*pi
		massgrains=0d0
		Agrains=0d0
		do ii=1,ngrains
			massgrains=massgrains+(4d0*pi*Grain(ii)%rv**3*Grain(ii)%rho(1)/3d0)*C(i,j)%w(ii)
c			Agrains=Agrains+Grain(ii)%rv**2*pi*C(i,j)%w(ii)
			do iopac=1,Grain(ii)%nopac
			   Agrains=Agrains+Grain(ii)%Kpstar(iopac)*(4d0*pi*Grain(ii)%rv**3*Grain(ii)%rho(1)/3d0)*C(i,j)%w(ii)
			enddo
		enddo
		mass=0d0
		do jj=j,1,-1
			mass=mass+(C(i,jj)%gasdens*gas2dust+C(i,jj)%dens)*D%R_av(i)*(D%Thet(jj+1)-D%Thet(jj))
		enddo

		mindens=(massgrains-mass*Agrains)/(Agrains*D%R_av(i)*pi)
		if(mindens.gt.C(i,j)%dens) then
			if(mindens.lt.C(i,j+1)%dens) then
				C(i,j)%dens=mindens
			else
				C(i,j)%dens=C(i,j+1)%dens
			endif
		endif
	enddo
	masstot=0d0
	do j=1,D%nTheta-1
		masstot=masstot+C(i,j)%dens*C(i,j)%V
	enddo
	do j=1,D%nTheta-1
		C(i,j)%dens=C(i,j)%dens*masstot0/masstot
		C(i,j)%dens0=C(i,j)%dens
		C(i,j)%mass=C(i,j)%dens*C(i,j)%V
	enddo
	enddo
	
	return
	end
	
	


	subroutine compareCells(C1,C2,errT,errRho,errE)
	use Parameters
	IMPLICIT NONE
	type(Cell) C1(0:D%nR,0:D%nTheta),C2(0:D%nR,0:D%nTheta)
	real*8 errT,errRho,errE
	integer i,j,iT,ii
	
	errT=0d0
	errRho=0d0
	errE=0d0
	do i=1,D%nR-1
	do j=1,D%nTheta-1
		if(C1(i,j)%dT.gt.1d-200.and.C2(i,j)%dT.gt.1d-200) then
		errT=errT+
     & (C1(i,j)%T-C2(i,j)%T)**2/
     &      (C1(i,j)%dT**2+C2(i,j)%dT**2)
		errRho=errRho+
     & (C1(i,j)%dens-C2(i,j)%dens)**2/
     &      (C1(i,j)%ddens**2+C2(i,j)%ddens**2)
c		errE=errE+
c     & (C1(i,j)%mass*C1(i,j)%EJv-C2(i,j)%mass*C2(i,j)%EJv)**2/
c     &      ((C1(i,j)%mass*C1(i,j)%dEJv)**2+(C2(i,j)%mass*C2(i,j)%dEJv)**2)
		errE=errE+
     & (C1(i,j)%EJv-C2(i,j)%EJv)**2/
     &      (C1(i,j)%dEJv**2+C2(i,j)%dEJv**2)
		endif
	enddo
	enddo
	errT=errT/(real(D%nR-1)*real(D%nTheta-1)-1d0)
	errRho=sqrt(errRho/(real(D%nR-1)*real(D%nTheta-1)-1d0))
	errE=errE/(real(D%nR-1)*real(D%nTheta-1)-1d0)

	return
	end
	

	subroutine destroymaxtau(ii0)
	use Parameters
	real*8 r(D%nR),ct,tau,lr0,lr1
	real*8 Mtot,Mtot2,Fsc,p,p0,p1,Kext,wav,tauR(D%nR)
	logical escape,hitstar
	type(photon) phot
	integer i,j,l,n,ii,nl,iopac
	real*8 mintau(D%nR,D%nTheta)

c GFORTRAN syntax error (should also happen in ifort I guess)
c but that seems to be dead code anyway!
c	if(Grain(ii).maxtau.lt.0d0) return
	if(Grain(ii)%maxtau.lt.0d0) return

c	wav=0.55
c	phot%nr=1
c	tautot=0d0
c	tauR=0d0
c	do i=1,D%nR-1
c		r(i)=D%R_av(i)/AU
c		phot%x=r(i)
c		phot%y=0d0
c		phot%z=sqrt(D%R(D%nR)**2-phot%x**2)
c		ct=phot%z/D%R(D%nR)
c		do j=1,D%nTheta-1
c			if(ct.lt.D%Theta(j).and.ct.gt.D%Theta(j+1)) then
c				phot%i=D%nR-1
c				phot%j=j
c			endif
c		enddo
c		phot%edgenr=2
c		phot%onEdge=.true.
c		phot%vx=0d0
c		phot%vy=0d0
c		phot%vz=-1d0
c		tau=1d0
c		phot%E=1d0
c		call TraceStellar(phot,tau,escape,hitstar,.false.)
c		if(escape) phot%z=-1d0 * (D%R(D%nR)**2-phot%x**2)
c		tau1(i,l)=phot%z
c		enddo
c		tauR(i+1)=tauR(i)
c		do ii=1,ngrains
c		   do iopac=1,Grain(ii)%nopac
c			tauR(i+1)=tauR(i+1)+Grain(ii)%Kpstar(iopac)*C(i,D%nTheta-1)%dens*C(i,D%nTheta-1)%w(ii)*(D%R(i+1)-D%R(i))*AU
c		   enddo
c		enddo
c	enddo
c	
	
	return
	end
	

c-----------------------------------------------------------------------
c This subroutine outputs the cumulative radial tau for stellar
c radiation to the file filename
c-----------------------------------------------------------------------
	subroutine radialtaumidplane(filename)
	use Parameters
	IMPLICIT NONE
	character*500 filename
	real*8 tau(0:D%nR-1,0:ntau1_lam),lam0,wl1,wl2
	integer ilam1,ilam2,l,i,j,ii,iopac,jj

	jj=D%nTheta-1
	do l=0,ntau1_lam
		if(l.ne.0) then
			lam0=tau1_lam(l)
			do j=1,nlam-1
				if(lam0.ge.lam(j).and.lam0.le.lam(j+1)) then
					ilam1=j
					ilam2=j+1
					wl1=(lam(j+1)-lam0)/(lam(j+1)-lam(j))
					wl2=(lam0-lam(j))/(lam(j+1)-lam(j))
				endif
			enddo
			do ii=1,ngrains
				do iopac=1,Grain(ii)%nopac
					Grain(ii)%KextL(iopac)=Grain(ii)%Kext(iopac,ilam1)*wl1+
     &     Grain(ii)%Kext(iopac,ilam2)*wl2
				enddo
			enddo
			tau(0,l)=0d0
			do i=1,D%nR-1
				tau(i,l)=tau(i-1,l)
				do ii=1,ngrains
					do iopac=1,Grain(ii)%nopac
						tau(i,l)=tau(i,l)+(D%R(i+1)-D%R(i))*AU*
     &					C(i,jj)%dens*C(i,jj)%wopac(ii,iopac)*C(i,jj)%w(ii)*Grain(ii)%KextL(iopac)
					enddo
				enddo
			enddo
		else
			tau(0,l)=0d0
			do i=1,D%nR-1
				tau(i,l)=tau(i-1,l)
				do ii=1,ngrains
					do iopac=1,Grain(ii)%nopac
						tau(i,l)=tau(i,l)+(D%R(i+1)-D%R(i))*AU*
     &				C(i,jj)%dens*C(i,jj)%wopac(ii,iopac)*C(i,jj)%w(ii)*Grain(ii)%Kpstar(iopac)
					enddo
				enddo
			enddo
		endif
	enddo

	open(unit=80,file=filename,RECL=6000)
	do i=1,D%nR-1
		write(80,*) D%R(i+1),(tau(i,0:ntau1_lam))
	enddo
	close(unit=80)

	return
	end


c-----------------------------------------------------------------------
c This subroutine outputs the tau=1 height of the disk for stellar
c radiation to the file filename
c-----------------------------------------------------------------------
	subroutine tau1height(filename)
	use Parameters
	IMPLICIT NONE
	character*500 filename
	real*8 tau1(D%nR*10,0:10),r(D%nR*10),ct,tau,lr0,lr1
	real*8 Mtot,Mtot2,Fsc,p,p0,p1,Kext,wav(10),tauR(D%nR*10)
	logical escape,hitstar
	type(photon) phot
	integer i,j,l,n,ii,nl,iopac
	
	nl=ntau1_lam
	do i=1,nl
	   wav(i)=tau1_lam(i)
	enddo

	phot%nr=1
	tautot=0d0
	tauR=0d0
	do i=1,D%nR-1
		r(i)=D%R_av(i)/AU
		do l=0,nl
		phot%x=r(i)
		phot%y=0d0
		phot%z=sqrt(D%R(D%nR)**2-phot%x**2)
		ct=phot%z/D%R(D%nR)
		do j=1,D%nTheta-1
			if(ct.lt.D%Theta(j).and.ct.gt.D%Theta(j+1)) then
				phot%i=D%nR-1
				phot%j=j
			endif
		enddo
		phot%edgenr=2
		phot%onEdge=.true.
		phot%vx=0d0
		phot%vy=0d0
		phot%vz=-1d0
		tau=1d0
		if(l.ne.0) then
			phot%lam=wav(l)
			do j=1,nlam-1
				if(phot%lam.ge.lam(j).and.phot%lam.le.lam(j+1)) then
					phot%ilam1=j
					phot%ilam2=j+1
					phot%wl1=(lam(j+1)-phot%lam)/(lam(j+1)-lam(j))
					phot%wl2=(phot%lam-lam(j))/(lam(j+1)-lam(j))
				endif
			enddo
			phot%nu=nu(phot%ilam1)*phot%wl1+nu(phot%ilam2)*phot%wl2
			phot%E=1d0
			do ii=1,ngrains
				do iopac=1,Grain(ii)%nopac
					Grain(ii)%KabsL(iopac)=Grain(ii)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kabs(iopac,phot%ilam2)*phot%wl2
					Grain(ii)%KextL(iopac)=Grain(ii)%Kext(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kext(iopac,phot%ilam2)*phot%wl2
				enddo
			enddo
			call trace2d(phot,tau,escape,hitstar,.false.)
		else
			phot%E=1d0
			call TraceStellar(phot,tau,escape,hitstar,.false.)
		endif
c		if(escape) phot%z=0d0!-1d0/0d0!sqrt(D%R(D%nR)**2-phot%x**2)
		if(escape) phot%z=-1d0 * (D%R(D%nR)**2-phot%x**2)
		tau1(i,l)=phot%z
		enddo
		tauR(i+1)=tauR(i)
		do ii=1,ngrains
		   do iopac=1,Grain(ii)%nopac
			tauR(i+1)=tauR(i+1)+Grain(ii)%Kpstar(iopac)*C(i,D%nTheta-1)%dens*C(i,D%nTheta-1)%w(ii)*(D%R(i+1)-D%R(i))*AU
		   enddo
		enddo
	enddo
	open(unit=80,file=filename,RECL=6000)
	do i=1,D%nR-1
		write(80,*) r(i),(tau1(i,0:nl)),tauR(i)
	enddo
	close(unit=80)
		
	return
	end
	


c-----------------------------------------------------------------------
c This subroutine outputs the tau=1 height of the disk for stellar
c radiation when viewed from the central star
c to the file filename
c-----------------------------------------------------------------------
	subroutine tau1heightR(filename)
	use Parameters
	IMPLICIT NONE
	character*500 filename,file
	real*8 tau1z(D%nTheta,0:ntau1_lam),tau1x(D%nTheta,0:ntau1_lam),r(D%nTheta),tau1r(D%nTheta)
	real*8 ct,tau,lr0,lr1,dHrdr
	real*8 Mtot,Mtot2,Fsc,p,p0,p1,Kext
	logical escape,hitstar
	type(photon) phot
	integer i,j,l,n,ii,iopac

	integer it
	real*8 tt
	
	tt=1d0
	phot%nr=1
	tautot=0d0
	do i=1,D%nTheta
		do l=0,ntau1_lam
		r(i)=D%R(1)
		if(i.eq.D%nTheta) then
			ct=D%Theta(i-1)*0.01
		else
			ct=D%Theta(i)*0.99+D%Theta(i+1)*0.01
		endif
		if(ct.gt.1d0) ct=1d0
		phot%x=sqrt(1d0-ct**2)*r(i)
		phot%y=0d0
		phot%z=ct*r(i)

		phot%edgenr=1
		phot%onEdge=.true.
		phot%vx=phot%x/r(i)
		phot%vy=0d0
		phot%vz=phot%z/r(i)
		tau=tt
		if(l.ne.0) then
			phot%lam=tau1_lam(l)
		else
			phot%lam=0.55d0
		endif
		phot%i=1
		if(i.eq.D%nTheta) then
			phot%j=D%nTheta-1
		else
			phot%j=i
		endif
		do j=1,nlam-1
			if(phot%lam.ge.lam(j).and.phot%lam.le.lam(j+1)) then
				phot%ilam1=j
				phot%ilam2=j+1
				phot%wl1=(lam(j+1)-phot%lam)/(lam(j+1)-lam(j))
				phot%wl2=(phot%lam-lam(j))/(lam(j+1)-lam(j))
			endif
		enddo
		phot%nu=nu(phot%ilam1)*phot%wl1+nu(phot%ilam2)*phot%wl2
		phot%E=1d0
		do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				Grain(ii)%KabsL(iopac)=Grain(ii)%Kabs(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kabs(iopac,phot%ilam2)*phot%wl2
				Grain(ii)%KextL(iopac)=Grain(ii)%Kext(iopac,phot%ilam1)*phot%wl1+
     &     Grain(ii)%Kext(iopac,phot%ilam2)*phot%wl2
			enddo
		enddo
		if(filename.ne.' '.and.l.ne.0) then
			call trace2d(phot,tau,escape,hitstar,.false.)
		else
			call TraceStellar(phot,tau,escape,hitstar,.false.)
			tau1r(i)=sqrt(phot%x**2+phot%z**2)
		endif
		if(escape) phot%z=0d0!-1d0/0d0!sqrt(D%R(D%nR)**2-phot%x**2)
		tau1x(i,l)=phot%x
		tau1z(i,l)=phot%z
		enddo
	enddo

	if(.not.allocated(muRad)) allocate(muRad(D%nTheta))
	muRad(D%nTheta)=1d0
	do i=D%nTheta-1,1,-1
		if(tau1z(i+1,0).gt.0d0.and.tau1z(i,0).gt.0d0) then
			dHrdr=(tau1z(i+1,0)/tau1x(i+1,0)-tau1z(i,0)/tau1x(i,0))/(tau1x(i+1,0)-tau1x(i,0))
		else
			dHrdr=0d0
		endif
		if(tau1z(i+1,0).gt.0d0.and.tau1z(i,0).gt.0d0) then
			muRad(i)=abs(sin(atan(tau1x(i,0)/tau1z(i,0))+atan(dHrdr)-pi/2d0))
		else
			muRad(i)=0d0
		endif
	enddo

c	write(file,'(a,i4)') filename(1:len_trim(filename)),it+1000
	if(filename.ne.' ') then
		open(unit=80,file=filename,RECL=6000)
		do i=1,D%nTheta
			write(80,*) tau1x(i,0),tau1z(i,0),muRad(i),(tau1x(i,l),tau1z(i,l),l=1,ntau1_lam)
		enddo
		close(unit=80)
	endif

	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	subroutine determine_inc_tau()
	use Parameters
	IMPLICIT NONE
	real*8 tautheta(D%nTheta)
	real*8 lam0,wl1,wl2,Kext
	integer i,j,l,n,ii,iopac,il1,il2

	tautheta=0d0
	lam0=0.55
	do j=1,nlam-1
		if(lam0.ge.lam(j).and.lam0.le.lam(j+1)) then
			il1=j
			il2=j+1
			wl1=(lam(j+1)-lam0)/(lam(j+1)-lam(j))
			wl2=(lam0-lam(j))/(lam(j+1)-lam(j))
		endif
	enddo
	do j=1,D%nTheta
		do i=1,D%nR-1
		Kext=0d0
		do ii=1,ngrains
			do iopac=1,Grain(ii)%nopac
				Kext=Kext+C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)*
     &					(Grain(ii)%Kext(iopac,il1)*wl1+Grain(ii)%Kext(iopac,il2)*wl2)
			enddo
		enddo
		tautheta(j)=tautheta(j)+Kext*(D%R(i+1)-D%R(i))*AU*C(i,j)%dens
		enddo
	enddo

	D%IA=0.9999d0*pi/2d0
	do j=1,D%nTheta-2
		if(tau_inc.gt.tautheta(j).and.tau_inc.le.tautheta(j+1)) then
			D%IA=D%thet(j)+(tau_inc-tautheta(j))*(D%thet(j+1)-D%thet(j))/(tautheta(j+1)-tautheta(j))
		endif
	enddo
	D%IA=D%IA*180d0/pi
	write(*,'("Sightline optical depth:       ",f6.2)') tau_inc
	write(9,'("Sightline optical depth:       ",f6.2)') tau_inc
	write(*,'("Adjusted inclination angle to: ",f6.2)') D%IA
	write(9,'("Adjusted inclination angle to: ",f6.2)') D%IA
	
	return
	end
	





c-----------------------------------------------------------------------
c This subroutine outputs the pressure and density scaleheight of the
c disk using the temperature and density structure. It also stores the
c Toomre Q parameter. Note that the temperature is needed for
c this. Output is written to filename. 
c
c The scaleheight is defined as the height at which the pressure/density
c drops off by a factor e^-0.5
c-----------------------------------------------------------------------
	subroutine scaleheight(filename)
	use Parameters
	IMPLICIT NONE
	character*500 filename
	real*8 sh(2,D%nR),Sig(D%nR),Q,shgas,shdust
	real*8 ToomreQ ! function
	logical shset(2)
	integer i,j
	
	do i=1,D%nR-1
	   call calcscaleheight(i,shgas,shdust)
	   sh(1,i)=shgas
	   sh(2,i)=shdust

	   Sig(i)=0d0
	   do j=1,D%nTheta-1
	      Sig(i)=Sig(i)+C(i,j)%mass/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)
	   enddo
	enddo

	open(unit=80,file=filename,RECL=1000)
	write(80,'("# R[AU], pressure scaleheight, density scaleheight, Toomre Q")')
	do i=1,D%nR-1
!		Q=(sh(1,i)*AU*D%Mstar)/(D%R_av(i)**3*pi*Sig(i)*gas2dust)
	   Q=ToomreQ(i)
	   write(80,*) D%R_av(i)/AU,sh(1,i),sh(2,i),Q
	enddo
	close(unit=80)
	return
	end

c-----------------------------------------------------------------------
c This subroutine calculates the pressure scaleheight of the disk using
c the temperature and density structure at one radius.
c The scaleheight is defined as the height at which the pressure/density
c drops off by a factor e^-0.5
c-----------------------------------------------------------------------
	subroutine calcscaleheight(i,shgas,shdust)
	use Parameters
	IMPLICIT NONE

	real*8 z(D%nTheta)
	real*8 p,p0,p1,d0,d1,dens,shgas,shdust
	logical shset(2)
	integer i,j
	
	p0=C(i,D%nTheta-1)%gasdens*C(i,D%nTheta-1)%T
	p1=p0
	p0=p0/exp(0.5d0)	! p(z)=p0*e^{-z^2/2H_p^2} -> p(H_p)=p0*e^-0.5
	d0=C(i,D%nTheta-1)%gasdens
	d1=d0
	d0=d0/exp(0.5d0)	! p(z)=p0*e^{-z^2/2H_p^2} -> p(H_p)=p0*e^-0.5
	shset=.false.
	z(D%nTheta-1)=D%R_av(i)*cos(D%theta_av(D%nTheta-1))/AU
	do j=D%nTheta-2,1,-1
	   z(j)=D%R_av(i)*cos(D%theta_av(j))/AU
	   p=C(i,j)%gasdens*C(i,j)%T
	   if(p.lt.p0.and..not.shset(1)) then
	      shgas=(z(j+1)+(z(j)-z(j+1))*(p1-p0)/(p1-p))
	      shset(1)=.true.
	   endif
	   p1=p
	   dens=C(i,j)%gasdens
	   if(dens.lt.d0.and..not.shset(2)) then
	      shdust=(z(j+1)+(z(j)-z(j+1))*(d1-d0)/(d1-dens))
	      shset(2)=.true.
	   endif
	   d1=dens
	enddo
	
	return
	end

c-----------------------------------------------------------------------
c This subroutine calculate the Toomre Q parameter for the disk.
c It requires the temperature to be set! (through the scaleheight)
c-----------------------------------------------------------------------
	real*8 function ToomreQ(i)
	use Parameters
	IMPLICIT NONE
	integer i,j
	real*8 sh(2,D%nR),Sig(D%nR),shgas,shdust
	
	! integrate surface density
	Sig(i)=0d0
	do j=1,D%nTheta-1
	   Sig(i)=Sig(i)+C(i,j)%mass/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)
	enddo
	
	call calcscaleheight(i,shgas,shdust)
	sh(1,i)=shgas
	sh(2,i)=shdust

        ! Calc Toomre Q parameter
	ToomreQ=(sh(1,i)*AU*D%Mstar)/(D%R_av(i)**3*pi*Sig(i)*gas2dust)

	return
	end

c-----------------------------------------------------------------------
c  This subroutine returns the temperature from the module diskstruct if
c  it is defined, otherwise it returns C(*,*)%T
c-----------------------------------------------------------------------
	
	subroutine GetTemp(ReturnTemp)
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	real*8 ReturnTemp(1:D%nR-1,1:D%nTheta-1)
		
	if(allocated(Temp)) then
	   ReturnTemp(1:D%nR-1,1:D%nTheta-1)=Temp(1:D%nR-1,1:D%nTheta-1)
	else
	   do i=1,D%nR-1
	      do j=1,D%nTheta-1
		 ReturnTemp(i,j)=C(i,j)%T
	      enddo
	   enddo
	endif

	end subroutine

c-----------------------------------------------------------------------
c This subroutine computes the optical depth through the midplane of
c the disk. It does this at lambda=0.55 micron.
c It also computes the radial optical depth in the first cell.
c In order for the radiative transfer to be accurate this needs to
c be small.
c-----------------------------------------------------------------------
	real*8 function radialtau(lam0,tau1,ntau1,j0)
	use Parameters
	IMPLICIT NONE
	real*8 tau,Kext,wl1,wl2,lam0,tau1,tt(D%nR,D%nTheta)
	integer i,j,ii,ilam1,ilam2,ntau1,j0,iopac
	
	do j=1,nlam-1
		if(lam0.ge.lam(j).and.lam0.le.lam(j+1)) then
			ilam1=j
			ilam2=j+1
			wl1=(lam(j+1)-lam0)/(lam(j+1)-lam(j))
			wl2=(lam0-lam(j))/(lam(j+1)-lam(j))
		endif
	enddo
	j=D%nTheta-1
	
	do j=1,D%nTheta-1
	ntau1=0
	tau=0d0
	tau1=-1d0
	do i=1,D%nR-1
		Kext=0d0
		do ii=1,ngrains
		   do iopac=1,Grain(ii)%nopac
			Kext=Kext+(Grain(ii)%Kext(iopac,ilam1)*C(i,j)%wopac(ii,iopac)*C(i,j)%w(ii))*wl1+
     &				  (Grain(ii)%Kext(iopac,ilam2)*C(i,j)%wopac(ii,iopac)*C(i,j)%w(ii))*wl2
		   enddo
		enddo
		tau=tau+(D%R(i+1)-D%R(i))*C(i,j)%dens*Kext*AU
		if(tau.gt.0.2.and.tau1.lt.0d0) tau1=(D%R(i+1)-D%R(i))*C(i,j)%dens*Kext*AU
		if(tau.gt.0.2d0.and.(tau-(D%R(i+1)-D%R(i))*C(i,j)%dens*Kext*AU).lt.1d0) then
c			write(*,'("Cell nr. ",i3,": ",f10.3)') i,tau
c			write(9,'("Cell nr. ",i3,": ",f10.3)') i,tau
			ntau1=ntau1+1
		endif
		tt(i,j)=tau
	enddo
	if(j.eq.j0) radialtau=tau
	enddo

	return
	end

c-----------------------------------------------------------------------
c This subroutine computes the density structure so that the scaleheight
c corresponds to a powerlaw in (cylindrical) radius.
c The surface density has to be set for this.
c-----------------------------------------------------------------------
	subroutine SetScaleHeight(h0,h1,pow)
	use Parameters
	IMPLICIT NONE
	real*8 M0,M1,r,z
	real*8 h0,h1,h,pow
	real*8 rho(D%nTheta-1)
	integer i,j

	do i=1,D%nR-1
		M0=0d0
		do j=1,D%nTheta-1
			M0=M0+C(i,j)%dens*C(i,j)%V
		enddo
		do j=1,D%nTheta-1
			r=D%R_av(i)*sin(D%theta_av(j))/AU
			z=D%R_av(i)*cos(D%theta_av(j))/AU
			h=(h1-h0*D%R(1)**pow)*exp(-sqrt(10d0*abs(r-D%R(1))/D%R(1)))+h0*r**pow
			h=h0*r**pow+h1*(D%R_av(i)/D%R_av(1))**(-5d0)
			rho(j)=exp(-(z/h)**2)
		enddo
		M1=0d0
		do j=1,D%nTheta-1
			M1=M1+rho(j)*C(i,j)%V
		enddo
		do j=1,D%nTheta-1
			C(i,j)%dens=rho(j)*M0/M1
			if(C(i,j)%dens.lt.1d-50) C(i,j)%dens=1d-50
			C(i,j)%mass=C(i,j)%dens*C(i,j)%V
			C(i,j)%dens0=C(i,j)%dens
		enddo
	enddo

	return
	end



c-----------------------------------------------------------------------
c This subroutine destroys all dust with temperatures above Tdes
c If no dust is destroyed tdes_iter is set to .false.
c 2008-03-26: MM added a check for the case that w0=0d0
c
c 2008-09-05: We could use the most accurate solution of the evaporation:
c						ln(rho)=A+B/T+ln(T)
c             if we turn to the partial pressure in stead of using the
c			  evaporation temperature. (see paper on my desk!!)
c-----------------------------------------------------------------------
	subroutine DestroyDustT(C1,Rdes,Rmin,maxdtau,BW)
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	integer ii,k,l,ntau,itemp,iT,i2,itau,iter,viter,niter,iopac
	real*8 Rdes,Tpowdes,tot,Rmin,f,T,R,determineT,fmax,Kext,Tevap
	logical dust,jumpdepend,BW
	real*8 f_new(0:D%nR,ngrains),f_prev(0:D%nR,ngrains)
	real*8 dens1(0:D%nR),w1(0:D%nR,ngrains),maxdtau,fw(0:D%nR,ngrains)
	real*8 f_gas(0:D%nR,1:D%nTheta,ngrains)
c MK: type Cell is defined in Modules.f
	type(Cell) C1(0:D%nR,0:D%nTheta)
	real*8 spec(nlam),kp,E,frac,lr,gasfrac,dens,maxT,minT,pow,powmin,powmax
	real*8 Rtau1,tau,tau0
	type(photon) phot
	real*8 fitA,fitB,determinegasfrac,Tthin,W,determineTP
c Mihkelexp - introduce and set Tzerocounter to count number of zero Monte Carlo temperature points:
	integer Tzerocounter,gradientcount,istrongest,number_invalid
	real*8 gf,gfstrongest
	logical gradienterror,destroy
	
	real*8 minf(ngrains)

	Tzerocounter=0

	if(.not.allocated(iRfirst)) allocate(iRfirst(1:ngrains,0:D%nTheta))
	iRfirst(1:ngrains,0:D%nTheta)=D%nR-1

	call DestroyDustR()		
	
	write(*,'("Destroying high temperature dust")')
	write(9,'("Destroying high temperature dust")')

	call determineRdes()
	
	dust=.false.
	Rdes=D%R(D%nR)
	do j=D%nTheta-1,1,-1

	call tellertje(D%nTheta-j,D%nTheta-1)

	pow=1d0
	powmin=1d0
	powmax=-1d0
	ntau=0
	tau=0d0
	Rtau1=D%R(D%nR)
	do i=1,D%nR-1
		tau0=C1(i,j)%KextLRF*C1(i,j)%dens*AU*(D%R(i+1)-D%R(i))
		if((tau+tau0).gt.1d0) then
			Rtau1=D%R(i)+(D%R(i+1)-D%R(i))*(1d0-tau)/tau0
			goto 1
		endif
		tau=tau+tau0
	enddo
1	continue
	ntau=ntau+1
	do i=0,D%nR-1
c Mihkelexp
c If Monte Carlo temperature is less than 1K, set gas fraction to zero:
		phot%i=i
		phot%j=j
		istrongest=1
		gfstrongest=1d200
		if(i.ne.0) then
			r=D%R(i)*AU
		else
			r=D%R(1)*AU
		endif
		W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))
		do ii=1,ngrains
			spec(1:nlam)=0d0
			do iopac=1,Grain(ii)%nopac
			   spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)*C1(i,j)%wopac(ii,iopac)/(pi*D%Rstar**2)
			enddo
			call integrate(spec,phot%E)
			phot%E=W*phot%E
			Tthin=determineTP(phot,ii)
			gf=determinegasfrac(Tthin,i,j,ii)
			if(gf.lt.gfstrongest) then
				istrongest=ii
				gfstrongest=gf
			endif
		enddo
		do ii=1,ngrains
			if(Grain(ii)%qhp) then
				T=C1(i,j)%Tqhp(Grain(ii)%qhpnr)
			else if(tcontact.or.C1(i,j)%diff) then
				T=C1(i,j)%T
			else
				T=C1(i,j)%TP(ii)
			endif	
			Tevap=T
			minT=0d0
			maxT=real(TMAX)*dT
10			continue!			do iter=1,100
				if(determinegasfrac(Tevap,i,j,ii).gt.1d0) then
					maxT=Tevap
					Tevap=(Tevap+minT)/2d0
				else
					minT=Tevap
					Tevap=(Tevap+maxT)/2d0
				endif
			if(abs(maxT-minT).gt.1d0) goto 10!			enddo
c			niter=1
c			if(viscous.and.T.ge.Tevap) then
c				if(C1(i,j)%EviscDirect(ii).ne.0d0) niter=2
c			endif
c			do viter=1,niter
			if(C1(i,j)%w0(ii).ne.0d0) then
			f_prev(i,ii)=1d0 - C1(i,j)%dens*C1(i,j)%w(ii) / (C1(i,j)%dens0*C1(i,j)%w0(ii))
			if(Grain(ii)%qhp) then
				T=C1(i,j)%Tqhp(Grain(ii)%qhpnr)
			else! if(.not.viscous.or.T.lt.Tevap) then
				if(tcontact.or.C1(i,j)%diff) then
					T=C1(i,j)%T
				else
					T=C1(i,j)%TP(ii)
				endif	
c			else
c				if(C1(i,j)%EviscDirect(ii).eq.0d0) then
c					if(tcontact.or.C1(i,j)%diff) then
c						T=C1(i,j)%T
c					else
c						T=C1(i,j)%TP(ii)
c					endif
c				else
c				if(tcontact) then
c					if(viter.eq.1) then
c						phot%E=C1(i,j)%EJv
c     &	-C1(i,j)%EviscDirect(ii)/(C1(i,j)%mass*C(i,j)%w(ii))
c     &	+C1(i,j)%EviscDirect(ii)/(C1(i,j)%dens0*C(i,j)%w0(ii)*C1(i,j)%V)
c					else
c						phot%E=C(i,j)%EJv
c     &	-C1(i,j)%EviscDirect(ii)/(C1(i,j)%mass*C1(i,j)%w(ii))
c     &	+C1(i,j)%EviscDirect(ii)/((1d0-f_new(i,ii))*C1(i,j)%dens0*C(i,j)%w0(ii)*C1(i,j)%V)
c					endif
c					T=determineT(phot)
c				else
c					if(viter.eq.1) then
c						phot%E=C1(i,j)%EJvP(ii)
c     &	-C1(i,j)%EviscDirect(ii)/(C1(i,j)%mass*C1(i,j)%w(ii))
c     &	+C1(i,j)%EviscDirect(ii)/(C1(i,j)%dens0*C1(i,j)%w0(ii)*C1(i,j)%V)
c					else
c						phot%E=C(i,j)%EJvP(ii)
c     &	-C1(i,j)%EviscDirect(ii)/(C1(i,j)%mass*C1(i,j)%w(ii))
c     &	+C1(i,j)%EviscDirect(ii)/((1d0-f_new(i,ii))*C1(i,j)%dens0*C1(i,j)%w0(ii)*C1(i,j)%V)
c					endif
c					T=determineTP(phot,ii)
c				endif	
c				endif
			endif
		
			if(T.lt.1) then
				f_new(i,ii)=0d0
				Tzerocounter=Tzerocounter+1
			else
				f_new(i,ii)=determinegasfrac(T,i,j,ii)
			endif
			else
				f_new(i,ii)=0d0
				f_prev(i,ii)=0d0
			endif
			if(f_new(i,ii).gt.1d0.and.viter.ne.niter) f_new(i,ii)=1d0
c			enddo
			if(Grain(ii)%qhp) then
				T=C1(i,j)%Tqhp(Grain(ii)%qhpnr)
			else if(tcontact.or.C1(i,j)%diff) then
				T=C1(i,j)%T
			else
				T=C1(i,j)%TP(ii)
			endif	
			if(f_new(i,ii).lt.1d0.and.D%R_av(i)/AU.lt.Rdes) Rdes=D%R_av(i)/AU

			if(i.ne.0) then
				r=D%R(i)*AU
			else
				r=D%R(1)*AU
			endif
			W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))
			phot%E=0d0
			do iopac=1,Grain(ii)%nopac
			   phot%E=phot%E+W*Grain(ii)%Kpabsstar(iopac)*C1(i,j)%wopac(ii,iopac)*D%Lstar/(pi*D%Rstar**2)
			enddo
			Tthin=determineTP(phot,ii)

			fw(i,ii)=abs((Tevap-T)/(Tevap+T))
c			if(ii.eq.istrongest) 
			fw(i,ii)=fw(i,ii)*abs((Tevap-Tthin)/(Tevap+Tthin))
			fw(i,ii)=f_weight*(fw(i,ii)**pow)
			if(f_new(i,ii).lt.0d0) f_new(i,ii)=0d0
			if(T.gt.Tevap) then
				fw(i,ii)=f_weight*(1d0-exp(-(20d0*(Tevap-T)/Tevap)**2))
				f_new(i,ii)=1d0
c			else
c				fw(i,ii)=f_weight*(1d0-exp(-((Tevap-T)/200d0)**2))
			endif
			if(fw(i,ii).gt.f_weight) fw(i,ii)=f_weight
			if(Grain(ii)%tdes_fast.gt.fw(i,ii)) 
     &	fw(i,ii)=(1d0-Grain(ii)%tdes_fast)*fw(i,ii)+Grain(ii)%tdes_fast
		enddo
	enddo
	f_new(0,1:ngrains)=1d0
	f_prev(0,1:ngrains)=1d0

	do i=D%nR-1,0,-1
		do ii=1,ngrains
			f_gas(i,j,ii)=fw(i,ii)*f_new(i,ii)+(1d0-fw(i,ii))*f_prev(i,ii)
			if(D%R(i+1).lt.Rthindes(j)) f_gas(i,j,ii)=1d0
			if(f_gas(i,j,ii).gt.1d0) f_gas(i,j,ii)=1d0
			if(f_gas(i,j,ii).lt.0d0) f_gas(i,j,ii)=0d0
			if(i.lt.D%nR-1.and.Grain(ii)%force_vert_gf) then
				if(f_gas(i,j,ii).lt.f_gas(i+1,j,ii).and.(D%R_av(i+1)/AU).lt.Grain(ii)%maxrad) f_gas(i,j,ii)=f_gas(i+1,j,ii)
			endif
			if(j.lt.D%nTheta-1.and.Grain(ii)%force_vert_gf) then
				if(f_gas(i,j,ii).lt.f_gas(i,j+1,ii)) f_gas(i,j,ii)=f_gas(i,j+1,ii)
			endif
			w1(i,ii)=(1d0-f_gas(i,j,ii))*C1(i,j)%w0(ii)
		enddo
		tot=0d0
		do ii=1,ngrains
c			Sum of all the solid density fractions:
			tot=tot+w1(i,ii)
		enddo
c----------------------
		dens1(i)=C1(i,j)%dens+(tot*C1(i,j)%dens0-C1(i,j)%dens)
		if(dens1(i).gt.C1(i,j)%dens0) dens1(i)=C1(i,j)%dens0
		if(dens1(i).lt.1d-50) then
			dens1(i)=1d-50
		endif
	enddo

	if(ntau.gt.500) goto 2
	tau=0d0
	do i=1,D%nR-1
		tau0=C1(i,j)%KextLRF*dens1(i)*AU*(D%R(i+1)-D%R(i))
		if(D%R(i+1).gt.Rtau1) then
			tau=tau+C1(i,j)%KextLRF*dens1(i)*AU*(Rtau1-D%R(i))
			if(tau.gt.(1d0+maxdtau)) then
				powmin=pow
				if(powmax.gt.0d0) then
					pow=(pow+powmax)/2d0
				else
					pow=pow*2d0
				endif
				goto 1
			else if(pow.gt.(powmin*1.01)) then
				powmax=pow
				pow=(pow+powmin)/2d0
				goto 1
			else
				goto 2
			endif
		endif
		tau=tau+tau0
	enddo

2	continue

	do i=D%nR-1,0,-1
		do ii=1,ngrains
			f_gas(i,j,ii)=fw(i,ii)*f_new(i,ii)+(1d0-fw(i,ii))*f_prev(i,ii)
			if(D%R(i+1).lt.Rthindes(j)) f_gas(i,j,ii)=1d0
			if(f_gas(i,j,ii).gt.1d0) f_gas(i,j,ii)=1d0
			if(f_gas(i,j,ii).lt.0d0) f_gas(i,j,ii)=0d0
			if(i.lt.D%nR-1.and.Grain(ii)%force_vert_gf) then
				if(f_gas(i,j,ii).lt.f_gas(i+1,j,ii).and.(D%R_av(i+1)/AU).lt.Grain(ii)%maxrad) f_gas(i,j,ii)=f_gas(i+1,j,ii)
			endif
			if(j.lt.D%nTheta-1.and.Grain(ii)%force_vert_gf) then
				if(f_gas(i,j,ii).lt.f_gas(i,j+1,ii)) f_gas(i,j,ii)=f_gas(i,j+1,ii)
			endif
			C1(i,j)%w(ii)=(1d0-f_gas(i,j,ii))*C1(i,j)%w0(ii)
			if(f_gas(i,j,ii).lt.1d0.and.i.lt.iRfirst(ii,j)) iRfirst(ii,j)=i
		enddo
		tot=0d0
		do ii=1,ngrains
c			Sum of all the solid density fractions:
			tot=tot+C1(i,j)%w(ii)
		enddo
c----------------------
		if(tot.gt.1d-50) then
			C1(i,j)%w(1:ngrains)=C1(i,j)%w(1:ngrains)/tot
		else
			C1(i,j)%w(1:ngrains)=C1(i,j)%w0(1:ngrains)
		endif
c----------------------
		C1(i,j)%dens=C1(i,j)%dens+(tot*C1(i,j)%dens0-C1(i,j)%dens)
		if(C1(i,j)%dens.gt.C1(i,j)%dens0) C1(i,j)%dens=C1(i,j)%dens0
		if(C1(i,j)%dens.lt.1d-50) then
			C1(i,j)%w(1:ngrains)=C1(i,j)%w0(1:ngrains)
			C1(i,j)%dens=1d-50
			C1(i,j)%gasfrac=1d0
		else
			C1(i,j)%gasfrac=1d0-C1(i,j)%dens/C1(i,j)%dens0
		endif

		if(number_invalid(C1(i,j)%dens).ne.0) then
			C1(i,j)%dens=C1(i,j)%dens0
			C1(i,j)%gasfrac=0d0
		endif
		do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			if(number_invalid(C1(i,j)%wopac(ii,iopac)).ne.0) print*,i,j,ii,iopac
		enddo
		enddo
c	if(j.eq.D%nTheta-1) print*,C1(i,j)%gasfrac
c scale the mass in the cell
		C1(i,j)%mass=C1(i,j)%dens*C1(i,j)%V
	enddo

	enddo

c Mihkelexp
c Output number of zero-temperature Monte Carlo points:
	if(Tzerocounter.gt.0) then
	write(*,'("Number of Monte Carlo temperature points below threshold:")')
	write(*,'(i10)') Tzerocounter
	write(9,'("Number of Monte Carlo temperature points below threshold:")')
	write(9,'(i10)') Tzerocounter
	else
	write(*,'("No Monte Carlo temperature points below threshold value.")')
	write(9,'("No Monte Carlo temperature points below threshold value.")')
	endif
c end Mihkelexp

	Rmin=Rdes
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	call CheckCells()

	minf(1:ngrains)=1d0
	j=D%nTheta-1
	do i=1,D%nR-1
		do ii=1,ngrains
			if(tcontact) then
				T=C1(i,j)%T
			else
				T=C1(i,j)%TP(ii)
			endif
			f_gas(i,j,ii)=1d0 - C1(i,j)%dens*C1(i,j)%w(ii) / (C1(i,j)%dens0*C1(i,j)%w0(ii))
			if(f_gas(i,j,ii).lt.minf(ii).and.T.gt.1000d0) minf(ii)=f_gas(i,j,ii)
		enddo
	enddo
c	print*,(1d0-minf(1:ngrains))

	call DestroyDustR()		

	return
	end
	
	



c-----------------------------------------------------------------------
c This subroutine destroys all dust with temperatures above Tdes
c If no dust is destroyed tdes_iter is set to .false.
c In contrast with the above subroutine this subroutine does direct destruction
c and recondensation without taking care with the changes.
c 2008-03-26: MM added a check for the case that w0=0d0
c
c 2008-09-05: We could use the most accurate solution of the evaporation:
c						ln(rho)=A+B/T+ln(T)
c             if we turn to the partial pressure in stead of using the
c			  evaporation temperature. (see paper on my desk!!)
c-----------------------------------------------------------------------
	subroutine DestroyDustTDirect(C1,Rdes,Rmin,maxdtau,BW)
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	integer ii,k,l,ntau,itemp,iT,i2,itau,iter,iopac
	real*8 Rdes,Tpowdes,tot,Rmin,f,T,R,determineT,fmax,Kext,Tevap
	logical dust,jumpdepend,BW
	real*8 f_new(0:D%nR,ngrains),f_prev(0:D%nR,ngrains)
	real*8 dens1(0:D%nR),w1(0:D%nR,ngrains),maxdtau,fw(0:D%nR,ngrains)
	real*8 f_gas(0:D%nR,1:D%nTheta,ngrains)
c MK: type Cell is defined in Modules.f
	type(Cell) C1(0:D%nR,0:D%nTheta)
	real*8 spec(nlam),kp,E,frac,lr,gasfrac,dens,maxT,minT,pow,powmin,powmax
	real*8 Rtau1,tau,tau0
	type(photon) phot
	real*8 fitA,fitB,determinegasfrac,Tthin,W,determineTP
c Mihkelexp - introduce and set Tzerocounter to count number of zero Monte Carlo temperature points:
	integer Tzerocounter,gradientcount
	logical gradienterror,destroy

	Tzerocounter=0

	if(.not.allocated(iRfirst)) allocate(iRfirst(1:ngrains,0:D%nTheta))
	iRfirst(1:ngrains,0:D%nTheta)=D%nR-1

	call DestroyDustR()		
	
	write(*,'("Destroying high temperature dust")')
	write(9,'("Destroying high temperature dust")')

	call determineRdes()
	
	dust=.false.
	Rdes=D%R(D%nR)
	do j=D%nTheta-1,1,-1

	call tellertje(D%nTheta-j,D%nTheta-1)

	pow=1d0
	powmin=1d0
	powmax=-1d0
	ntau=0
	tau=0d0
	Rtau1=D%R(D%nR)
	do i=1,D%nR-1
		tau0=C1(i,j)%KextLRF*C1(i,j)%dens*AU*(D%R(i+1)-D%R(i))
		if((tau+tau0).gt.1d0) then
			Rtau1=D%R(i)+(D%R(i+1)-D%R(i))*(1d0-tau)/tau0
			goto 1
		endif
		tau=tau+tau0
	enddo
1	continue
	ntau=ntau+1
	do i=0,D%nR-1
c Mihkelexp
c If Monte Carlo temperature is less than 1K, set gas fraction to zero:
		do ii=1,ngrains
			if(C1(i,j)%w0(ii).ne.0d0) then
			f_prev(i,ii)=1d0 - C1(i,j)%dens*C1(i,j)%w(ii) / (C1(i,j)%dens0*C1(i,j)%w0(ii))
			if(Grain(ii)%qhp) then
				T=C1(i,j)%Tqhp(Grain(ii)%qhpnr)
			else if(tcontact.or.C1(i,j)%diff) then
				T=C1(i,j)%T
			else
				T=C1(i,j)%TP(ii)
			endif	
		
			if(T.lt.1) then
				f_new(i,ii)=0d0
				Tzerocounter=Tzerocounter+1
			else
				f_new(i,ii)=determinegasfrac(T,i,j,ii)
			endif
			else
				f_prev(i,ii)=0d0
				f_new(i,ii)=0d0
			endif
			if(f_new(i,ii).lt.1d0.and.D%R_av(i)/AU.lt.Rdes) Rdes=D%R_av(i)/AU

			if(f_new(i,ii).lt.0d0) f_new(i,ii)=0d0
			if(f_new(i,ii).gt.1d0) f_new(i,ii)=1d0
		enddo
	enddo
	f_new(0,1:ngrains)=1d0
	f_prev(0,1:ngrains)=1d0

	do i=D%nR-1,0,-1
		do ii=1,ngrains
c			f_gas(i,j,ii)=f_weight*f_new(i,ii)+(1d0-f_weight)*f_prev(i,ii)
			f_gas(i,j,ii)=f_new(i,ii)
c			if(D%R(i+1).lt.Rthindes(j)) f_gas(i,j,ii)=1d0
			if(f_gas(i,j,ii).gt.1d0) f_gas(i,j,ii)=1d0
			if(f_gas(i,j,ii).lt.0d0) f_gas(i,j,ii)=0d0
c			if(i.lt.D%nR-1) then
c				if(f_gas(i,j,ii).lt.f_gas(i+1,j,ii)) f_gas(i,j,ii)=f_gas(i+1,j,ii)
c			endif
c			if(j.lt.D%nTheta-1) then
c				if(f_gas(i,j,ii).lt.f_gas(i,j+1,ii)) f_gas(i,j,ii)=f_gas(i,j+1,ii)
c			endif
			C1(i,j)%w(ii)=(1d0-f_gas(i,j,ii))*C1(i,j)%w0(ii)
			if(f_gas(i,j,ii).lt.1d0.and.i.lt.iRfirst(ii,j)) iRfirst(ii,j)=i
		enddo
		tot=0d0
		do ii=1,ngrains
c			Sum of all the solid density fractions:
			tot=tot+C1(i,j)%w(ii)
		enddo
c----------------------
		if(tot.gt.1d-50) then
			C1(i,j)%w(1:ngrains)=C1(i,j)%w(1:ngrains)/tot
		else
			C1(i,j)%w(1:ngrains)=C1(i,j)%w0(1:ngrains)
		endif
c----------------------
		C1(i,j)%dens=C1(i,j)%dens+(tot*C1(i,j)%dens0-C1(i,j)%dens)
		if(C1(i,j)%dens.gt.C1(i,j)%dens0) C1(i,j)%dens=C1(i,j)%dens0
		if(C1(i,j)%dens.lt.1d-50) then
			C1(i,j)%w(1:ngrains)=C1(i,j)%w0(1:ngrains)
			C1(i,j)%dens=1d-50
			C1(i,j)%gasfrac=1d0
		else
			C1(i,j)%gasfrac=1d0-C1(i,j)%dens/C1(i,j)%dens0
		endif

c	if(j.eq.D%nTheta-1) print*,C1(i,j)%gasfrac
c scale the mass in the cell
		C1(i,j)%mass=C1(i,j)%dens*C1(i,j)%V
	enddo

	enddo

c Mihkelexp
c Output number of zero-temperature Monte Carlo points:
	if(Tzerocounter.gt.0) then
	write(*,'("Number of Monte Carlo temperature points below threshold:")')
	write(*,'(i10)') Tzerocounter
	write(9,'("Number of Monte Carlo temperature points below threshold:")')
	write(9,'(i10)') Tzerocounter
	else
	write(*,'("No Monte Carlo temperature points below threshold value.")')
	write(9,'("No Monte Carlo temperature points below threshold value.")')
	endif
c end Mihkelexp

	Rmin=Rdes
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')

	call CheckCells()

	call DestroyDustR()		

	return
	end
	
	



c-----------------------------------------------------------------------
c This subroutine destroys all dust outside the given maxrad and minrad radii.
c It also rounds of the innner rim around minrad, using either:
c   + shaperad: a concave ( <= 0), vertical (1) or rounded (>= 1) rim 
c   + roundwidth/softedge 
c-----------------------------------------------------------------------
	subroutine DestroyDustR()
	use Parameters
	IMPLICIT NONE
	integer ii,k,i,j,ii2
	real*8 tot,f,T,r1,r2,theta
	real*8 scale,RoundOff ! Gijsexp, for rounding SDP

c	MassTot0=0d0
c	do i=0,D%nR-1
c	do j=1,D%nTheta-1
c		MassTot0=MassTot0+C(i,j)%dens*C(i,j)%V
c	enddo
c	enddo

	do j=1,D%nTheta-1

	call tellertje(j,D%nTheta-1)

	theta=D%theta_av(j)*180d0/pi

	do i=0,D%nR-1

	   !  C(i,j)%w(ii) goes from density fraction to partial density
	   do ii=1,ngrains	
	      C(i,j)%w(ii)=C(i,j)%w(ii)*C(i,j)%dens
	   enddo

	   !  Destroy dust using minrad#,maxrad#,softedge,shaperad# and roundwidth# keywords 
	   do ii=1,ngrains
	      r1=D%R_av(i)/AU*sin(D%thet(j))**Grain(ii)%shaperad
	      r2=D%R_av(i)/AU	!*sin(D%thet(j+1))**Grain(ii)%shaperad/AU

	      !  destroy outside [minrad,maxrad], round off using shaperad
	      if(r2.gt.Grain(ii)%maxrad.or.r1.lt.Grain(ii)%minrad
     &			.or.theta.gt.Grain(ii)%maxtheta.or.theta.lt.Grain(ii)%mintheta) then
		 C(i,j)%w(ii)=0d0

	      !  round rim at minrad?
	      else if (Grain(ii)%shaperad.eq.0d0.and.
     &                r1.lt.Grain(ii)%minrad+Grain(ii)%roundwidth) then

		 scale=RoundOff(r1,Grain(ii)%minrad+Grain(ii)%roundwidth,
     &                          Grain(ii)%roundtype,Grain(ii)%roundpow,1d-300)

                 ! scale from dens0 if no significant dust destruction
		 if(C(i,j)%w(ii).gt.(C(i,j)%w0(ii)*C(i,j)%dens0*scale)) then
		    C(i,j)%w(ii)=C(i,j)%w0(ii)*C(i,j)%dens0*scale
		 endif

	      !  don't destroy dust
	      else if(.not.tdes_iter) then
		 C(i,j)%w(ii)=C(i,j)%w0(ii)*C(i,j)%dens0
	      endif
	   enddo

	   !  rebuild density array from partial density
	   C(i,j)%dens=0d0
	   do ii=1,ngrains
	      C(i,j)%dens=C(i,j)%dens+C(i,j)%w(ii)
	   enddo
	   tot=C(i,j)%dens
c----------------------
	   ! Set C(i,j)%w(ii) back from partial density to density fraction
	   if(tot.gt.1d-50) then 
	      C(i,j)%w(1:ngrains)=C(i,j)%w(1:ngrains)/tot
	   else
	      C(i,j)%w(1:ngrains)=C(i,j)%w0(1:ngrains)
	      C(i,j)%dens=1d-50
	   endif
c----------------------
	   if(C(i,j)%dens.gt.C(i,j)%dens0) C(i,j)%dens=C(i,j)%dens0
	   if(C(i,j)%dens.lt.1d-50) then
	      C(i,j)%w(1:ngrains)=C(i,j)%w0(1:ngrains)
	      C(i,j)%dens=1d-50
	      C(i,j)%gasfrac=1d0
	   else
	      C(i,j)%gasfrac=1d0-C(i,j)%dens/C(i,j)%dens0
	   endif
	   C(i,j)%mass=C(i,j)%dens*C(i,j)%V
	enddo

	enddo

	call CheckCells()

	return
	end
	
	






	subroutine denstempFits(IMDIM,Rmax)
	use Parameters
	IMPLICIT NONE
	integer IMDIM
	real*8 imT(IMDIM,IMDIM),imD(IMDIM,IMDIM),Rmax,theta,x,y,R
	integer i,j,k,l

	imT=0d0
	imD=0d0

	do i=1,IMDIM
	call tellertje(i,IMDIM)
	do j=1,IMDIM
		x=Rmax*real(i-IMDIM/2)/real(IMDIM/2)
		y=Rmax*real(j-IMDIM/2)/real(IMDIM/2)
		R=sqrt(x**2+y**2)
		if(R.gt.D%R(D%nR).or.R.lt.D%R(1)) goto 1
		theta=acos(x/R)
		if(theta.gt.pi/2d0) theta=pi-theta
		do k=1,D%nR-1
			if(R.ge.D%R(k).and.R.lt.D%R(k+1)) goto 2
		enddo
2		do l=1,D%nTheta-1
			if(theta.ge.D%thet(l).and.theta.le.D%thet(l+1)) goto 3
		enddo
3		if(k.gt.D%nR-1) k=D%nR-1
		if(l.gt.D%nTheta-1) l=1
		if(l.lt.1) l=1
		imT(i,j)=C(k,l)%T
		imD(i,j)=C(k,l)%dens

1		continue
	enddo
	enddo


	call FITSImageStruct(imT,IMDIM,Rmax,'T')
	call FITSImageStruct(imD,IMDIM,Rmax,'D')

	return
	end



	subroutine FITSImageStruct(im,IMDIM,Rmax,add)
	use Parameters
	IMPLICIT NONE
	integer IMDIM
	character*500 filename,cmd
	character add
	real*8 flux,w1,w2,il1,i1,ir1,im(IMDIM,IMDIM)
	real*8 array(IMDIM,IMDIM),Rmax,cornerpix,pixscale
	integer i,k,i10,il10,il100,il1000,il10000
	integer ir10,ir100,ir1000,ir10000

      integer status,unit,blocksize,bitpix,naxis,naxes(2)
      integer j,group,fpixel,nelements
      logical simple,extend


	ir10000=2d0*Rmax/10000d0
	ir1000=(2d0*Rmax-10000d0*ir10000)/1000d0
	ir100=(2d0*Rmax-1000d0*ir1000-10000d0*ir10000)/100d0
	ir10=(2d0*Rmax-100d0*ir100-1000d0*ir1000-10000d0*ir10000)/10d0
	ir1=(2d0*Rmax-10d0*ir10-100d0*ir100-1000d0*ir1000-10000d0*ir10000)

	if(ir1.ge.9.95d0) then
		ir10=ir10+1
		ir1=ir1-10d0
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
	write(filename,'(a,"Image",a1,"_fov",i1,i1,i1,i1,f3.1,".fits")') outdir(1:len_trim(outdir))
     &			,add,ir10000,ir1000,ir100,ir10,ir1

      status=0
C     Get an unused Logical Unit Number to use to create the FITS file
      call ftgiou(unit,status)
C     create the new empty FITS file
      blocksize=1
      call ftinit(unit,filename,blocksize,status)
1	if(status.eq.105) then
		write(*,'("FITS file already exists, overwriting")')
		write(9,'("FITS file already exists, overwriting")')
		write(cmd,'("rm ",a)') filename(1:len_trim(filename))
		call system(cmd)
      call ftclos(unit, status)
		status=0
      call ftfiou(unit, status)
		status=0
      call ftgiou(unit,status)
		status=0
      call ftinit(unit,filename,blocksize,status)
      goto 1
	endif	

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


	subroutine RegridThetaOld(n1)
	use Parameters
	IMPLICIT NONE
	integer n,n1,i,j,ii,j0,iopac
	real*8 Kext,tau,theta0,Vold
	character*500 thetagridfile

	n=n1
	if(n.gt.(D%nTheta/2)) n=D%nTheta/2
	
	j0=0
	
	do j=D%nTheta-1,1,-1
	tau=0
	do i=1,D%nR-1
		Kext=0d0
		do ii=1,ngrains
		   do iopac=1,Grain(ii)%nopac
			Kext=Kext+C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)*Grain(ii)%Kpstar(iopac)
		   enddo
		enddo
		tau=tau+(D%R(i+1)-D%R(i))*C(i,j)%dens*Kext*AU
	enddo
	if(tau.gt.0.25d0) j0=j-1
	enddo

	if(j0.eq.0) return
	
1	continue

	theta0=D%theta_av(j)
	if((theta0/real(n)).lt.(pi/2d0-theta0)/real(D%nTheta-n)) then
		n=0
		theta0=(pi/2d0)/real(D%nTheta)
	endif
	do j=1,n
c		D%Theta(j)=1d0-cos(theta0)*real(j-1)/real(n-1)
		D%Theta(j)=cos(theta0*real(j-1)/real(n-1))
	enddo
	do j=n+1,D%nTheta
c		D%Theta(j)=theta0+(pi/2d0-theta0)*real(j-n-1)/real(D%nTheta-n-1)
c		D%Theta(j)=cos(D%Theta(j))
		D%Theta(j)=cos(theta0)-cos(theta0)*(real(j-n-1)/real(D%nTheta-n-1))**0.5
	enddo
	call sort(D%Theta(1:D%nTheta),D%nTheta)
	do j=1,D%nTheta-1
		D%theta_av(j)=acos((D%Theta(D%nTheta-j+1)+D%Theta(D%nTheta-j))/2d0)
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

	write(thetagridfile,'(a,"thetagrid.dat")') outdir(1:len_trim(outdir))
	open(unit=60,file=thetagridfile,RECL=100)
	do i=1,D%nTheta
		write(60,*) acos(D%Theta(i))
	enddo
	close(unit=60)

	
	do j=1,D%nTheta-1
		do i=0,D%nR-1
			C(i,j)%xedge(3)=D%Theta(j)
			C(i,j)%xedge(4)=D%Theta(j+1)
			C(i,j)%xedge2(3)=D%Theta(j)**2
			C(i,j)%xedge2(4)=D%Theta(j+1)**2
			Vold=C(i,j)%V
			C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
			C(i,j)%dens0=C(i,j)%dens0*Vold/C(i,j)%V
			C(i,j)%dens=C(i,j)%dens*Vold/C(i,j)%V
			C(i,j)%gasdens=C(i,j)%gasdens*Vold/C(i,j)%V
			call CheckMinimumDensity(i,j)
		enddo
	enddo
	
	return
	end




	subroutine RegridTheta(n1)
	use Parameters
	IMPLICIT NONE
	integer n,n1,i,j,ii,j0,iopac,jj,iter,nii
	real*8 Kext,tau,theta0,Vold
	character*500 thetagridfile
	
	real*8 d0,d1,z(D%nTheta),dens,t_sh,sh
	real*8,allocatable :: sh0(:),sh_tmp(:)

	allocate(sh0(ngrains*(D%nR-1)))
	allocate(sh_tmp(ngrains*(D%nR-1)))

	n=n1

	do i=1,D%nR-1
		do ii=1,ngrains
			d0=C(i,D%nTheta-1)%dens*C(i,D%nTheta-1)%w(ii)
			d1=d0
			d0=d0/exp(0.5d0) ! p(z)=p0*e^{-z^2/2H_p^2} -> p(H_p)=p0*e^-0.5
			if(d0.gt.1d-50) then
				z(D%nTheta-1)=D%R_av(i)*cos(D%theta_av(D%nTheta-1))/AU
				do j=D%nTheta-2,1,-1
					z(j)=D%R_av(i)*cos(D%theta_av(j))/AU
					dens=C(i,j)%dens*C(i,j)%w(ii)
					if(dens.lt.d0) then
						sh=(z(j)+(z(j+1)-z(j))*(d0-dens)/(d0-d1))
						goto 2
					else if(dens.le.1d-50) then
						sh=z(j)
						goto 2
					endif
					d1=dens
				enddo
			else
				z(D%nTheta/2)=D%R_av(i)*cos(D%theta_av(D%nTheta/2))/AU
				sh=z(D%nTheta/2)
			endif
2			continue
			sh_tmp((i-1)*ngrains+ii)=sh/(D%R_av(i)/AU)
		enddo
	enddo

	sh_tmp=abs(sh_tmp)
	call sort(sh_tmp,(D%nR-1)*ngrains)
	
	nii=ngrains*(D%nR-1)

	sh0(1)=sh_tmp(1)
	sh_tmp(1)=-1d0
	sh0(2)=sh_tmp(nii)
	sh_tmp(nii)=-1d0
	sh0(3)=sh_tmp(nii/2)
	sh_tmp(nii/2)=-1d0
	sh0(4)=sh_tmp(nii/4)
	sh_tmp(D%nR/4)=-1d0
	sh0(5)=sh_tmp(3*nii/4)
	sh_tmp(3*nii/4)=-1d0

	ii=5
	do i=1,nii
		if(sh_tmp(i).gt.0d0) then
			ii=ii+1
			sh0(ii)=sh_tmp(i)
		endif
	enddo

	nii=n
	
	do iter=1,2

	if(n.gt.(D%nTheta/2)) n=D%nTheta/2
	
	j0=0
	
	do j=D%nTheta-1,1,-1
	tau=0
	do i=1,D%nR-1
		Kext=0d0
		do ii=1,ngrains
		   do iopac=1,Grain(ii)%nopac
			Kext=Kext+C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)*Grain(ii)%Kpstar(iopac)
		   enddo
		enddo
		tau=tau+(D%R(i+1)-D%R(i))*C(i,j)%dens*Kext*AU
	enddo
	if(tau.lt.0.1d0) goto 4
	enddo

4	j0=j-1

	if(j0.le.0.or.j0.ge.D%nTheta-1) then
		deallocate(sh0)
		deallocate(sh_tmp)
		return
	endif
	
1	continue

	theta0=D%theta_av(j0)/1.25d0
	if((theta0/real(n)).lt.(pi/2d0-theta0)/real(D%nTheta-n)) then
		n=0
		theta0=(pi/2d0)/real(D%nTheta)
	endif
	do j=1,n
c		D%Theta(j)=1d0-cos(theta0)*real(j-1)/real(n-1)
		D%Theta(j)=cos(theta0*real(j-1)/real(n-1))
	enddo
	do j=n+1,D%nTheta-nii
c		D%Theta(j)=theta0+(pi/2d0-theta0)*real(j-n)/real(D%nTheta-n)
c		D%Theta(j)=cos(D%Theta(j))
		D%Theta(j)=cos(theta0)-cos(theta0)*(real(j-n)/real(D%nTheta-n-nii))
	enddo


	ii=D%nTheta-nii

	call sort(D%Theta(1:ii),ii)

	nii=0
	do i=1,ngrains*(D%nR-1)
		sh=sh0(i)*0.1
		j0=1
3		t_sh=atan(1d0/sh)
		t_sh=cos(t_sh)
		do j=j0,ii
			if(sh.gt.8d0*sh0(i)) exit
			if(D%Theta(j).gt.t_sh) then
				if(j.eq.j0) then
					ii=ii+1
					nii=nii+1
					if(ii.le.D%nTheta) then
						D%Theta(ii)=t_sh
						call sort(D%Theta(1:ii),ii)
					endif
				endif
				sh=sh*2d0
				j0=j
				if(ii.lt.D%nTheta) goto 3
			endif
		enddo				
	enddo
	
	if(nii.gt.D%nTheta-n) nii=D%nTheta-n
	
	enddo

	do j=ii+1,D%nTheta
		D%Theta(j)=cos(theta0*(real(j-ii)/real(D%nTheta-ii+1)))
	enddo

	call sort(D%Theta(1:D%nTheta),D%nTheta)
	do j=1,D%nTheta-1
		D%theta_av(j)=acos((D%Theta(D%nTheta-j+1)+D%Theta(D%nTheta-j))/2d0)
	enddo

5	continue
	do j=1,D%nTheta-2
		if(D%theta_av(j).eq.D%theta_av(j+1)) then
			if(j.eq.1) then
				D%theta_av(j)=D%theta_av(j)/2d0
			else
				D%theta_av(j)=(D%theta_av(j)+D%theta_av(j-1))/2d0
			endif
			goto 5
		endif
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
	
	do j=1,D%nTheta-1
		do i=0,D%nR-1
			C(i,j)%xedge(3)=D%Theta(j)
			C(i,j)%xedge(4)=D%Theta(j+1)
			C(i,j)%xedge2(3)=D%Theta(j)**2
			C(i,j)%xedge2(4)=D%Theta(j+1)**2
			Vold=C(i,j)%V
			C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
			C(i,j)%dens0=C(i,j)%dens0*Vold/C(i,j)%V
			C(i,j)%dens=C(i,j)%dens*Vold/C(i,j)%V
			C(i,j)%gasdens=C(i,j)%gasdens*Vold/C(i,j)%V
			C(i,j)%mass=C(i,j)%dens*C(i,j)%V
			call CheckMinimumDensity(i,j)
		enddo
	enddo

	deallocate(sh0)
	deallocate(sh_tmp)
	
	return
	end
	

	real*8 function determinegasfrac(T,i,j,ii)
	use Parameters
	IMPLICIT NONE
	real*8 A,B,T,dens
	integer i,j,ii,k
	real*8 mu,Rgas,mole,dsdt,Pv,G,Omega
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
	parameter(Rgas=8.314e7) ! erg/mol/K
	parameter(mole=6.022e23)
	parameter(G=6.67300d-8) ! in cm^3/g/s

	if(Grain(ii)%material.eq.'UNKNOWN') then
		dens=1d0
	else
		dens=0d0
		do k=1,ngrains
			if(Grain(ii)%material.eq.Grain(k)%material) dens=dens+C(i,j)%w0(k)
		enddo
	endif
	dens=dens*C(i,j)%dens0

	A=Grain(ii)%TdesA
	B=Grain(ii)%TdesB

	if(abs(B).gt.1d-15) then
		if(T.ge.3d0) then
			determinegasfrac=10d0**(A/B-1d4/(T*B)-log10(T))/dens
		else
			determinegasfrac=0d0
		endif
	else
		if(T.gt.A) then
			determinegasfrac=2d0
		else
			determinegasfrac=0d0
		endif
	endif

c ------------ for the IN05 sublimation law --------------------
c	determinegasfrac=(T/2000d0)**(1d0/1.95d-2)/(100d0*dens)
c ------------ for the IN05 sublimation law --------------------

c------------- determine the time for evaporation --------------
	Pv=10d0**(A/B-1d4/(T*B)-log10(T))*Rgas*T/(mu*mole)
	dsdt=(Pv/Grain(ii)%rho(1))*sqrt(mu/(2d0*pi*kb*T))

	Omega=2d0*pi*sqrt(D%R_av(i)**3/(G*D%Mstar))

	if((Grain(ii)%rv/dsdt).gt.Omega) determinegasfrac=0d0
	
	return
	end
	

	subroutine destroyQHP(TqhpMax)
	use Parameters
	IMPLICIT NONE
	integer i,j,ii,itemp,iopac
	real*8 TqhpMax,MassTot,MassHigh,tot,tot0
	character*500 output
	
	if(TqhpMax.lt.0d0) return

	do i=0,D%nR-1
		do ii=1,ngrains
		if(Grain(ii)%qhp) then
		if(D%R(i+1).le.Grain(ii)%maxrad.and.D%R(i).ge.Grain(ii)%minrad) then
			MassTot=0d0
			MassHigh=0d0
			do j=1,D%nTheta-1
				MassTot=MassTot+C(i,j)%dens0*C(i,j)%w0(ii)*C(i,j)%V
				do itemp=1,NTQHP
					if(tgridqhp(itemp).gt.TqhpMax) then
						MassHigh=MassHigh+C(i,j)%tdistr(Grain(ii)%qhpnr,itemp)*C(i,j)%dens0*C(i,j)%w0(ii)*C(i,j)%V
					endif
				enddo
			enddo
			if((MassHigh/MassTot).gt.1d-8) then
				do j=1,D%nTheta-1
					C(i,j)%w(ii)=0d0
				enddo
			else
				do j=1,D%nTheta-1
					C(i,j)%w(ii)=C(i,j)%w0(ii)
				enddo
			endif
		endif
		endif
		enddo
		do j=1,D%nTheta-1
			tot=0d0
			do ii=1,ngrains
				tot=tot+C(i,j)%w(ii)
			enddo
			if(tot.gt.1d-100) then
				C(i,j)%w(1:ngrains)=C(i,j)%w(1:ngrains)/tot
			endif
			if(tot0.lt.1d-100) tot0=1d0
			C(i,j)%dens=C(i,j)%dens*tot/tot0
			if(C(i,j)%dens.gt.C(i,j)%dens0) C(i,j)%dens=C(i,j)%dens0
			if(C(i,j)%dens.lt.1d-50) then
				C(i,j)%w(1:ngrains)=C(i,j)%w0(1:ngrains)
				C(i,j)%dens=1d-50
				C(i,j)%gasfrac=1d0
			else
				C(i,j)%gasfrac=1d0-C(i,j)%dens/C(i,j)%dens0
			endif
		enddo
	enddo

	write(output,'(a,"composition.dat")') outdir(1:len_trim(outdir))
	open(unit=20,file=output,RECL=6000)
	do i=1,D%nR-1
		write(20,*) D%R_av(i)/AU,(C(i,D%nTheta-1)%w(ii),ii=1,ngrains)
	enddo
	close(unit=20)
	
	return
	end
	
	
	subroutine determineRdes()
	use Parameters
	use Diskstruct
	IMPLICIT NONE
	integer ii,Rthin,Rbw,iter,iopac
	real*8 Tevap,Tthin,Tbw,W,r,spec(nlam)
	real*8 determineT,determineTP,fwmin,tau,Kext
	real*8 f,determinegasfrac,dR,dRprev,dRnew
	type(photon) phot
	logical BBGrains

	if(.not.allocated(Rthindes)) then
		allocate(Rthindes(D%nTheta))
		allocate(Rbwdes(D%nTheta))
	endif

	fwmin=1d-5
	
	do j=1,D%nTheta-1
	Rthin=D%nR-1
	Rbw=D%nR-1
	do i=D%nR-1,0,-1
		if(i.ne.0) then
			r=D%R(i)*AU
		else
			r=D%R(1)*AU
		endif
		W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))
		phot%i=i
		phot%j=j
		if(tcontact) then
			spec(1:nlam)=0d0
			BBGrains=.false.
			do ii=1,ngrains
				if(.not.Grain(ii)%qhp.and.C(i,j)%w(ii).gt.0d0) then
				   do iopac=1,Grain(ii)%nopac
				      spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)*C(i,j)%wopac(ii,iopac)
     &			*C(i,j)%w(ii)/(pi*D%Rstar**2)
				   enddo
				   BBGrains=.true.
				endif
			enddo
			if(BBGrains) then
				call integrate(spec,phot%E)
				phot%E=W*phot%E
				Tthin=determineT(phot)
				phot%E=phot%E*fBW(1,j)
				Tbw=determineT(phot)
			else
				Tthin=D%Tstar*W**(0.25)
				Tbw=D%Tstar*(W*fBW(1,j))**(0.25)
			endif
			do ii=1,ngrains
				f=determinegasfrac(Tthin,i,j,ii)
				if(f.le.1d0.or.Tthin.lt.50d0) then
					Rthin=i
				endif
				f=determinegasfrac(Tbw,i,j,ii)
				if(f.le.1d0.or.Tbw.lt.50d0) then
					Rbw=i
				endif
			enddo
		else
			do ii=1,ngrains
			if(.not.Grain(ii)%qhp.and.C(i,j)%w0(ii).gt.0d0) then
			        spec(1:nlam)=0d0
				do iopac=1,Grain(ii)%nopac
				   spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)*C(i,j)%wopac(ii,iopac)/(pi*D%Rstar**2)
				enddo
				call integrate(spec,phot%E)
				phot%E=W*phot%E
				Tthin=determineTP(phot,ii)
				phot%E=phot%E*fBW(ii,j)
				Tbw=determineTP(phot,ii)
				f=determinegasfrac(Tthin,i,j,ii)
				if(f.le.1d0.or.Tthin.lt.50d0) then
					Rthin=i
				endif
				f=determinegasfrac(Tbw,i,j,ii)
				if(f.le.1d0.or.Tbw.lt.50d0) then
					Rbw=i
				endif
			endif
			enddo
		endif
	enddo

	if(Rthin.ge.D%nR-1) Rthin=1
	if(Rbw.ge.D%nR-1) Rbw=1

	if(Rthin.gt.1) then
	dR=0.5
	dRprev=0d0
	do iter=1,100
		r=((1d0-dR)*D%R(Rthin-1)+dR*D%R(Rthin))*AU
		W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))
		phot%i=Rthin-1
		phot%j=j
		if(tcontact) then
			spec(1:nlam)=0d0
			BBGrains=.false.
			do ii=1,ngrains
				if(.not.Grain(ii)%qhp.and.C(Rthin-1,j)%w(ii).gt.0d0) then
				   do iopac=1,Grain(ii)%nopac
				      spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)*
     &            C(Rthin-1,j)%wopac(ii,iopac)*C(Rthin-1,j)%w(ii)/(pi*D%Rstar**2)
				   enddo
				   BBGrains=.true.
				endif
			enddo
			if(BBGrains) then
				call integrate(spec,phot%E)
				phot%E=W*phot%E
				Tthin=determineT(phot)
			else
				Tthin=D%Tstar*W**(0.25)
			endif
			do ii=1,ngrains
				f=determinegasfrac(Tthin,Rthin-1,j,ii)
				if(f.le.1d0) then
					if(dRprev.lt.dR) then
						dRnew=(dR+dRprev)/2d0
					else
						dRnew=dR/2d0
					endif
					goto 2
				endif
			enddo
			if(dRprev.gt.dR) then
				dRnew=(dR+dRprev)/2d0
			else
				dRnew=(1d0+dR)/2d0
			endif
		else
			do ii=1,ngrains
			if(.not.Grain(ii)%qhp.and.C(Rthin-1,j)%w0(ii).gt.0d0) then
			        spec(1:nlam)=0d0
				do iopac=1,Grain(ii)%nopac
				   spec(1:nlam)=spec(1:nlam)+D%Fstar(1:nlam)*Grain(ii)%Kabs(iopac,1:nlam)*C(Rthin-1,j)%wopac(ii,iopac)/(pi*D%Rstar**2)
				enddo
				call integrate(spec,phot%E)
				phot%E=W*phot%E
				Tthin=determineTP(phot,ii)
				f=determinegasfrac(Tthin,Rthin-1,j,ii)
				if(f.le.1d0) then
					if(dRprev.lt.dR) then
						dRnew=(dR+dRprev)/2d0
					else
						dRnew=dR/2d0
					endif
					goto 2
				endif
				if(dRprev.gt.dR) then
					dRnew=(dR+dRprev)/2d0
				else
					dRnew=(1d0+dR)/2d0
				endif
			endif
			enddo
		endif
2		continue
		dRprev=dR
		dR=dRnew
	enddo
	Rthindes(j)=((1d0-dR)*D%R(Rthin-1)+dR*D%R(Rthin))
	else
	Rthindes(j)=D%Rin
	endif


	if(Rbw.gt.1) then
	dR=0.5
	dRprev=0d0
	do iter=1,100
		r=((1d0-dR)*D%R(Rbw-1)+dR*D%R(Rbw))*AU
		W=0.5d0*(1d0-sqrt(1d0-(D%Rstar/r)**2))
		phot%i=Rbw-1
		phot%j=j
		if(tcontact) then
			spec(1:nlam)=0d0
			BBGrains=.false.
			phot%E=0d0
			do ii=1,ngrains
				if(.not.Grain(ii)%qhp.and.C(Rbw-1,j)%w(ii).gt.0d0) then
				   do iopac=1,Grain(ii)%nopac
						phot%E=phot%E+Grain(ii)%Kpabsstar(iopac)*D%Lstar*
     &       C(Rbw-1,j)%wopac(ii,iopac)*C(Rbw-1,j)%w(ii)/(pi*D%Rstar**2)
				   enddo
				   BBGrains=.true.
				endif
			enddo
			if(BBGrains) then
				phot%E=phot%E*fBW(1,j)
				Tbw=determineT(phot)
			else
				Tbw=D%Tstar*(W*fBW(1,j))**(0.25)
			endif
			do ii=1,ngrains
				f=determinegasfrac(Tbw,Rbw-1,j,ii)
				if(f.le.1d0) then
					if(dRprev.lt.dR) then
						dRnew=(dR+dRprev)/2d0
					else
						dRnew=dR/2d0
					endif
					goto 3
				endif
			enddo
			if(dRprev.gt.dR) then
				dRnew=(dR+dRprev)/2d0
			else
				dRnew=(1d0+dR)/2d0
			endif
		else
			do ii=1,ngrains
			if(.not.Grain(ii)%qhp.and.C(Rbw-1,j)%w0(ii).gt.0d0) then
				phot%E=0d0
		        do iopac=1,Grain(ii)%nopac
					phot%E=phot%E+Grain(ii)%Kpabsstar(iopac)*D%Lstar*C(Rbw-1,j)%wopac(ii,iopac)/(pi*D%Rstar**2)
				enddo
				phot%E=W*phot%E
				Tbw=determineTP(phot,ii)
				f=determinegasfrac(Tbw,Rbw-1,j,ii)
				if(f.le.1d0) then
					if(dRprev.lt.dR) then
						dRnew=(dR+dRprev)/2d0
					else
						dRnew=dR/2d0
					endif
					goto 3
				endif
				if(dRprev.gt.dR) then
					dRnew=(dR+dRprev)/2d0
				else
					dRnew=(1d0+dR)/2d0
				endif
			endif
			enddo
		endif
3		continue
		dRprev=dR
		dR=dRnew
	enddo
	Rbwdes(j)=((1d0-dR)*D%R(Rbw-1)+dR*D%R(Rbw))
	else
	Rbwdes(j)=D%Rin
	endif

	enddo

	write(*,'("Optically thin desctruction:    ",f10.3," AU")') Rthindes(D%nTheta-1)
	write(*,'("Backwarming desctruction:       ",f10.3," AU")') Rbwdes(D%nTheta-1)
	write(9,'("Optically thin desctruction:    ",f10.3," AU")') Rthindes(D%nTheta-1)
	write(9,'("Backwarming desctruction:       ",f10.3," AU")') Rbwdes(D%nTheta-1)

	return
	end
	
	subroutine CheckCells()
	use Parameters
	IMPLICIT NONE
	integer i,j,ii
	real*8 massmax(ngrains),mass(ngrains)
	logical rescale

	do i=0,D%nR-1
	do j=1,D%nTheta-1
		rescale=.false.
		do ii=1,ngrains
			mass(ii)=C(i,j)%dens*C(i,j)%w(ii)
			massmax(ii)=C(i,j)%dens0*C(i,j)%w0(ii)
			if(mass(ii).gt.massmax(ii)) then
				mass(ii)=massmax(ii)
				rescale=.true.
			endif
		enddo
		if(rescale) then
			C(i,j)%dens=0d0
			do ii=1,ngrains
				C(i,j)%dens=C(i,j)%dens+mass(ii)
			enddo
			if(C(i,j)%dens.gt.1d-50) then
				do ii=1,ngrains
					C(i,j)%w(ii)=mass(ii)/C(i,j)%dens
				enddo
				C(i,j)%gasfrac=1d0-C(i,j)%dens/C(i,j)%dens0
			else
				C(i,j)%dens=1d-50
				C(i,j)%gasfrac=1d0
				C(i,j)%w(1:ngrains)=C(i,j)%w0(1:ngrains)
			endif
			C(i,j)%mass=C(i,j)%dens*C(i,j)%V
		endif
	enddo
	enddo

	return
	end
	


	subroutine TraceStellar(phot,tau0,escape,hitstar,locfield)
	use Parameters
	IMPLICIT NONE
	type(Photon) phot
	real*8 tau0,tau,tau1,v,Kext,ct,EJv,Kabs,MinDist,dmin,wlam
	real*8 KabsBBGrains,nRWinteract,vx,vy,vz,inp
	logical escape,hitstar,locfield,RandomWalk
	integer inext,jnext,ntrace,j,ii,iopac
c	character*500 s

	escape=.false.
	hitstar=.false.
	tau=0d0
	ntrace=0

1	continue

	Kext=0d0
	do ii=1,ngrains
	   do iopac=1,Grain(ii)%nopac
		Kext=Kext+Grain(ii)%Kpstar(iopac)*C(phot%i,phot%j)%wopac(ii,iopac)*C(phot%i,phot%j)%w(ii)
	   enddo
	enddo

	call Trace2edge(phot,v,inext,jnext)

	tau1=Kext*v*C(phot%i,phot%j)%dens*AU
	if((tau+tau1).gt.tau0) then
		v=v*(tau0-tau)/tau1
		tau1=tau0-tau
		phot%x=phot%x+phot%vx*v
		phot%y=phot%y+phot%vy*v
		phot%z=phot%z+phot%vz*v
		tautot=tautot+(tau0-tau)
		phot%onEdge=.false.
		return
	endif

	phot%x=phot%x+phot%vx*v
	phot%y=phot%y+phot%vy*v
	phot%z=phot%z+phot%vz*v

	tau=tau+tau1

	if(inext.ge.D%nR) then
		escape=.true.
		return
	endif
	if(inext.lt.0) then
		escape=.true.
		hitstar=.true.
		return
	endif

	phot%i=inext
	phot%j=jnext

	ntrace=ntrace+1
	if(ntrace.gt.D%nR*D%nTheta*2) then
		print*,'raar!',phot%i,phot%j,phot%x,phot%y,phot%z
		return
	endif

	goto 1

	return
	end

c-----------------------------------------------------------------------
c  Compute the temperature in the midplane according to a
c  Shakura-Sunyaev disk (alpha-disk)
c  This is a test routine, it doesn't affect disk temperature
c-----------------------------------------------------------------------

	subroutine ShakuraSunyaev()
	use Parameters
	IMPLICIT NONE
	integer i,j,iopac
	real*8 T,Surfdens,WeightedAlpha,mu,G,tau,T2,T3
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
	parameter(G=6.67300d-8) ! in cm^3/g/s^2
	character*500 file
	type(photon) phot
	real*8 v,KR,Kabs(nlam),spec(nlam),dBB(nlam),int1,int2,Kext(1:nlam),Kp
	integer inext,jnext,iT,iter,ii

	return
	write(file,'(a,"ShakuraSunyaev.dat")') outdir(1:len_trim(outdir))
	open(unit=42,file=file,RECL=6000)
	
	do i=1,D%nR-1
		iT=int(C(i,D%nTheta-1)%T/dT)+1
		if(iT.lt.1) iT=1
		if(iT.gt.TMAX-1) iT=TMAX-1

		do iter=1,10

		j=D%nTheta-1
		tau=0d0
		phot%x=(D%R_av(i)/AU)*sin(D%theta_av(j))*sin(0.2d0)
		phot%y=(D%R_av(i)/AU)*sin(D%theta_av(j))*cos(0.2d0)
		phot%z=(D%R_av(i)/AU)*cos(D%theta_av(j))
		phot%vx=sin(1d-5)*cos(0.1d0)
		phot%vy=sin(1d-5)*sin(0.1d0)
		phot%vz=cos(1d-5)
		phot%i=i
		phot%j=D%nTheta-1
		phot%onEdge=.false.
		tau=0d0

		Kext=0d0
		do ii=1,ngrains
		   do iopac=1,Grain(ii)%nopac
			Kext=Kext+C(phot%i,phot%j)%w(ii)*C(phot%i,phot%j)%wopac(ii,iopac)*Grain(ii)%Kext(iopac,1:nlam)
		   enddo
		enddo
		dBB(1:nlam)=BB(1:nlam,iT+1)-BB(1:nlam,iT)
		call integrate(dBB/Kext,int1)
		call integrate(dBB,int2)

		KR=int2/int1

1		continue

		call Trace2edge(phot,v,inext,jnext)
		phot%onEdge=.true.

		phot%x=phot%x+phot%vx*v
		phot%y=phot%y+phot%vy*v
		phot%z=phot%z+phot%vz*v

		tau=tau+v*C(phot%i,phot%j)%dens*AU*KR

		if(inext.ge.D%nR.or.inext.lt.0) goto 2

		phot%i=inext
		phot%j=jnext

		goto 1
2		continue

		!  Calculate surface density and density-weighted viscosity (in case of deadzones)
		Surfdens=0d0
		WeightedAlpha=0d0
		do j=1,D%nTheta-1
			Surfdens=Surfdens+C(i,j)%gasdens*gas2dust*C(i,j)%V
			WeightedAlpha=WeightedAlpha+C(i,j)%gasdens*gas2dust*C(i,j)%V*C(i,j)%alphavis
		enddo
		WeightedAlpha=WeightedAlpha/Surfdens
c		print*,i,D%R_av(i),WeightedAlpha
		Surfdens=Surfdens/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)
		

		T=( tau*27d0*kb*WeightedAlpha*Surfdens*sqrt(G*D%Mstar/D%R_av(i)**3)/(32d0*mu*sigma) )**(1d0/3d0)
		
		T2=(3d0*tau*D%Mdot*(G*D%Mstar/D%R_av(i)**3)/(32d0*pi*sigma))**0.25d0
		
		KR=2d0*tau/Surfdens

		T3=(D%Mdot**2*KR*mu/(pi**2*64d0*sigma*WeightedAlpha*kb))**(1d0/5d0) * (G*D%Mstar/D%R_av(i)**3)**(3d0/10d0)
		
		iT=int(T/dT)+1
		if(iT.lt.1) iT=1
		if(iT.gt.TMAX-1) iT=TMAX-1

		enddo

		write(42,*) D%R_av(i)/AU,T,C(i,D%nTheta-1)%T,T2,T3,tau,Surfdens
		
	enddo

c	print*,(G*D%Mstar)**(1d0/3d0)*(mu/(pi**2*64d0*sigma*kb))**(10d0/45d0)
c     &			/AU*(Msun/(365d0*24d0*3600d0))**(20d0/45d0)/(160d0**(10d0/9d0))
	
	close(unit=42)
	
	return
	end
	
	
	
	real*8 function ShakuraSunyaevIJ(i,j)
	use Parameters
	IMPLICIT NONE
	integer i,j
	real*8 T,Surfdens,mu,G,tau
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
	parameter(G=6.67300d-8) ! in cm^3/g/s
	character*500 file
	type(photon) phot
	real*8 v,KR,Kabs(nlam),spec(nlam),dBB(nlam),int1,int2
	integer inext,jnext,iT,iter,ii,jj,iopac

	tau=0d0
	phot%x=(D%R_av(i)/AU)*sin(D%theta_av(j))*sin(0.2d0)
	phot%y=(D%R_av(i)/AU)*sin(D%theta_av(j))*cos(0.2d0)
	phot%z=(D%R_av(i)/AU)*cos(D%theta_av(j))
	phot%vx=sin(1d-5)*cos(0.1d0)
	phot%vy=sin(1d-5)*sin(0.1d0)
	phot%vz=cos(1d-5)
	
	phot%i=i
	phot%j=j
	phot%onEdge=.false.
	tau=0d0

1	continue

	call Trace2edge(phot,v,inext,jnext)
	phot%onEdge=.true.

	phot%x=phot%x+phot%vx*v
	phot%y=phot%y+phot%vy*v
	phot%z=phot%z+phot%vz*v

	tau=tau+v*C(phot%i,phot%j)%dens*AU

	if(inext.ge.D%nR.or.inext.lt.0) goto 2

	phot%i=inext
	phot%j=jnext
	
	goto 1
2	continue

	Surfdens=0d0
	do jj=1,j
		Surfdens=Surfdens+C(i,jj)%gasdens*gas2dust*C(i,jj)%V
	enddo
	Surfdens=Surfdens/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)

	iT=int(C(i,j)%T/dT)+1

	do iter=1,5

	if(iT.lt.1) iT=1
	if(iT.gt.TMAX-1) iT=TMAX-1

	dBB(1:nlam)=BB(1:nlam,iT+1)-BB(1:nlam,iT)
	call integrate(dBB,int1)
	Kabs(1:nlam)=0d0
	do ii=1,ngrains
	   do iopac=1,Grain(ii)%nopac
	      Kabs(1:nlam)=Kabs(1:nlam)+Grain(ii)%Kabs(iopac,1:nlam)*C(i,j)%wopac(ii,iopac)*C(i,j)%w(ii)
	   enddo
	enddo
	spec(1:nlam)=dBB(1:nlam)/Kabs(1:nlam)
	call integrate(spec,int2)
	KR=int1/int2


	T=(3d0*tau*KR*D%Mdot*(G*D%Mstar/D%R_av(i)**3)/(32d0*pi*sigma))**0.25d0

	iT=int(T/dT)+1

	enddo

	ShakuraSunyaevIJ=T
	
	return
	end
	


	subroutine CheckMinimumDensity(i,j)
	use Parameters
	IMPLICIT NONE
	integer i,j,ii
	real*8 tot(ngrains)

c Set minimum density of condensable stuff
	if(C(i,j)%dens0.le.1d-50) then
		C(i,j)%dens0=1d-50
	endif
c Set minimum dust density
	if(C(i,j)%dens.le.1d-50) then
		C(i,j)%dens=1d-50
	endif
c Set minimum gas density
	if(C(i,j)%gasdens.le.1d-50) then
		C(i,j)%gasdens=1d-50
	endif
c Make sure density of condensable stuff is larger than dust density
	do ii=1,ngrains
		tot(ii)=C(i,j)%dens*C(i,j)%w(ii)
	enddo
	do ii=1,ngrains
		if(tot(ii).lt.1d-50) then
			tot(ii)=1d-50
		endif
		if(tot(ii).gt.C(i,j)%dens0*C(i,j)%w0(ii)) then
			tot(ii)=C(i,j)%dens0*C(i,j)%w0(ii)
		endif
	enddo
	C(i,j)%dens=sum(tot(1:ngrains))
c Set minimum dust density
	if(C(i,j)%dens.le.1d-50) then
		C(i,j)%dens=1d-50
		do ii=1,ngrains
		   C(i,j)%w(ii)=C(i,j)%w0(ii)
		enddo
	else
		do ii=1,ngrains
		   C(i,j)%w(ii)=tot(ii)/C(i,j)%dens
		enddo
	endif
	C(i,j)%gasfrac=1d0-C(i,j)%dens/C(i,j)%dens0
	C(i,j)%mass=C(i,j)%dens*C(i,j)%V
	
	return
	end
	

	subroutine CheckNaN(string)
	use Parameters
	IMPLICIT NONE
	integer i,j,ii
	character*10 string
	logical cnan
	cnan=.false.

	do i=1,D%nR-1
	do j=1,D%nTheta-1
		if(IsNaN(C(i,j)%dens)) then
			cnan=.true.
			write(*,'(a,a,i0,i0)') trim(string),' dens ',i,j
		endif
		if(IsNaN(C(i,j)%dens0)) then
			cnan=.true.
			write(*,'(a,a,i0,i0)') trim(string),' dens0',i,j
		endif
c		do ii=1,ngrains
c			if(IsNaN(C(i,j)%w(ii))) write(*,'(a,a,i4,i,i)') trim(string),' w',ii,i,j
c		enddo
	enddo
	enddo
	
	if(cnan) stop 'I found a NaN'
	
	return
	end
	
	
