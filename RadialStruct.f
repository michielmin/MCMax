c-----------------------------------------------------------------------
c This routine calculates the radial structure of the disk, for a given
c viscous alpha. The surface density can be capped off at Toomre=2 with
c gravstable=.true.
c
c If getalpha=.true., this radial structure is ignored, and only used to
c calculate what alpha would be necessary to obtain the current density
c structure (radialalpha.dat)	
c
c Note: the exponential tail for MdotR does not work because
c DiskStructure() is called twice in a row without a call to
c MakeDeadZone() in between
c
c-----------------------------------------------------------------------

	subroutine RadialStruct()
	use Parameters
	use DiskStruct
	IMPLICIT NONE
	real*8 dens(D%nR,D%nTheta),Eint,Er(1:D%nR),mu,G,DiskMass,alpha,Mdot,Eint0
	real*8 Sig(1:D%nR),surfscale(1:D%nR)
	real*8 Sig_old(1:D%nR),Sig_lin(1:D%nR),Sig_dead(1:D%nR),Sig_used(1:D%nR) ! check only
	real*8 Q,Q0,fac,Eact,Sigact,Edead,Sigdead
	real*8 ToomreQ ! function
	real*8 prevmass(D%nR-1,D%nTheta-1),prevgasdens(D%nR-1,D%nTheta-1)
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
	parameter(G=6.67300d-8) ! in cm^3/g/s^2
	real*8 infall_rho,infall_mu,infall_mu0,infall_tot
	integer ii,imin,imax
	character*100 file

	logical testall,testscale
	integer itest
	character*10 type

	testall=.false.
	testscale=.false.
	itest=77

	! Loop to calculate surface density scaling
	do i=1,D%nR-1
	   if(gravstable.or.D%Rexp.lt.1d10) then
c	      Mdot=D%MdotR(i)
	      Mdot=D%Mdot  ! bug?
	   else
	      Mdot=D%Mdot
	   endif
	   do j=1,D%nTheta-1
	      dens(i,j)=C(i,j)%gasdens
	   enddo

	   !  Calculate viscous energy (Er) in a column corresponding to the mass accretion rate.
	   Er(i)=(2d0*sqrt(D%Rstar/(AU*D%R(i+1)))-3d0)/(3d0*D%R(i+1)*AU)
	   Er(i)=Er(i)-(2d0*sqrt(D%Rstar/(AU*D%R(i)))-3d0)/(3d0*D%R(i)*AU)
	   Er(i)=2d0*pi*Er(i)*3d0*G*D%Mstar*Mdot/(4d0*pi)

	   !  Calculate viscous energy (Eint) in a column for given surface
	   !   density, alpa and temperature
	   call VisHeatColumn(i,dens(i,1:D%nTheta-1),Temp(i,1:D%nTheta-1),Eint,Sig(i),fac,Er(i))

	   !  Rescale the density so Er=Eint
	   !
	   !  First for Eint > Er in a deadzone (scale down)
	   if (deadzone.and.D%MPdead(i).and.Eint.gt.Er(i)) then
	      surfscale(i)=fac
	      type="down"

	      if (Eint.gt.Er(i)*1.05d0) then 
c		 write(*,'("i=",i3," r=",f6.2," Eint=",e10.3," > Er=",e10.3," fac=",e10.3)') i,D%R_av(i)/AU,Eint,Er(i),fac
c		 print*,"downscaling not yet implemented!"
	      endif
	      
	   else if (deadzone.and.D%MPdead(i)) then
              !  MP is dead, assume extra mass ends up in deadzone
	      surfscale(i)=1d0 + (Er(i)-Eint)/(fac*Sig(i))
	      type="dead"

	   else if (deadzone.and.Er(i)/Eint*Sig(i)*gas2dust/2d0.gt.deadcolumn.and.Temp(i,D%nTheta-1).lt.deadtemp) then
!	   else if (deadzone.and.Er(i)/Eint*Sig(i)*gas2dust/2d0.gt.deadcolumn.and.D%MPzombie(i)) then
	      !  MP becomes dead after rescaling
	      !  Energy/Surface density in rest of active layer
	      Sigact=2d0*deadcolumn/gas2dust - Sig(i)
	      Eact=Sigact *Eint/Sig(i)

	      !  Energy/Surface density of deadzone
	      Edead=Er(i)-Eint-Eact
	      Sigdead=Edead/fac

	      !  Calculate scale factor from active and dead layer:
!	      surfscale(i)=1d0+(Sigact+Sigdead)/Sig(i)
	      surfscale(i)=(2d0*deadcolumn/gas2dust+Sigdead)/Sig(i)

	      type="zombie"
	      if (i.eq.itest) then
		 print*
		 print*,"  Sig:           ",Sig(i)
		 print*,"  Sig new:       ",Er(i)/Eint*Sig(i)
		 print*,"  New col:       ",Er(i)/Eint*Sig(i)*gas2dust/2d0
		 print*,"  deadcolumn:    ",deadcolumn
		 print*
	      endif

	   else 
              ! midplane is active
	      surfscale(i)=Er(i)/Eint ! equal to dens(i,j)/C(i,j)%gasdens
	      type="alive"

	   endif

	   ! Store the old and new surface density profile (check only)
	   if (testscale) then
	      Sig_old(i)=Sig(i)
	      Sig_lin(i)=Sig(i)*(Er(i)/Eint)
	      Sig_dead(i)=Sig(i)*(1d0 + (Er(i)-Eint)/(fac*Sig(i)))
	      Sig_used(i)=Sig(i)*surfscale(i)
	   endif

	   ! diagnostics
	   if (testall) then
	      write(*,'("i=",i3," r=",f6.2," Er=",e10.3," Eint=",e10.3," Sig=",e10.3,x,a)') 
     &          i,D%R_av(i)/AU,Er(i),Eint,Sig(i),type
	      if(i.eq.D%nR-1) print*
	   endif
	   if (testscale.and.i.eq.itest) then
	      write(*,'("1st: i=",i3," r=",f6.2," Er=",e10.3," Eint=",e10.3," Sig=",e10.3,x,a)') 
     &          i,D%R_av(i)/AU,Er(i),Eint,Sig(i),type
	   endif

	enddo ! i=1,D%nR-1

	!  Check/recalculate the surface density scaling in case of a deadzone
	!
	if(deadzone) then
	   
	   ! set density and calculate deadzone location
	   do i=1,D%nR-1
	      do j=1,D%nTheta-1
		 prevmass(i,j)=C(i,j)%mass
		 prevgasdens(i,j)=C(i,j)%gasdens
		 C(i,j)%gasdens=dens(i,j)*surfscale(i) ! used by MakeDeadZone
		 C(i,j)%mass=C(i,j)%mass*surfscale(i) ! used by VisHeatColumn
	      enddo
	   enddo

	   call MakeDeadZone(.false.) ! set to false (not printing location)

	   do i=1,D%nR-1

	      ! Recalculate viscous energy (Eint) in a column for new density and deadzone location
	      call VisHeatColumn(i,dens(i,1:D%nTheta-1)*surfscale(i),Temp(i,1:D%nTheta-1),Eint,Sig(i),fac,Er(i))

	      ! Eint > Er(i), scale down
	      if (deadzone.and.D%MPdead(i).and.Eint.gt.Er(i)) then
		 surfscale(i)=surfscale(i)*fac
		 type="down"

		 if (Eint.gt.Er(i)*1.05d0) then 
c		    write(*,'("i=",i3," r=",f6.2," Eint=",e10.3," > Er=",e10.3," fac=",e10.3)') i,D%R_av(i)/AU,Eint,Er(i),fac
c		    print*,"downscaling not yet implemented!"
		 endif

	      ! MP is dead, assume extra mass ends up in deadzone
	      else if (deadzone.and.D%MPdead(i)) then
		 surfscale(i)=surfscale(i)* (1d0 + (Er(i)-Eint)/(fac*Sig(i))) 
		 type="dead"

	      !  MP becomes dead after rescaling
	      else if (deadzone.and.Er(i)/Eint*Sig(i)*gas2dust/2d0.gt.deadcolumn.and.Temp(i,D%nTheta-1).lt.deadtemp) then
!	      else if (deadzone.and.Er(i)/Eint*Sig(i)*gas2dust/2d0.gt.deadcolumn.and.D%MPzombie(i)) then

	         !  Energy/Surface denisty in rest of active layer
		 Sigact=2d0*deadcolumn/gas2dust - Sig(i)
		 Eact=Sigact *Eint/Sig(i)

	         !  Energy/Surface density of deadzone
		 Edead=Er(i)-Eint-Eact
		 Sigdead=Edead/fac

	         !  Calculate scale factor from active and dead layer:
!	         surfscale(i)=1d0+(Sigact+Sigdead)/Sig(i)
		 surfscale(i)=surfscale(i)*(2d0*deadcolumn/gas2dust+Sigdead)/Sig(i)

		 type="zombie"

	      else ! midplane is active 
		 surfscale(i)=surfscale(i)* (Er(i)/Eint)
		 type="alive"
	      endif


	      ! diagnostics
	      !
	      if (testall) then
		 write(*,'("i=",i3," r=",f6.2," Er=",e10.3," Eint=",e10.3," Sig=",e10.3,x,a)') 
     &                 i,D%R_av(i)/AU,Er(i),Eint,Sig(i),type
		 if(i.eq.D%nR-1) print*
	      endif
	      if (testscale.and.i.eq.itest) then
		 write(*,'("2nd: i=",i3," r=",f6.2," Er=",e10.3," Eint=",e10.3," Sig=",e10.3,x,a)') 
     &                i,D%R_av(i)/AU,Er(i),Eint,Sig(i),type
	      endif

	      if (testscale.and.i.eq.itest) then
		 print*
		 write(*,'("Surface density was:     ",f16.4," gr/cm^2")') Sig_old(i)
		 write(*,'("     linear scaling:     ",f16.4," gr/cm^2")') Sig_lin(i)
		 write(*,'("   deadzone scaling:     ",f16.4," gr/cm^2")') Sig_dead(i)
		 write(*,'("       used scaling:     ",f16.4," gr/cm^2")') Sig_used(i)
		 write(*,'("       recalculated:     ",f16.4," gr/cm^2")') Sig_old(i)*surfscale(i)
	      endif

	      if (.true..and.Eint.gt.Er(i)*1.05d0) then
		 print*,"Houston, we have a problem!"
		 print*,"This shouldn't happen..."
		 write(*,'("i=",i3," r=",f6.2," Eint=",e14.3," <? Er=",e14.3)') i,D%R_av(i)/AU,Eint,Er(i)
		 stop
	      endif

	   enddo ! i=1,D%nR-1

	   ! Restore gasdens, mass
	   do i=1,D%nR-1
	      do j=1,D%nTheta-1
		 C(i,j)%mass=prevmass(i,j) !*surfscale(i) ! used by ToomreQ
		 C(i,j)%gasdens=prevgasdens(i,j) !*surfscale(i) ! used by ToomreQ
	      enddo
	   enddo

	endif ! deadzone

	! Lower surface density and (local) Mdot in case of gravitational instabilities
	!
	if (gravstable) then

           ! Radii where disk is gravitationally unstable
           imin=D%nR
           imax=0
	
	   do i=1,D%nR-1

	      !  Keep density above Q=2 for stability
	      Q=ToomreQ(i) / surfscale(i)
	      Q0=Q

	      if (Q.ge.2d0) then ! stable
c		 D%MdotR(i)=D%MdotR(i)
		 D%MdotR(i)=D%Mdot ! bug?
	      else ! unstable
		 surfscale(i)=surfscale(i)*(Q/2d0)
		 imin=min(i,imin)
		 imax=max(i,imax)

c		 D%MdotR(i)=D%MdotR(i)*(Q/2d0)
		 D%MdotR(i)=D%Mdot*(Q/2d0) ! bug?
	      endif	   	   
c	      write(*,'("r: ",f6.2," Q before: ",f12.3," Q after:",f12.3," Q stable:",f12.3)') 
c     &              D%R_av(i)/AU,ToomreQ(i),Q0,ToomreQ(i) / surfscale(i)
c	   write(*,'("r=",f6.2," Q=",f14.3," Mdot=",e14.3)') D%R_av(i)/AU,Q,D%MdotR(i)/Msun*(365.25d0*24d0*60d0*60d0)
	   enddo

	   !  Print if grav unstable
           if(imin.le.imax) then
	      write(*,'(a,f6.2," and ",f6.2," AU")') 
     &             "Grav unstable between:    ",D%R_av(imin)/AU,D%R_av(imax)/AU
	      write(9,'(a,f6.2," and ",f6.2," AU")') 
     &             "Grav unstable between:    ",D%R_av(imin)/AU,D%R_av(imax)/AU
           endif

	   !  Lower mass accretion rate interior to the deadzone?
c	   if (reducemdot) then
	   if (deadzone.and.reducemdot.and.imin.le.imax) then

c	      do i=D%nR-2,1,-1
	      do i=imax,1,-1

	         ! higher accretion rate then further out -> scale down dens, mdot
		 fac= D%MdotR(i+1)/D%MdotR(i)
		 if (fac.lt.1d0) then
		    surfscale(i)=surfscale(i)*fac
		    D%MdotR(i)=D%MdotR(i) *fac
		    if(i.eq.1) then
		       write(*,'("Reduced accretion rate: ",e14.3,"Msun/yr")') 
     &		             D%MdotR(i)/Msun*(365.25d0*24d0*60d0*60d0)
		       write(9,'("Reduced accretion rate: ",e14.3,"Msun/yr")') 
     &                       D%MdotR(i)/Msun*(365.25d0*24d0*60d0*60d0)
		    endif
		 endif
	      enddo

	   endif ! reducemdot
	endif ! gravstable
	
	! Update the density
	do i=1,D%nR-1
	   dens(i,1:D%nTheta-1)=dens(i,1:D%nTheta-1)*surfscale(i)
	enddo

	if (.not.getalpha) then ! Set the radial structure

	   DiskMass=0d0
	   do i=1,D%nR-1
	      do j=1,D%nTheta-1
		 if(C(i,j)%gasdens.gt.1d-50) then
		    C(i,j)%dens=C(i,j)%dens*surfscale(i)
		    C(i,j)%dens0=C(i,j)%dens0*surfscale(i)
		 else
		    C(i,j)%dens=1d-60
		    C(i,j)%dens0=1d-60
		    C(i,j)%w=C(i,j)%w0
		    C(i,j)%gasdens=1d-60
		 endif
		 C(i,j)%mass=C(i,j)%dens*C(i,j)%V
		 C(i,j)%gasdens=dens(i,j)
		 call CheckMinimumDensity(i,j)
		 DiskMass=DiskMass+C(i,j)%V*(gas2dust*C(i,j)%gasdens+C(i,j)%dens0)
	      enddo
	      
	   enddo

	   write(*,'("Total disk mass: ",f10.3," Msun")') DiskMass/Msun
	   write(9,'("Total disk mass: ",f10.3," Msun")') DiskMass/Msun

	else			! only store the radial structure
	   
	   write(file,'(a,"radialalpha.dat")') outdir(1:len_trim(outdir))
	   open(unit=66,file=file,RECL=1000)
	   write(66,*) "# Value of alpha for current surface density"
	   write(66,*) "# It prints radius, alpha"
	   do i=1,D%nR-1
	      write(66,*) D%R_av(i)/AU,surfscale(i) * alphavis
	   enddo
	   close(unit=66)

	endif
	
c================================
c Add a possible infalling cloud
c================================
c first remove cloud particles 
c from the disk
c================================
	do ii=1,ngrains
		if(Grain(ii)%shtype.eq.'INFALL') then
			do i=1,D%nR-1
				do j=1,D%nTheta-1
					C(i,j)%w(ii)=0d0
					infall_tot=sum(C(i,j)%w(1:ngrains))
					if(infall_tot.gt.1d-60) then
						C(i,j)%w=C(i,j)%w/infall_tot
					endif
					C(i,j)%w0(ii)=0d0
					infall_tot=sum(C(i,j)%w0(1:ngrains))
					if(infall_tot.gt.1d-60) then
						C(i,j)%w0=C(i,j)%w0/infall_tot
					endif
				enddo
			enddo
		endif
	enddo
c================================
c================================

	do ii=1,ngrains
		if(Grain(ii)%shtype.eq.'INFALL') then
			do i=1,D%nR-1
				do j=1,D%nTheta-1
					infall_mu=cos(D%theta_av(j))
					call solve_mu0(D%R_av(i)/(AU*D%Rexp),infall_mu,infall_mu0)
					if(infall_mu0.lt.D%mu0max) then
						infall_rho=((D%Mdot/(4d0*pi*sqrt(6.67300d-8*D%Mstar*D%R_av(i)**3)))
     &					*(1d0+infall_mu/infall_mu0)**(-0.5d0)*(infall_mu/infall_mu0+2d0*infall_mu0**2*D%Rexp*AU/D%R_av(i))**(-1d0))
					else
						infall_rho=1d-60
					endif
					infall_rho=infall_rho/gas2dust
					C(i,j)%dens=C(i,j)%dens*(1d0-C(i,j)%w(ii))
					C(i,j)%dens0=C(i,j)%dens0*(1d0-C(i,j)%w0(ii))
					C(i,j)%gasdens=C(i,j)%gasdens*(1d0-C(i,j)%w0(ii))
					C(i,j)%w=C(i,j)%w*C(i,j)%dens
					C(i,j)%w0=C(i,j)%w0*C(i,j)%dens0
					C(i,j)%w(ii)=infall_rho
					C(i,j)%w0(ii)=infall_rho
					C(i,j)%dens=C(i,j)%dens+infall_rho
					C(i,j)%dens0=C(i,j)%dens0+infall_rho
					C(i,j)%gasdens=C(i,j)%gasdens+infall_rho
					C(i,j)%mass=C(i,j)%dens*C(i,j)%V
					infall_tot=sum(C(i,j)%w(1:ngrains))
					C(i,j)%w=C(i,j)%w/infall_tot
					infall_tot=sum(C(i,j)%w0(1:ngrains))
					C(i,j)%w0=C(i,j)%w0/infall_tot
				enddo
			enddo
		endif
	enddo
c================================
c================================

	return
	end
	
c-----------------------------------------------------------------------
c A helper routine for RadialStruct() to calculate viscous energy (Eint)
c in a column (i) for a given density, alpa and temperature
c It returns 
c - Eint
c - Sig (surface density)
c - fac (if Eint < Er: Energy release per surface density in bottom cell
c        if Eint >= Er: scale factor)
c-----------------------------------------------------------------------
	subroutine VisHeatColumn(i,densR,TempR,Eint,Sig,fac,Er)
	use Parameters
	IMPLICIT NONE
	integer i,j
	real*8 densR(1:D%nTheta-1),TempR(1:D%nTheta-1)
	real*8 Eint,mu,G,alpha,Eint0,Er
	real*8 Sig,Mtot,Mtot0,fac
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
	parameter(G=6.67300d-8) ! in cm^3/g/s^2
	real*8 lininterpol,number_invalid ! function
	
	Eint=0d0
	Sig=0d0
	Mtot=0d0
	fac=0d0
	do j=1,D%nTheta-1
	   alpha=min(1d0,C(i,j)%alphavis)
	   Sig=Sig+densR(j)*(D%R_av(i)*(D%Theta(j)-D%Theta(j+1)) )*2d0
           Mtot0=densR(j)*C(i,j)%V*gas2dust
	   Mtot=Mtot+Mtot0
	   Eint0=sqrt(G*D%Mstar/D%R_av(i)**3)*9d0
     &		     *densR(j)*C(i,j)%V*gas2dust*alpha*kb*TempR(j)/(4d0*mu)
	   Eint=Eint+Eint0

	   ! store column where Er = Eint
	   if (fac.eq.0d0.and.Eint.gt.Er) then
	      fac=lininterpol(Er,Eint-Eint0,Eint,Mtot-Mtot0,Mtot)
	   endif

	   if(.false..and.i.eq.77) then
	      write(*,'("j=",i3," th=",f6.2," Eint=",e10.3," Sig=",e10.3)') 
     &                j,pi/2-D%Theta_av(j),Eint,Sig
	      if(j.eq.D%nTheta-1) then
c		 print*,G,D%Mstar,D%R_av(i),densR(j),C(i,j)%V,gas2dust,alpha,kb,TempR(j),mu,sqrt(G*D%Mstar/D%R_av(i)**3)*9d0
c     &		     *densR(j)*C(i,j)%V*gas2dust*alpha*kb*TempR(j)/(4d0*mu)
		 print*,j,densR(j),C(i,j)%V,alpha,TempR(j),sqrt(G*D%Mstar/D%R_av(i)**3)*9d0
     &		     *densR(j)*C(i,j)%V*gas2dust*alpha*kb*TempR(j)/(4d0*mu)
c		 write(*,'(" Eint=",e10.3," Sig=",e10.3)')
	      endif
	   endif

	   ! return paramter fac for density scaling
	   if(j.eq.D%nTheta-1) then

	      if (Eint.gt.Er) then
		 !  this fraction of the surface density generates Er.
		 fac=fac/Mtot
	      	 
	      else	 
	         ! heating per unit surface density in bottom cell, assuming it is a deadzone		
		 alpha=deadalpha/prandtl
		 Eint0=sqrt(G*D%Mstar/D%R_av(i)**3)*9d0
     &		    *densR(j)*C(i,j)%V*gas2dust*alpha*kb*TempR(j)/(4d0*mu)
		 Mtot0=densR(j)*C(i,j)%V*gas2dust
		 fac=Eint0/Mtot0*(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)*gas2dust ! why?????
c	      fac=Eint0/(densR(j)*(D%Theta(j)-D%Theta(j+1))) ! wrong
	      endif
	   endif
c	   print*,D%Mstar,D%R_av(i),densR(j),C(i,j)%V,gas2dust,alpha,kb,TempR(j)
	enddo

c	print*,Eint/Mtot,fac
c	stop66167

	end subroutine

