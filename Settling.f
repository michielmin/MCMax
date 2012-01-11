c     --------------------------------------------------------------
c     This routine computes the vertical dust structure selfconsistently
c     using an _equilibrium_ between dust settling and turbulent mixing.
c     It is based on the f90 routine written by Kees Dullemond.
c     It uses implicit integration (backward Euler) to get to the final 
c     solution in one timestep.
c
c     NOTE: The zgrid used here goes in the opposite direction as the 
c           thetagrid used in mcmax.
c     NOTE: The diffusion equation assumes cylindrical coordinates, this
c           might introduce an extra error.
c     
c     TODO: - clean up: remove print statements and clear obsolete parts
c
c     --------------------------------------------------------------

      subroutine Settling(rho,nz,ns,ir)
      use Parameters
      implicit none
      integer nz,ns,ir
      doubleprecision rho(nz,0:ns)
      !
      !  Defined elsewhere in MCMax
      !
      doubleprecision radius,mstar
      doubleprecision zi(nz+1),zc(nz) ! i=interface, c=center of cell
      doubleprecision tgas(nz)
      doubleprecision mgrain(ns),agrain(ns)
      !
      !  Local variables
      !                         
      integer iz,is,ith
      integer its,nts,isteady
      integer logtstart, logtend
      integer ibottom
      doubleprecision f(nz,ns),fold(nz,ns)
      doubleprecision rhogas(nz)
      doubleprecision omk,vkep,diffcoef(1:ns,1:nz+1),vsett(1:ns,1:nz+1),alpha(nz)
      doubleprecision soundspeed(1:nz),soundspeed2(1:nz),stokes(1:ns,1:nz)
      doubleprecision diffnew(1:ns,1:nz+1)
      doubleprecision y(1:nz),yold(1:nz),dzc(1:nz),dzi(1:nz+1),mugas,rhog,tempg,sigfr
      doubleprecision rhs(1:nz),eps,eps1,force
      doubleprecision sigmadust,sigmadust_old,mcol,mcol_old,cellv(nz)
      doubleprecision tr_a(1:nz),tr_b(1:nz),tr_c(1:nz) !c interferes with MCMax
      doubleprecision dtime,convergetime(ns),dif,difmax,stoptolyear
      doubleprecision period,tsettmin,ageoftheuniverse,logstep
      doubleprecision scalebottom
      doubleprecision, allocatable :: timestep(:)
      doubleprecision, PARAMETER :: GG=6.6720000d-08
      doubleprecision, PARAMETER :: amu=1.66053886d-24
      doubleprecision, PARAMETER :: year=3.1556926d7
      character*100 outputfile, rdisk
      logical makeplot, condition
      logical hitbottom
      !
      !  testing only
      !
      logical stupid_test
      !
      !  Functions
      !
      doubleprecision vdrift
      integer number_invalid
      !
      !  testing only
      !
      stupid_test=.false.        ! test leveling of at low densities
      !
      !  Set dust parameters
      !
      do is=1,ns
         agrain(is)=Grain(is)%rv
         mgrain(is)=Grain(is)%m
      enddo
      goto 333  ! skip some print statements
      !
      !  Print dust parameters
      !
      write(*,'("  ")')
      write(*,'("ns= ",I3)') ns
      write(*,'("  ")')
      do is=1,ns
         write(*,'("agrain[",I0,"]= ",F10.4)') is,Grain(is)%rv * 1d4
         write(*,'("mgrain[",I0,"]= ",E16.5)') is,Grain(is)%m
      enddo
 333  continue                  ! half evil
      !
      !  Set global disk parameters
      !
      radius=D%R_av(ir)
      mstar=D%mstar
      !
      !  Set disk parameters at grid interface
      !
      do iz=1,nz+1
         zi(iz)=radius * (pi/2d0 - D%thet((nz+1)-(iz-1)) ) ! reverse order     
      enddo
      !
      !  Now make the dz over the interfaces
      !
      do iz=1,nz+1
         if (iz.eq.1) then 
            dzi(1)    = 0.5d0 * ( zi(2) - zi(1) )
         else if (iz.eq.nz+1) then 
            dzi(nz+1) = 0.5d0 * ( zi(nz+1) - zi(nz) )
         else
            dzi(iz) = 0.5d0 * ( zi(iz+1)-zi(iz-1) )
         endif
c         write(*,'("zi[",I0,"]  = ",E16.5,"  | dzi[",I0,"]  = ",E16.5)') 
c     1        iz,zi(iz),iz,dzi(iz)   
      enddo
      !
      !  Set disk parameters at grid center
      !
      do iz=1,nz
         zc(iz) = 0.5d0 * ( zi(iz) + zi(iz+1)  )
         dzc(iz) = zi(iz+1)-zi(iz)
         ith=(nz)-(iz-1)
         cellv(iz)=C(ir,ith)%V  ! ith: reverse order
         tgas(iz)=C(ir,ith)%T   ! ith: reverse order
         rhogas(iz)=rho(ith,0) * gas2dust  ! ith: reverse order
         !
         !  Check for zero gas densities (seems to be fixed)
         !
c         if (rhogas(iz).eq.0d0) write(*,*) "Aarg at ir=",ir," and iz=", iz
c         if (rhogas(iz).eq.0d0) stop 
c     1        "Encountered zero gas density from Diskstructure.f"
         if (rhogas(iz).eq.0d0) rhogas(iz)=1d-300
         do is=1,ns
            f(iz,is)=rho(ith,is) ! ith: reverse order
            fold(iz,is)=f(iz,is)
         enddo
         !
         !  Alpha: global or some gradient (deadzone)
         !
         if (deadzone .and. zc(iz) / radius.le.deadheight) then
c$$$               write(*,'("Deadzone in cell [",i3,",",i2,"]")') ir,iz 
            alpha(iz)=deadalpha
         else
            alpha(iz)=alphaturb ! no radial dependence alfa, for now
         endif
      enddo
      goto 666                  ! skip some print statements
      !
      !  Print disk parameters
      !
      write(*,'("  ")')
      write(*,'("Calculating at r= ",F6.3)') radius / AU
      write(*,'("M= ",F6.3)') mstar / Msun
      write(*,'("  ")')
      write(*,'("nz= ",I3)') nz
      write(*,'("gas2dust= ",F6.1)') gas2dust
      write(*,'("  ")')
c      do iz=1,nz+1
c         write(*,'("Theta_edge[",I0,"]= ",F9.6)') 
c     &        iz,(pi/2. - D%thet((nz+1)-(iz-1)))          ! reverse order
c         write(*,'("zi[",I0,"]= ",E16.5)') iz,zi(iz)
c      enddo
      do iz=1,nz
c         write(*,'("Link ",I0," to ",I0)') iz,ith
c         write(*,'("z[",I0,"]  = ",E16.5)') iz,z(iz)
         write(*,'("z/r= ",F4.2)') zc(iz)/ radius
         write(*,'("Tgas[",I0,"]= ",F6.1)') iz,tgas(iz)  !constant before iter
         write(*,'("rhodust[",I0,"]= ",E16.5)') 1, f(iz,1)
         write(*,'("rhodust[",I0,"]= ",E16.5)') 2, f(iz,2)
         write(*,'("rhodust[",I0,"]= ",E16.5)') 3, f(iz,3)
         write(*,'("rhogas= ",E16.5)') rhogas(iz)
         write(*,'("  ")')
      enddo
 666  continue
      !
      ! Make the gravitational force constant
      !
      omk    = sqrt( GG * mstar / radius**3d0 )
      vkep   = sqrt( GG * mstar / radius )
c      write(*,'("omk= ",E16.5)') omk
c      write(*,'("vkep= ",E16.5)') vkep
      mugas=2.3
      !
      !  Settle to equilibrium or do time-dependent?
      !
      if (scseteq) then
c$$$         nts=1
c$$$         allocate(timestep(0:1))
c$$$         timestep(0)=0d0
c$$$         timestep(1)=1d10           ! make sure we reach an equilibrium
c$$$         convergetime(1:ns)=1d10
         !
         nts=10                            ! Always do a few timesteps
         ageoftheuniverse=13.73d9          ! need some large nr to  reach an equilibrium
         !
         allocate(timestep(0:nts))
         timestep(0)=0d0
         do its=1,nts
            timestep(its)=10d0**(its-nts) * ageoftheuniverse
         enddo
         convergetime(1:ns)=ageoftheuniverse
      else
         !
         !  Defaults to lifetime
         !
         convergetime(1:ns)=lifetime
         !
         !  Make the timegrid for timedependent settling:
         !
         !
         !  First calculate orbital period (2 pi r / vkep = 2 pi / omk)
         !
c         write(*,*) radius / AU, 2. * pi / (omk * year)
         period = 2d0 * pi / (omk * year)
         !
         !  Time to travel from surface to midplane at maximum speed (z/r vkep)
         !
         tsettmin = period / 4.
         tsettmin = tsettmin
c         write(*,*) radius / AU, tsettmin
         !
         !  Size of the timegrid (log!), using 10 steps per order of magnitude
         !
         !  NOTE:  extra timesteps for:
         !         Rounding of decimal for negative nr (-1)
         !         Sample the interval [0, ..., tsettmin] (-11)
         !         tend > lifetime (+0)   ; already in nts i.o nts-1
         !         convergence (+3)
c         write(*,*) tsettmin
c         write(*,*) log10(tsettmin)
         logtstart=  int(10d0*log10(tsettmin))-1
         logtend=    int(10d0*log10(lifetime))+3
c         write(*,*) logtstart,logtend
         nts = logtend - (logtstart - 11)
         !
         allocate(timestep(0:nts))
         timestep(0)=0d0
         !
         !  Sample the interval [0,tsettmin] properly
         !    using 10 steps per order of magnitude
         !    and making sure first step is the biggest
         !
         do its=1,9
c         do its=0,10
            timestep(its) = ((10d0**(its / 10d0) -1d0) / 9d0 ) * 10d0**(logtstart / 10d0)
c            write(*,*) timestep(its) 
         enddo
         !
         !  Sample the interval [tsettmin, lifetime] + extra steps for convergence
         !
         do its=10,nts
            timestep(its) = 10d0**( (its-10+logtstart)   / 10d0)   ! *10**logtstart
c            write(*,*) timestep(its) 
         enddo
         !
         !  Stupidity checks Kees-style
         !
         if (timestep(10).le.timestep(9)) stop "grid doesnt connect"
         if (timestep(nts).le.log10(lifetime)) stop "grid too small"
         !
         !  Skip this part? (only for debugging eps=0.5, eps=0)
         !
         goto 616
         !
         !  Print the time-grid 
         !
         write(outputfile,'(a,"settle_timegrid.dat")') outdir(1:len_trim(outdir))
c         if (ir.eq.100) then
         if (ir.eq.1) then
            open(unit=66,file=outputfile,RECL=1000)
            !
            !  Write some file info
            !
            write(66,*) "# scset: time grid"
            write(66,*) "# It prints its,time"
            do its=0,nts
               write(66,88) its,timestep(its)
 88            format(I3,' ',D16.6)
            enddo
         endif
         !
         !  Close the ouput file
         !
         close(unit=66)
         !
         !  If we skipped outputting the time-grid
         !
 616     continue               ! truly evil
      endif
      !
      !  When is model converged?
      !    From Kees code: 1d8
      !    Based on radial tau=1 convergence in TTs : 
      !    1d5:     too short
      !    1d6:     ok
      !    1d7-1d8  way too long
      !
      stoptolyear=1d6
      !
      !  Now do the settling in a vertical slab, for each particle
      !
      do is=1,ns
         !
         !  Do not use settling if particle abundance=0 
         !  NOTE: (check is only in midplane !)
         !
         if (C(ir,nz)%w(is).eq.0d0) then
c            write(*,*) "Species", is, "is zero!"
            convergetime(is)=0d0
            goto 222
         endif
         !
         !  Compute the frictional cross section
         !
         sigfr     = pi*agrain(is)**2d0
         !
         !  Compute Soundspeed and Stokes number
         !
         do iz=1,nz
            soundspeed2(iz) = kb*tgas(iz) / (mugas*amu)
            soundspeed(iz) = sqrt(soundspeed2(iz))
            !
            stokes(is,iz) = 0.75d0* mgrain(is)/sigfr * 
     1           omk/(rhogas(iz)*soundspeed(iz)) *
     2           alpha(iz)**(2d0*qturb - 1d0)
            !
c            if (iR.ge.2.and.D%R(iR-1).lt.1.0.and.D%R(iR).ge.1.0) then
c               if (iz.eq.1) then
c                  write(*,*)
c                  write(*,*) 'is= ',is
c                  write(*,*) alpha(iz),(2*qturb - 1),alpha(iz)**(2*qturb - 1)
c                  write(*,*)
c               endif
c               write(*,*) 'z/r,stokes:',zc(iz)/radius,stokes(is,iz)
c               write(*,*) 'z/r,1/(1+stokes):',zc(iz)/radius,1/(1+stokes(is,iz))
c            endif
         enddo
         !
         ! Compute the turbulent diffusion coefficient D 
         !   at the interface of the cells
         !
         ! Note: diffcoef is D * rhogas / (1 + St)
         ! D = alpha * soundspeed^2 / kepplerfrequency
         !
         do iz=2,nz
            !
            !  Old diffusion coefficient
            !
c            diffcoef(is,iz) = 0.5 * kb * ( rhogas(iz) * alpha(iz) * tgas(iz) +
c     1           rhogas(iz-1) * alpha(iz-1) * tgas(iz-1) ) / 
c     2           ( mugas*amu*omk )
            !
            !  And the new one
            !
            diffcoef(is,iz)= 0.5d0* ( rhogas(iz)  * alpha(iz)  * soundspeed2(iz)  / omk 
     1                              / ( 1d0 + stokes(is,iz)**2d0 )                           )
     2           +           0.5d0* ( rhogas(iz-1)* alpha(iz-1)* soundspeed2(iz-1)/ omk
     3                              / ( 1d0 + stokes(is,iz-1)**2d0 )                         )
            !
            !  print some stuff
            !
c            if (iR.ge.2.and.D%R(iR-1).lt.1.0.and.D%R(iR).ge.1.0) then
c               if (iz.eq.2) then
c                  write(*,*)
c                  write(*,*) 'is= ',is
c                  write(*,*)
c               endif
c               if ((diffcoef(is,iz)/diffnew(is,iz)-1).ge.1) then
c                  write(*,*) 'z/r,diffcoef old,new, stokes:',
c     1                 zc(iz)/radius,diffcoef(is,iz),diffnew(is,iz),stokes(is,iz)
c               endif
c            endif
         enddo
         diffcoef(is,1)    = 0.d0  ! Flux at the inner boundary zero
         diffcoef(is,nz+1) = 0.d0  ! Flux at the outer boundary zero
         !
         ! Do a loop over interfaces of the cells
         !
         do iz=2,nz
            !
            ! Compute temp and density of gas at the interface 
            ! using interpolation
            !
            tempg        = 0.5d0 * ( tgas(iz-1) + tgas(iz) )
            rhog         = 0.5d0 * ( rhogas(iz-1) + rhogas(iz) )
c            write(*,*) rhogas(iz-1), rhogas(iz)
c            write(*,*) rhog
            !
            ! Compute the force on the grain at the interface
            !
            force        = mgrain(is)*omk*omk*zi(iz)*(1d0+C(ir,(nz)-(iz-1))%FradZ*gas2dust)
            !
            ! Compute the equilibrium settling velocity at the interface
            !
c            write(*,*) tempg,rhog,mugas,force,sigfr
            vsett(is,iz) = -vdrift(tempg,rhog,mugas,force,sigfr)
            !
            ! NOTE: Warning: for very large particles this could be larger
            !       than (z/r)*v_Kepler, which is unphysical. The physics 
            !       here is that for such particles the equilibrium velocity 
            !       is never reached    !
            !       Trick: just limit it to (z/r)*v_kepler.
            !
            if (-vsett(is,iz).ge.vkep*(zi(iz)/radius)) then
                vsett(is,iz)=-vkep * (zi(iz)/radius)
c               write(*,'("vsett[",I0,",",I0"]= ",E16.5,A)') 
c     1              is,iz,vsett(is,iz),"  (=-z/r * vkep)"
            else 
c               write(*,'("vsett[",I0,",",I0"]= ",E16.5)') is,iz,vsett(is,iz)
            endif
            !
            ! Check validity of this number 
            !
            if(number_invalid(vsett(is,iz)).ne.0.and.
     &			Grain(is)%shtype.eq.'DISK'.and.Grain(is)%shscale(ir).gt.0d0) then
               write(*,*) 'ERROR: NaN or Inf value detected (2)'
               stop 
            endif
         enddo
         !
         ! Do the boundaries
         !
         vsett(is,1)    = 0.d0
         vsett(is,nz+1) = vsett(is,nz)
         !
         !  Now start the settling calculation
         !
         !  First copy the vertical density structure into a 1-D array
         !
         do iz=1,nz
c            y(iz) = f(is,iz) ! found the bug...
            y(iz) = f(iz,is)
         enddo
         !
         !  What kind of settling? (equilibrium or not)
         !
         if (scseteq) then
            !
            !  Don't use Crank-Nicholson method! ( i.e. eps=0.5)
            !  Implicit works better when we take one big timestep.
            !  (eps=0 --> explicit ; eps=1 --> implicit ; eps=0.5 --> crank-nich)
            !
            eps  = 1d0
         else
            !
            !  Can use Crank-Nicholson here? (not yet, need smaller timesteps?)
            !
c            eps  = 0.5d0            ! --> Works! But convergence criterium fails...
c            eps  = 0d0              ! --> Doesnt work, probably need smaller timesteps...
            eps  = 1d0
         endif
         !
         eps1 = 1.d0-eps
         !
         !  For checking convergence
         !
         isteady=0
         !
         !  Start the timeloop
         !
         do its=1,nts
            !
            !  Save the old vertical structure (for checking convergence)
            !
            do iz=1,nz
               yold(iz) = y(iz)
            enddo
            !
            !  Calculate dt for _this_ timestep
            !
            dtime=(timestep(its)-timestep(its-1)) * year
c            write(*,*) its, timestep(its), dtime / year
            !
            ! Fill the tridiagonal matrix elements a,b,c 
            !
            tr_a(1) = 0.d0
            do iz=2,nz
               tr_a(iz) = diffcoef(is,iz)/(dzc(iz)*dzi(iz)*rhogas(iz-1))
c               tr_a(iz) = diffcoef(is,iz)/(dzc(iz-1)*dzi(iz)*rhogas(iz-1))
            enddo
            do iz=1,nz
               tr_b(iz) = -(1.d0/(dzc(iz)*rhogas(iz))) * 
     1              ( diffcoef(is,iz)/dzi(iz) + diffcoef(is,iz+1)/dzi(iz+1)  )
     2              + vsett(is,iz)/dzc(iz)
            enddo
            tr_c(nz) = 0.d0
            do iz=1,nz-1
               tr_c(iz) = diffcoef(is,iz+1)/(dzc(iz)*dzi(iz+1)*rhogas(iz+1))
c               tr_c(iz) = diffcoef(is,iz+1)/(dzc(iz+1)*dzi(iz+1)*rhogas(iz+1))
     1           - vsett(is,iz+1)/dzc(iz)
            enddo
            !
            ! Now fill the right hand side of the matrix equation...
            !
            do iz=2,nz-1
               rhs(iz) = eps1*dtime*tr_a(iz)*y(iz-1) + 
     1              (1.d0+eps1*dtime*tr_b(iz))*y(iz) +
     2              eps1*dtime*tr_c(iz)*y(iz+1)
            enddo
            rhs(1)  =     (1.d0+eps1*dtime*tr_b(1))*y(1) + 
     1           eps1*dtime*tr_c(1)*y(2)
            rhs(nz) =     eps1*dtime*tr_a(nz)*y(nz-1) + 
     1           (1.d0+eps1*dtime*tr_b(nz))*y(nz) 
            !
            ! Now compute Mnew = 1-eps*dtime*M, which is the matrix we need to
            ! 'invert'.
            !
            do iz=1,nz
               tr_a(iz) = -eps*dtime*tr_a(iz)
               tr_b(iz) = 1.d0-eps*dtime*tr_b(iz)
               tr_c(iz) = -eps*dtime*tr_c(iz)
            enddo
            !
            ! Renormalize all
            !
            do iz=1,nz 
               tr_a(iz)   = tr_a(iz)/tr_b(iz)
               tr_c(iz)   = tr_c(iz)/tr_b(iz)
               rhs(iz)    = rhs(iz)/tr_b(iz)
               tr_b(iz)   = 1.d0
            enddo
            !
            ! Now call tridag to solve the matrix equation
            !
            call tridag(tr_a,tr_b,tr_c,rhs,y,nz)
            !
            !  Is solution sane? --> no!
            !
            if (stupid_test .and. ir .eq. 87 .and. is.eq.5) then
               write(*,'("Calculating at r= ",F6.3)') radius / AU
               open(unit=66,file="stupid_test",RECL=2000)
               write(66,'("# stupid test: directly after tridag ")')
               do iz=1,nz
                  write(66,'(e16.8,e16.8)') zi(iz) / radius, y(iz)
               enddo
               close(unit=66)
            endif
            !
            ! Check for negative densities
            !
            do iz=1,nz
               if(y(iz).lt.0.d0) y(iz)=0.d0
            enddo
            !
            !  Renormalize the dust density using cell volume (conserves mass)
            !  Using the surface density does not work in a spherical grid.
            !  The max difference in cell width is only 10% at z/r = 0.5
            !  so this should not affect the settling routine much
            !
            mcol=0d0
            mcol_old=0d0
            do iz=1,nz
               mcol=      mcol     + cellv(iz)*y(iz)
               mcol_old=  mcol_old + cellv(iz)*fold(iz,is)
c               write(*,*) zc(iz)/radius, cellv(iz) / dzc(iz) ! check cellwidth
            enddo
            do iz=1,nz
               y(iz)=y(iz)*mcol_old/mcol
               !
               !  Set zero densities to 1d-100
               !  TODO: - check that this doesn't screw up the mass conservation
               !
               !  NOTE: MCMax should handle this
               !
c               if (f(iz,is) .eq. 0d0) f(iz,is)=1d-100
            enddo
            !
            !  Calculate convergence criterium
            !
            difmax=0.d0
            do iz=1,nz
c               if (yold(iz).gt.1d-100.and.y(iz).gt.1d-100) then
               dif = abs(yold(iz)-y(iz)) / abs(yold(iz)+y(iz)+1d-99)
               if(dif.gt.difmax) then
                  difmax = dif
               endif
            enddo
            !
            !  Check for convergence (of one particle)
            !
            if(dtime/(difmax+1d-99).gt.stoptolyear*year) then
               isteady = isteady +1 
            endif
c            write(*,*) difmax,dtime/(difmax+1d-99),stoptolyear*year
            !
            !  Save convergence time and escape from loop
            !
            if(isteady.ge.3) then
c               write(*,*) 'Converged at ',timestep(its),' years'
               convergetime(is) = timestep(its-3)
               goto 999
            endif
            !
            !  End the timeloop
            !
         enddo
 999     continue               ! don't ask...
         !
         ! Copy back to f
         !
         do iz=1,nz
c            f(is,iz) = y(iz) ! found the bug...
            f(iz,is) = y(iz)
         enddo
         !
         !  Arrive here if abundance=0
         !
 222     continue
         !
         !  End the particle loop
         !
      enddo
      !
      !  Write convergence time at this distance
      !
      if (scsetsave) then
         !
         !  Open a new file for writing, if ir=0
         !
         write(outputfile,'(a,"settle_eq.dat")') outdir(1:len_trim(outdir))
         if (ir.eq.1) then
            open(unit=66,file=outputfile,RECL=1000)
            !
            !  Write some file info
            !
            write(66,*) "# Output file made if scsetsave=.true."
            write(66,*) "# It prints nr,ns"
            write(66,99) "# ", D%nR-1,ns
 99         format(A2,' ',I6,' ',I3)
            write(66,*) "# It prints r, time to reach equilibrium (ns)"
         else
            open(unit=66,file=outputfile,RECL=1000,ACCESS='append')
         endif
         !
         !  Write time to reach equilibrium
         !
c         write(66,*) radius / AU, convergetime(1:ns)  ! plot radius as well?
         write(66,*) convergetime(1:ns) 
         !
         !  Close the ouput file
         !
         close(unit=66)
      endif
      !
      !  Fix abundance at densities lower than 1d-50
      !  TODO: make sure this works with halo's etc
      !  TODO: make sure it also works at larger radii (normalization?)
      !
      if (thinparticle.ge.-1.and.thinparticle.le.ngrains) then
         hitbottom=.false.
!     condition=ir.le.135.and.ir.ge.130
         condition=0
         if (condition) write(*,*) "                                 ", radius / AU

         do iz=1,nz
            if (condition) write(*,'(e10.3,$)') sum(f(iz,1:ns))
            if (sum(f(iz,1:ns)).le.1d-50) then !  arrive at <1d-50 ?
               hitbottom=.true.
               ibottom=iz-1
               scalebottom=1.01d-50 / sum(f(iz-1,1:ns))
            endif
            
            if (hitbottom) then !  >1d-50 for everything above
               if (thinparticle.eq.-1) then ! use last density
                  f(iz,1:ns)=f(ibottom,1:ns) * scalebottom
               else if (thinparticle.eq.0) then ! use input abundance
                  f(iz,1:ns)=C(ir,nz)%w(1:ns)*1.01d-50
               else if (thinparticle.ge.1) then ! use particle # 
                  f(iz,1:ns)=0d0
                  f(iz,thinparticle)=1.01d-50
               endif
            endif
            if (condition) write(*,'(e10.3)') sum(f(iz,1:ns))
         enddo
      endif
      !
      !  Store the new dust distribution in rho
      !
      do iz=1,nz
         ith=(nz)-(iz-1)
         do is=1,ns
            if(Grain(is)%shtype.eq.'DISK'.and.Grain(is)%shscale(ir).gt.0d0) then
               rho(ith,is)=f(iz,is) ! ith: reverse order
            endif
         enddo
      enddo
      !
      !  Deallocate arrays
      !
      deallocate(timestep)
      !
      !  We did it!
      !
      end

      
!-------------------------------------------------------------------
!                          DRIFT VELOCITY
!
! This routine computes the drift velocity of a dust particle
! with respect to the gas. As input it needs the external force
! fext and the cross-section to gas drag sigfric. 
!-------------------------------------------------------------------
      doubleprecision function vdrift(tempg,rhog,mugas,force,sigfr)
      use Parameters
      implicit none
      doubleprecision tempg,rhog,mugas,force,sigfr
      !
      doubleprecision kk,mp,cc
      parameter(kk  = 1.3807d-16) ! Bolzmann's constant     [erg/K]
      parameter(mp  = 1.6726d-24) ! Mass of proton          [g]
      parameter(cc  = 2.9979d10) ! Light speed             [cm/s]
      doubleprecision cs,cst,dum1
      parameter(cst=1.13176848421)
      !
      if(rhog.lt.1d-60) then
         vdrift = cc            ! A ridiculously large velocity
         return
      endif
      cs     = sqrt(kk*tempg/(mugas*mp))
      vdrift = (3./4.) * force / (rhog*sigfr*cs) / sqrt(8d0/pi)
      if(vdrift.gt.cc) then
         vdrift=cc
      endif
      return
      end



!-------------------------------------------------------------------
!                 NUMERICAL RECIPES ROUTINE TRIDAG()
!      (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..
!-------------------------------------------------------------------
      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n,NMAX
      DOUBLEPRECISION a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=500)
      INTEGER j
      DOUBLEPRECISION bet,gam(NMAX)
      if(b(1).eq.0.)stop 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j)*gam(j)
         if(bet.eq.0.)stop 'tridag failed'
         u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
      return
c      END SUBROUTINE tridag
      end  ! weird errors

!-------------------------------------------------------------------
!                   FUNCTION: TEST VALIDITY OF NUMBER
!
!     0 = Number is okay
!     1 = Number is INF
!     2 = Number is NAN
!-------------------------------------------------------------------
      function number_invalid(a)
      implicit none
      integer number_invalid
      doubleprecision a,b,c
      logical div,sub
      !
      b=a*2.d0
      b=b/2.d0
      c=a-1.d100
      !
      div = (b.eq.a)
      sub = (c.lt.a)
      !
      if(div.and.sub) then
         number_invalid = 0
      elseif(div) then
         number_invalid = 1
      else
         number_invalid = 2
      endif
      return
      end 
