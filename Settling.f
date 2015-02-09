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
c     --------------------------------------------------------------

      subroutine Settling(rho,nz,ns,ir)
      use Parameters
      implicit none
      integer nz,ns,ir
      doubleprecision rho(nz,0:ns)

      !  Defined elsewhere in MCMax
      doubleprecision radius,mstar
      doubleprecision zi(nz+1),zc(nz) ! i=interface, c=center of cell
      doubleprecision tgas(nz)
      doubleprecision mgrain(ns,nz),agrain(ns,nz)

      !  Local variables
      integer iz,is,ith
      integer its,nts,isteady
      integer logtstart, logtend
      integer ibottom
      doubleprecision f(nz,ns)
      doubleprecision rhogas(nz)
      doubleprecision omk,vkep,diffcoef(1:ns,1:nz+1),vsett(1:ns,1:nz+1),alpha(nz)
      doubleprecision soundspeed(1:nz),soundspeed2(1:nz),stokes(1:ns,1:nz)
      doubleprecision diffnew(1:ns,1:nz+1)
      doubleprecision y(1:nz),yold(1:nz),dzc(1:nz),dzi(1:nz+1),mugas,rhog,tempg,sigfr(1:nz)
      doubleprecision rhs(1:nz),eps,eps1,force
      doubleprecision sigmadust,sigmadust_old,mcol,mcol_old,cellv(nz)
      doubleprecision tr_a(1:nz),tr_b(1:nz),tr_c(1:nz) !c interferes with MCMax
      doubleprecision dtime,settletime(ns),dif,difmax,stoptolyear
      doubleprecision period,tsettmin,logstep
      doubleprecision scalebottom
      doubleprecision, allocatable :: timestep(:)
      doubleprecision, PARAMETER :: GG=6.6720000d-08
      doubleprecision, PARAMETER :: amu=1.66053886d-24
      doubleprecision, PARAMETER :: year=3.1556926d7
      doubleprecision, PARAMETER :: ageoftheuniverse=13.73d9
      character*100 outputfile, rdisk
      logical makeplot, condition
      logical hitbottom, do_settling, do_settle(1:ns)

      !  Debugging/testing only
      logical stupid_test, print_dust, print_disk, print_time

      !  Functions
      doubleprecision vdrift
      integer number_invalid

      do iz=1,nz
         ith=(nz)-(iz-1)
         do is=0,ns
	            if(number_invalid(rho(ith,is)).ne.0) then
   	    	        write(*,*) 'ERROR: NaN or Inf value detected (2)'
					print*,ith,is,mgrain(is,iz),agrain(is,iz)
   	        	    stop 
   		         endif
			if(rho(ith,is).lt.1d-50) rho(ith,is)=1d-50
         enddo
      enddo

c-------------------------------------------------------------------
c     Here starts the initialisation
c-------------------------------------------------------------------

      !  Check which grains need to be settled
      !    (disk particle, not a wedge, abundance not zero)
      do is=1,ns
         do_settle(is)=(Grain(is)%shtype.eq.'DISK'.and.
     1        Grain(is)%shscale(ir).gt.0d0.and.C(ir,nz)%w(is).ne.0d0 )
         if (do_settle(is)) do_settling=.true.
      enddo
*      if (.not.do_settling) print*,' dont settle anything'
      if (.not.do_settling) goto 333     ! skip entire routine

      !  Debugging/testing only
      stupid_test=.false.        ! test leveling of at low densities
      print_dust=.false. 
      print_disk=.false.
      print_time=.false.

      !  Set dust size and mass
      !
	do iz=1,nz
      do is=1,ns
         call ComputeRM_nopac(is,ir,(nz+1)-(iz-1),agrain(is,iz),mgrain(is,iz))
      enddo
	enddo
      !  Print dust parameters ?
      if (print_dust) then
         write(*,'("  ")')
         write(*,'("ns= ",I3)') ns
         write(*,'("  ")')
         do is=1,ns
            write(*,'("agrain[",I0,"]= ",F10.4)') is,Grain(is)%rv * 1d4
            write(*,'("mgrain[",I0,"]= ",E16.5)') is,mgrain(is,1)
         enddo
      endif

      !  Set global disk parameters
      !
      radius=D%R_av(ir)
      mstar=D%mstar

      !  Set up the grid at cell interfaces
      !
      do iz=1,nz+1
         zi(iz)=radius * (pi/2d0 - D%thet((nz+1)-(iz-1)) ) ! reverse order     
      enddo

      !  Now make the dz over the interfaces
      do iz=1,nz+1
         if (iz.eq.1) then 
            dzi(1)    = 0.5d0 * ( zi(2) - zi(1) )
         else if (iz.eq.nz+1) then 
            dzi(nz+1) = 0.5d0 * ( zi(nz+1) - zi(nz) )
         else
            dzi(iz) = 0.5d0 * ( zi(iz+1)-zi(iz-1) )
         endif
      enddo

      !  Set up the grid at cell center and disk parameters
      !
      do iz=1,nz
         zc(iz) = 0.5d0 * ( zi(iz) + zi(iz+1)  )
         dzc(iz) = zi(iz+1)-zi(iz)
         ith=(nz)-(iz-1)
         cellv(iz)=C(ir,ith)%V  ! ith: reverse order
         tgas(iz)=C(ir,ith)%T   ! ith: reverse order
         rhogas(iz)=rho(ith,0) * gas2dust  ! ith: reverse order

         !  Check for zero gas densities
         if (rhogas(iz).eq.0d0) rhogas(iz)=1d-300

         !  Convert the input dust density from MCMax
         do is=1,ns
            f(iz,is)=rho(ith,is) ! ith: reverse order
         enddo
         
         !  Alpha: global or some gradient (deadzone)
c         if (deadzone) then
            alpha(iz)=C(ir,ith)%alphaturb
c         else
c            alpha(iz)=alphaturb ! no radial dependence on alpha, for now
c         endif
      enddo

      !  Print disk parameters ?
      !
      if (print_disk) then
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
            write(*,'("Tgas[",I0,"]= ",F6.1)') iz,tgas(iz) !constant before iter
            write(*,'("rhodust[",I0,"]= ",E16.5)') 1, f(iz,1)
            write(*,'("rhodust[",I0,"]= ",E16.5)') 2, f(iz,2)
            write(*,'("rhodust[",I0,"]= ",E16.5)') 3, f(iz,3)
            write(*,'("rhogas= ",E16.5)') rhogas(iz)
            write(*,'("  ")')
         enddo
      endif
      
      ! Make the gravitational force constant
      !
      omk    = sqrt(GG*mstar/radius**3d0)
      vkep   = sqrt(GG*mstar/radius)
      mugas=2.3

      !  Set up the time grid (equilibrium or time-dependent)
      !
      if (scseteq) then
         nts=10                 ! Always do a few timesteps
         allocate(timestep(0:nts))
         timestep(0)=0d0
         do its=1,nts
            timestep(its)=10d0**(its-nts) * ageoftheuniverse
         enddo
         settletime(1:ns)=ageoftheuniverse ! large nr to reach equilibrium
      else

         !  The time-dependent grid runs from a fraction of the orbital 
         !    period to the lifetime in logarhitmic steps  
         !
         settletime(1:ns)=lifetime !  Default upper limit

         !  First calculate orbital period (2 pi r / vkep = 2 pi / omk)
         period = 2d0 * pi / (omk * year)

         !  Time to travel from surface to midplane at maximum speed (z/r vkep)
         tsettmin = period / 4.
         tsettmin = tsettmin

         !  Size of the timegrid (log!), using 10 steps per order of magnitude
         !
         !  NOTE:  extra timesteps for:
         !         Rounding of decimal for negative nr (-1)
         !         Sample the interval [0, ..., tsettmin] (-11)
         !         tend > lifetime (+0)   ; already in nts i.o nts-1
         !         convergence (+3)
         !
         logtstart=  int(10d0*log10(tsettmin))-1
         logtend=    int(10d0*log10(lifetime))+3
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
            timestep(its) = ((10d0**(its/10d0)-1d0)/9d0)*10d0**(logtstart/10d0)
         enddo
         !
         !  Sample the interval [tsettmin,lifetime] +extra steps for convergence
         !
         do its=10,nts
            timestep(its) = 10d0**( (its-10+logtstart)   / 10d0)   ! *10**logtstart
         enddo
         !
         !  Stupidity checks Kees-style
         !
         if (timestep(10).le.timestep(9)) stop "grid doesnt connect"
         if (timestep(nts).le.log10(lifetime)) stop "grid too small"

         !  Save time grid to file
         if (print_time .and. ir.eq.100) then
            write(outputfile,'(a,"settle_timegrid.dat")') outdir(1:len_trim(outdir))
            open(unit=66,file=outputfile,RECL=1000)
            write(66,*) "# scset: time grid"
            write(66,*) "# It prints its,time"
            do its=0,nts
               write(66,88) its,timestep(its)
 88            format(I3,' ',D16.6)
            enddo
            close(unit=66)      !  Close the ouput file            
         endif

      endif ! finished setting up the time grid

      !
      !  When is model converged?
      !    From Kees code: 1d8
      !    Based on radial tau=1 convergence in TTs : 
      !    1d5:     too short
      !    1d6:     ok
      !    1d7-1d8  way too long
      !
      stoptolyear=1d6

c-------------------------------------------------------------------
c     Here starts the settling routine
c-------------------------------------------------------------------

      !  Now do the settling in a vertical slab, for each particle
      do is=1,ns

         !  Check if particle need to be settled
         !
         if(.not.do_settle(is)) then
!            print*,' dont settle ',is,' at r= ',radius/AU
            settletime(is)=0d0
            goto 222
         endif
         
         
         !  Compute Soundspeed and Stokes number
         !
         do iz=1,nz
         !  Compute the frictional cross section
         !
	         sigfr(iz)     = pi*agrain(is,iz)**2d0
            soundspeed2(iz) = kb*tgas(iz) / (mugas*amu)
            soundspeed(iz) = sqrt(soundspeed2(iz))
            !
            stokes(is,iz) = 0.75d0* mgrain(is,iz)/sigfr(iz) * 
     1           omk/(rhogas(iz)*soundspeed(iz)) *
     2           alpha(iz)**(2d0*qturb - 1d0)
            !
         enddo

         ! Compute the turbulent diffusion coefficient D 
         !   at the interface of the cells
         !
         ! Note: diffcoef is D * rhogas / (1 + St^2)
         ! D = alpha * soundspeed^2 / kepplerfrequency
         !
         do iz=2,nz

            diffcoef(is,iz)= 0.5d0*(rhogas(iz)  *alpha(iz)  *soundspeed2(iz)   
     1                              / omk / ( 1d0 + stokes(is,iz)**2d0 )                           )
     2           +           0.5d0*(rhogas(iz-1)*alpha(iz-1)*soundspeed2(iz-1)
     3                              / omk / ( 1d0 + stokes(is,iz-1)**2d0 )                         )
         enddo
         diffcoef(is,1)    = 0.d0  ! Flux at the inner boundary zero
         diffcoef(is,nz+1) = 0.d0  ! Flux at the outer boundary zero
 
         !  Calculate settling speeds at cell interfaces 
         !
         do iz=2,nz
 
            !  Interpolate temp and density of gas at the interface 
            !
            tempg        = 0.5d0 * ( tgas(iz-1) + tgas(iz) )
            rhog         = 0.5d0 * ( rhogas(iz-1) + rhogas(iz) )
            if(rhog.lt.1d-50) rhog=1d-50

            !  Compute the force on the grain at the interface
            !  NOTE: includes radiation pressure
            !
            force        = mgrain(is,iz)*omk*omk*zi(iz)*
     &                     (1d0+C(ir,(nz)-(iz-1))%FradZ*gas2dust)

            ! Compute the equilibrium settling velocity at the interface
            !
            vsett(is,iz) = -vdrift(tempg,rhog,mugas,force,sigfr(iz))

            ! NOTE: Warning: for very large particles this could be larger
            !       than (z/r)*v_Kepler, which is unphysical. The physics 
            !       here is that for such particles the equilibrium velocity 
            !       is never reached    !
            !       Trick: just limit it to (z/r)*v_kepler.
            !
            if (-vsett(is,iz).ge.vkep*(zi(iz)/radius)) then
                vsett(is,iz)=-vkep * (zi(iz)/radius)
            endif

            !  Check validity of this number 
            if(number_invalid(vsett(is,iz)).ne.0) then
               write(*,*) 'ERROR: NaN or Inf value detected (2)'
               stop 
            endif

         enddo ! Finished calculating settling speeds

         ! Settling speed at the boundaries
         vsett(is,1)    = 0.d0
         vsett(is,nz+1) = vsett(is,nz)

         !  Now start the settling calculation
         !  First copy the vertical density structure into a 1-D array
         !
         do iz=1,nz
            y(iz) = f(iz,is)
         enddo

         !  Use implicit integration
         eps  = 1d0
         eps1 = 1.d0-eps
         
         !  For checking convergence
         isteady=0

         !  Start the timeloop with the actual settling routine
         !
         do its=1,nts

            !  Save the previous vertical structure (for checking convergence)
            do iz=1,nz
               yold(iz) = y(iz)
            enddo

            !  Calculate dt for _this_ timestep
            dtime=(timestep(its)-timestep(its-1)) * year

            ! Fill the tridiagonal matrix elements a,b,c 
            !
            tr_a(1) = 0.d0
            do iz=2,nz
               tr_a(iz) = diffcoef(is,iz)/(dzc(iz)*dzi(iz)*rhogas(iz-1))
            enddo
            do iz=1,nz
               tr_b(iz) = -(1.d0/(dzc(iz)*rhogas(iz))) * 
     1              ( diffcoef(is,iz)/dzi(iz) + diffcoef(is,iz+1)/dzi(iz+1)  )
     2              + vsett(is,iz)/dzc(iz)
            enddo
            tr_c(nz) = 0.d0
            do iz=1,nz-1
               tr_c(iz) = diffcoef(is,iz+1)/(dzc(iz)*dzi(iz+1)*rhogas(iz+1))
     1           - vsett(is,iz+1)/dzc(iz)
            enddo

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

            ! Check for negative densities
            do iz=1,nz
               if(y(iz).lt.0.d0) y(iz)=0.d0
            enddo

            ! Check for zero density in bottom grid cell
            if(y(nz).lt.1d-50) y(nz)=1d-50

            !  Renormalize the dust density using cell volume (conserves mass)
            !  Using the surface density does not work in a spherical grid.
            !  The max difference in cell width is only 10% at z/r = 0.5
            !  so this should not affect the settling routine much
            !
            mcol=0d0
            mcol_old=0d0
            do iz=1,nz
               mcol=      mcol     + cellv(iz)*y(iz)
               mcol_old=  mcol_old + cellv(iz)*f(iz,is)
            enddo
            do iz=1,nz
               y(iz)=y(iz)*mcol_old/mcol
            enddo

            !  Calculate convergence criterium
            !
            difmax=0.d0
            do iz=1,nz
               dif = abs(yold(iz)-y(iz)) / abs(yold(iz)+y(iz)+1d-99)
               if(dif.gt.difmax) then
                  difmax = dif
               endif
            enddo

            !  Check for convergence (of one particle)
            !
            if(dtime/(difmax+1d-99).gt.stoptolyear*year) then
               isteady = isteady +1 
            endif

            !  Save convergence time and escape from loop
            !
            if(isteady.ge.3) then
               settletime(is) = timestep(its-3)
               goto 999
            endif

         enddo ! Reached settletime without converging
 999     continue ! Converged before settletime 

         !  Copy new density into to f
         !
         do iz=1,nz
            f(iz,is) = y(iz)
         enddo

 222     continue !  Arrive here if particle didn't need settling

      enddo !  Finish settling in a vertical slab, for each particle

c-------------------------------------------------------------------
c     Here ends the settling routine, write output and check densities
c-------------------------------------------------------------------

      !  Write convergence time at this distance to a file
      !  NOTE: This file+header is created at ir=0, and appended for other ir
      !
      if (scsetsave) then
         write(outputfile,'(a,"settle_eq.dat")') outdir(1:len_trim(outdir))

         !  Open a (new) file for writing
         if (ir.eq.1) then

            open(unit=66,file=outputfile,RECL=1000)

            !  Write header
            write(66,*) "# Convergence times for the settling routine"
            write(66,*) "# Per particle, made if scsetsave=.true."
            write(66,*) "# Number of: disk radii(nr), particles (ns)"
            write(66,99) "# ", D%nR-1,ns
 99         format(A2,' ',I6,' ',I3)
            write(66,*) "# It prints r, time to reach equilibrium (ns)"
         else
            open(unit=66,file=outputfile,RECL=1000,ACCESS='append')
         endif

         !  Write time to reach equilibrium
         !
         write(66,*) radius / AU, settletime(1:ns)  ! plot radius as well?
c         write(66,*) settletime(1:ns) 

         close(unit=66)         !  Close the ouput file
      endif

      !  Fix abundance at densities lower than 1d-45
      !  This is to avoid kinks in temperature contours when 
      !  densities approach the lower limit of 1d-50
      !  Useful for making plots only.
      !
      if (thinparticle.ge.-1.and.thinparticle.le.ngrains) then
         hitbottom=.false.
!         condition=ir.eq.120
         condition=0
         if (condition) write(*,*) "                                 ", radius / AU

         do iz=1,nz
            if (condition) write(*,'(e10.3,$)') sum(f(iz,1:ns))
            if (.not.hitbottom.and.sum(f(iz,1:ns)).le.1d-45) then !  arrive at <1d-45 ?
               hitbottom=.true.
               ibottom=max(iz-1,1)
               scalebottom=1.01d-45 / sum(f(iz-1,1:ns))
            endif
            if (hitbottom) then !  >1d-50 for everything above
               if (thinparticle.eq.-1) then ! use last density
                  f(iz,1:ns)=f(ibottom,1:ns) * scalebottom
               else if (thinparticle.eq.0) then ! use input abundance
                  f(iz,1:ns)=C(ir,nz)%w(1:ns)*1.01d-45
               else if (thinparticle.ge.1) then ! use particle # 
                  f(iz,1:ns)=1d-50
                  f(iz,thinparticle)=1.01d-45
               endif
            endif
         
            if (condition) write(*,'(e10.3)') sum(f(iz,1:ns))
         enddo
      endif
      
      !  Pass the new dust distribution back to MCMax
      !
      do iz=1,nz
         ith=(nz)-(iz-1)
         do is=1,ns
			if(do_settle(is)) then
	            rho(ith,is)=f(iz,is) ! ith: reverse order
			endif
         enddo
      enddo

      !  Deallocate arrays
      deallocate(timestep)


 333  continue                  ! arrive here if no settling at all
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



c-------------------------------------------------------------------------------
c This function returns the average radius and mass of grain ii at location i,j
c averaged over the iopac (i.e. the temperature dependend opacity)
c-------------------------------------------------------------------------------
	subroutine ComputeRM_nopac(ii,i,j,rv,mass)
	use Parameters
	IMPLICIT NONE
	real*8 rv,mass,tot
	integer ii,i,j,iopac
	
	rv=0d0
	do iopac=1,Grain(ii)%nopac
		rv=rv+C(i,j)%wopac(ii,iopac)*Grain(ii)%rv*Grain(ii)%rscale(iopac)
	enddo
	tot=sum(C(i,j)%wopac(ii,1:Grain(ii)%nopac))
	if(tot.gt.1d-100) then
		rv=rv/tot
	else
		rv=Grain(ii)%rv*Grain(ii)%rscale(iopac)
	endif
	
	mass=0d0
	do iopac=1,Grain(ii)%nopac
		mass=mass+Grain(ii)%rho(iopac)*C(i,j)%wopac(ii,iopac)*(Grain(ii)%rv*Grain(ii)%rscale(iopac))**3
	enddo
	mass=mass*4d0*pi/3d0

	return
	end

