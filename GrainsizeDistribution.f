c-----------------------------------------------------------------------
c This subroutine calculates the grainsize distribution according to a
c growth-fragmentation equilibrium in a vertical slab.
c The equilibrium value is dominated by collisions in the midplane, so 
c so the detailed vertical structure is not taken into account (similar to
c mpset=.true. and mpstr=.true.).
c
c Input:
c (gsd_rmin)   lower bound of the grain size grid.
c (gsd_rmax)   upper bound of the grain size grid.
c (gsd_rstep)  number of grain sizes per order of magnitude.
c              NOTE: The grains sizes supplied to MCMax (rgrain) should 
c                    fall within this grid.
c
c (gsd_vfrag)  velocity at which grains fragment. Increase for bigger grains
c (gsd_xi)     slope of the size distribution of grain fragments.
c              NOTE: has not yet been tested for xi =!= 1.8
c
c Output:
c "gsd.dat"        Radial grain size distribution for MCMax (Mass per bin)
c "gsd_fit.dat"    Radial grain size distribution from Tills fit (Mass per bin)
c "gsd_regime.dat" Regime sizes used in construct the fit
c
c Credits: Tilman Birnstiel (birnstiel@mpia.de)
c          Gijs Mulders     (mulders@uva.nl)
c
c-----------------------------------------------------------------------


      subroutine GrainsizeDistribution()
      use Parameters
      use fit_module
      implicit none
      integer, parameter :: nm=250 ! ~40 per order of magnitude (6)
      doubleprecision    :: fit(1:nm),m_grid(1:nm),a_grid(1:nm),a_dr(1:nm)
      doubleprecision    :: xi,T,alpha(1:D%nR-1),sigma_g,sigma_d,rho_s,m_star,R,v_frag
      doubleprecision    :: a_01,a_12,a_l,a_p,a_r,a_sett
      doubleprecision    :: conc,cond,m_min,m_max
      integer            :: i

      !  Extra variables
      integer ir,ith,j,nr,ii
      doubleprecision logstep,abun(1:ngrains),r_int(1:ngrains+1)
      doubleprecision surf0,surf(1:D%nR-1),tmp
      doubleprecision, PARAMETER :: micron=1d-4 ! confusing! (also MPSet)
      doubleprecision shdust,shgas,col,dcostheta(1:D%nTheta-1)

      !  Tests for diagnosing bugs
      logical runonce,printgrid

      !  For plotting fit(r)
      integer im
      character*100 ffit,fgsd,freg
      doubleprecision    :: fit2D(1:D%nR-1,1:nm),abun2D(1:D%nR-1,1:ngrains)
      doubleprecision    :: regime(1:D%nR-1,1:6)

      ! write simulation name
      write(*,'("Finding grain size distribution")')
      write(9,'("Finding grain size distribution")')

c-------------------------------------------------------------------
c     Here starts the initialisation
c-------------------------------------------------------------------

      !  Run tests / diagnostics?
      printgrid=.false.
      runonce=.false.
c     gsd_diag=50               ! can be set through command line      
c     gsd_plot=.true.           ! make plots of detailed gsd?

      !  
      ith=D%nTheta-1            ! in midplane
      nr=D%nR-1

      !  Fixed input parameters
      !
      xi      = gsd_xi
      rho_s   = Grain(1)%rho
      m_star  = D%mstar
      v_frag  = gsd_vfrag

      !  Deadzone: average turbulence over 1st scale height
      if(deadzone) then

         do i=1,nr
            call calcscaleheight(i,shgas,shdust)
            col=0d0
            alpha(i)=0d0

            do j=D%nTheta-1,1,-1
               dcostheta(j)=(D%Theta(j)-D%Theta(j+1)) 
               if(D%R_av(i)*cos(D%theta_av(j))/AU.lt.shgas) then
                  col=col+ C(i,j)%gasdens*gas2dust*(D%R_av(i)*dcostheta(j))
                  alpha(i)=alpha(i) + C(i,j)%gasdens*gas2dust*(D%R_av(i)*dcostheta(j))*C(i,j)%alphaturb

!                  col=col+ C(i,j)%mass
!                  alpha(i)=alpha(i) + C(i,j)%mass*C(i,j)%alphaturb
               endif
            enddo
            alpha(i)=alpha(i)/col
            print*,i,D%R_av(i),alpha(i)
         enddo
      else
         alpha(1:nr) = alphaturb         
      endif
      

      !   Grain size grid for the simulation, based on radius
      !   NOTE: lower limit MRN distribution = 0.025 micron
      !         at 0.01 micron, brownian motion reaches fragmentation velocity 
      !
      if (gsd_rmin.lt.2.5d-6) then
         write(*,'(/,"Minimum grain size below fragmentation limit,")')
         write(9,'(/,"Minimum grain size below fragmentation limit,")')
         write(*,'("  forcing gsd_rmin=0.025 micron",/)')
         write(9,'("  forcing gsd_rmin=0.025 micron",/)')
         gsd_rmin=2.5d-6
      endif

      !  the actual grid
      logstep=log10(gsd_rmax/gsd_rmin) / (nm - 1)
      do i=1,nm
         a_grid(i)=gsd_rmin*10d0**(logstep*(i-1))
         m_grid(i)=4d0/3d0*pi*rho_s*a_grid(i)**3d0
         a_dr(i)=gsd_rmin*(10d0**(logstep*(i-0.5))-10d0**(logstep*(i-1.5)))
      enddo

      !  check if grids comply (within 1d-6). If not, issue a warning
      if (Grain(1)%rvmin.gt.a_grid(1)*(1d0+1d-6)) then
         write(*,'("Grid smaller than smallest particle",f,f)') Grain(1)%rvmin,a_grid(1)
         write(9,'("Grid smaller than smallest particle",f,f)') Grain(1)%rvmin,a_grid(1)
      endif
      if (Grain(ngrains)%rvmax.lt.a_grid(nm)*(1d0-1d-6)) then
         write(*,'("Grid bigger than largest particle",f,f)') Grain(ngrains)%rvmax,a_grid(nm)
         write(9,'("Grid bigger than largest particle",f,f)') Grain(ngrains)%rvmax,a_grid(nm)
      endif

      ! print grid (debugging only)
      if (printgrid) then
         write(*,'(/,"Simulation from rmin= ",f9.2," to rmax= ",f9.2)') 
     1        a_grid(1)/micron,a_grid(nm)/ micron
         write(*,'("Simulation from mmin= ",e9.3," to mmax= ",e9.3,/)') 
     1        m_grid(1),m_grid(nm)
         do i=1,nm
            write(*,'(" i=",i3," r= ",f10.3,"  m= ",e9.3)') 
     1           i,a_grid(i)/ micron, m_grid(i)
         enddo
      endif

      !  Grain size grid used in MCMax
      !  NOTE: Lowest and highest grid cell are extended to encompass
      !        the entire gsd-grid 
      !
      r_int(1)        =min(Grain(1)%rvmin,a_grid(1))
      r_int(2:ngrains)=Grain(2:ngrains)%rvmin
      r_int(ngrains+1)=max(Grain(ngrains+1)%rvmax,a_grid(nm))
      

      !  Set up the surface density
      !
      do i=1,nr
         surf0=0d0              ! 0.5 mass between R(i) and R(i+1)
         do j=1,D%nTheta
            surf0=surf0+C(i,j)%mass
         enddo
         surf(i)=surf0/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)
      enddo

c-------------------------------------------------------------------
c     Here starts the loop that calculates the distributions
c-------------------------------------------------------------------

      !  Loop over all radii
      do ir=1,nr

         !  Always nice to have some dots to stare at
         if (gsd_diag.eq.-1) call tellertje(ir,D%nR-1)

         !  Surfacedensity, radius, temperature
         sigma_d=surf(ir)
         sigma_g=sigma_d * gas2dust
         R=D%R_av(ir)
         T=C(ir,ith)%T

         !  Run Til's code (default gsd_diag=-1)
         fit(1:nm)=0d0
         if (gsd_diag.le.0) then
            call fit_function(fit,a_01,a_12,a_l,a_p,a_r,a_sett,
     1       nm,xi,T,alpha(ir),sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
         endif

         !  Check for NaN values (shouldn't occur)
         !
         if (IsNaN(sum(fit(1:nm))).and.gsd_diag.le.0) then
            write(*,'("Found a NaN at ir= ",i3,", r= ",F8.3,
     1              " sigma_d= ",F10.5, " temp= ",F6.1)') ir,R / AU, sigma_d, T
            write(*,'("Regime, c=",i3,":",6F9.4)') ir+2, 
     1           a_01/micron,a_12/micron,a_l/micron,
     2           a_p/micron,a_r/micron,a_sett/micron
         endif

         !  Run some diagnostics (Debugging part)
         !
         if (gsd_diag.eq.0) then
            write(*,'("ir= ",i3,", r= ",F8.3," sum= ",F8.5,
     1        " fit(0)= ",f10.7," fit(nm)= ",f10.7)') 
     2        ir,R / AU, sum(fit(1:nm))/sigma_d,fit(1)/sigma_d,fit(nm)/sigma_d
         endif
         
         !  Run the code only once at ir=gsd_diag
         !
         if (ir .eq. gsd_diag) then
            write(*,'(/,"  M= ",F6.3)') m_star / Msun
            write(*,'("  gas2dust= ",F6.1)') gas2dust
            write(*,'("  alpha= ",e10.3)')  alpha(ir) 
            write(*,'(/,"  nm= ",I3)')  nm 
            write(*,'("  Grid from rmin= ",f10.3," to rmax= ",f10.3)') a_grid(1)/micron,a_grid(nm)/micron
            write(*,'("  Grid from mmin= ",e10.3," to mmax= ",e10.3)') m_grid(1),m_grid(nm)
            write(*,'("  dust density= ",F5.2)') rho_s
            write(*,'(/,"  xi= ",F6.1)') xi
            write(*,'("  vfrag= ",F6.1)') v_frag
            write(*,'(/,"  sigma dust= ",F10.5)') sigma_d
            write(*,'("  temp= ",F6.1,/)') T
            call fit_function(fit,a_01,a_12,a_l,a_p,a_r,a_sett,
     1           nm,xi,T,alpha(ir),sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
            write(*,'(/,"  im,  mass fraction= ",/)')
            do i=1,nm
               write(*,'(x,i4,x,f10.7)') i,fit(i) / sigma_d
            enddo
            write(*,'(/,"Made it through diagnostics -> stop 12343 ",/)')
            stop 12343
         endif

         !  Convert result to the MCMax grid
         !    fit(1:nm)     ->  abun(1:ngrains)    
         !    a_grid(1:nm)  ->  r_int(1:ngrains+1)
         !    mass per bin  ->  abundance per bin
         !
         abun(1:ngrains)=0d0
         do ii=1,ngrains

            ! calculate abundance
            do i=1,nm
               if (a_grid(i).ge.r_int(ii).and.a_grid(i).lt.r_int(ii+1)) then
                  abun(ii)=abun(ii)+fit(i)/sigma_d
               endif
            enddo

            !  check and store (MCMax crashes at low, non-zero abundances)
            do j=1,D%nTheta-1
               if (abun(ii) .le.1d-100) abun(ii)=0d0
               C(ir,j)%w(ii)=abun(ii)
               C(ir,j)%w0(ii)=abun(ii)
            enddo

         enddo
         
         !  check if they add up to 1
         if (abs(sum(abun(1:ngrains))-1d0).ge.1d-6.and.gsd_diag.le.0) then
            write(*,'("Something wrong with abundance normalization")')
            do ii=1,ngrains
               write(*,'(" abun",i02,"= ",f10.8)') ii,abun(ii)
            enddo
            write(*,'("Sum=",f10.8," =!= 1 -> stop 63636")')sum(abun(1:ngrains))
            stop 63636
         endif

         !  Store fit for plotting
         fit2D(ir,1:nm)=fit(1:nm) / sigma_d
         abun2D(ir,1:ngrains)=abun(1:ngrains)
         regime(ir,1)=a_01
         regime(ir,2)=a_12
         regime(ir,3)=a_l
         regime(ir,4)=a_p
         regime(ir,5)=a_r
         regime(ir,6)=a_sett
      enddo

c-------------------------------------------------------------------
c     Here starts the writing of output files
c-------------------------------------------------------------------

      !  First the gsd used in mcmax (abundance versus radius)
      !
      write(fgsd,'(a,"gsd.dat")') outdir(1:len_trim(outdir))
      open(unit=66,file=fgsd,RECL=2000)
      
      !  Wite header
      write(66,'("# Written from module GrainsizeDistribution.f ")')
      write(66,'("# It prints nr, ngrains ")')
      write(66,'("# ",i4," ",i4)') nr,ngrains
      write(66,'("# Column info: columns [4:nr+4] are R[AU]")')
      write(66,'("# radius[um] dr[um]",6x,"mass [g]",4x,150(F12.6,x))') 
     1     D%R_av(1:nr) / AU
      do ii=1,ngrains
         write(66, '(2(f12.3," "),e12.3," ",150(F12.10,x))') 
     1        Grain(ii)%rv / micron,(Grain(ii)%rvmax-Grain(ii)%rvmin)/micron,
     2        Grain(ii)%rho * 4d0/3d0*pi*Grain(ii)%rv**3d0,abun2D(1:nr,ii)
      enddo
      close(unit=66)

      !  Also store the regime changes
      !
      write(freg,'(a,"gsd_regime.dat")') outdir(1:len_trim(outdir))
      open(unit=66,file=freg,RECL=2000)

      !  Write header
      write(66,'("# Written from module GrainsizeDistribution.f ")')
      write(66,'("# It prints nr, nreg ")')
      write(66,'("# ",i4," ",i4)') nr,6
      write(66,'("# radius[AU]",6(7x,a6))') 
     1     "  a_01 [mu]","  a_12 [mu]","   a_l [mu]",
     2     "   a_p [mu]","   a_r [mu]","a_sett [mu]"
      do ir=1,nr
         write(66,'(f12.6,6(" ",f12.3))') D%R_av(ir)/AU,regime(ir,1:6)/micron
      enddo
      close(unit=66)

      !  Also plot the analytical fits (optional)
      !
      if (gsd_plot) then

         !  open file
         write(ffit,'(a,"gsd_fit.dat")') outdir(1:len_trim(outdir))
         open(unit=66,file=ffit,RECL=2000)

         !  Write header
         write(66,'("# Written from module GrainsizeDistribution.f ")')
         write(66,'("# It prints nr, nm ")')
         write(66,'("# ",i4," ",i4)') nr,nm
         write(66,'("# Column info: columns [4:nr+4] are R[AU]")')
        write(66,'("# radius[um] dr[um]",6x,"mass [g]",4x,150(F12.6,x))') 
     1        D%R_av(1:nr) / AU
         do im=1,nm
            write(66, '(2(f12.3," "),e12.3," ",150(F12.10,x))') 
     1           a_grid(im)/micron, a_dr(im)/micron,
     2           m_grid(im),fit2D(1:nr,im)
         enddo
         close(unit=66)

      endif
      
c-------------------------------------------------------------------
c     Made it!
c-------------------------------------------------------------------

      !  Debugging
      if (runonce) then
         write(*,'(/," Made it through all radii -> stop 12345 ",/)')
         stop 12345
      endif

      return
      end
