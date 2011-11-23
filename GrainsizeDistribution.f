c-----------------------------------------------------------------------
c This subroutine calculates the grainsize distribution according to a
c growth-fragmentation equilibrium. 
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
      doubleprecision    :: xi,T,alpha,sigma_g,sigma_d,rho_s,m_star,R,v_frag
      doubleprecision    :: a_01,a_12,a_l,a_p,a_r,a_sett
      doubleprecision    :: conc,cond,m_min,m_max
      integer            :: i
      !
      !  added
      !
      integer ir,ith,j,nr,ii,ngrains_gsd
      doubleprecision r_min,r_max,logstep,abun(1:ngrains),mdust
      doubleprecision surf0,surf(1:D%nR-1),tmp
      doubleprecision    :: r_int(1:ngrains+1) ! mass grid interface
      doubleprecision, PARAMETER :: micron=1d-4 ! confusing! (also MPSet)
      !
      !  Testing only
      !
      logical test1,test2,diagnose,oneiter,polycrash
      integer idiag
      !
      !  For plotting fit(r)
      !
      logical plot,plotallradii
      integer im
      integer ir1,ir2,ir3,ir4
      character*100 ffit,fgsd,freg
      doubleprecision    :: fit2D(1:D%nR-1,1:nm),abun2D(1:D%nR-1,1:ngrains)
      doubleprecision    :: regime(1:D%nR-1,1:6)
      !
      ! write simulation name
      !
      write(*,'("Finding grain size distribution")')
      write(9,'("Finding grain size distribution")')
      !
      !  Test1,test2, diagnose
      !
      test1=.false.
      test2=.false.
      oneiter=.true.
      oneiter=.false.
      if (gsd_diag .eq. 0) then
         polycrash=.false.
         diagnose=.false.
      else if (gsd_diag .gt. 0) then
         diagnose=.true.
         polycrash=.true.       ! (run after printing) -> replace by diagnose
         idiag=gsd_diag
      else
         polycrash=.false.
         diagnose=.false.
      endif
         
      !
      plot=.true.
      plotallradii=.true.
      plotallradii=.false.
      !
      ! SETUP test
      !
      if (.not.test1) goto 11                   ! skip test
      xi      = 1.8d0
      T       = 200d0
      alpha   = 1d-3
      sigma_g = 200d0
      sigma_d = 2d0
      rho_s   = 1.6d0
      m_star  = 0.5 * 1.989d33  ! 0.5 M_sun
      R       = 1d0 * 1.496d13  ! 1 AU
      v_frag  = 100d0
      m_min   = 6.7021d-15
      m_max   = 6.7021d3
      !
      ! produce mass grid
      !
      cond=1.d0/(1.d0-nm)*log10(m_min/m_max)
      conc=log10(m_min)-cond
      do i=1,nm
         m_grid(i)=10.d0**(conc+cond*i)
      end do
      !
      ! produce size grid
      !
      a_grid = (3d0*m_grid/(4d0*pi*rho_s))**(1d0/3d0)
      !
      ! call the fitting formula
      !
      write(*,*) 'test1 ...'
c      call fit_function(fit,nm,xi,T,alpha,sigma_g,sigma_d,
c     1     rho_s,m_grid,a_grid,m_star,R,v_frag)
      call fit_function(fit,a_01,a_12,a_l,a_p,a_r,a_sett,
     1     nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
      write(*,*) 'test1 finished'
 11   continue
      !
      !  Set parameters, make size grid
      !
      xi      = gsd_xi
      T       = 200d0
      alpha   = alphaturb
      sigma_g = 200d0
      sigma_d = 2d0
c      rho_s   = 1.6d0
      rho_s   = Grain(1)%rho
      m_star  = D%mstar
      R       = 1d0 * 1.496d13  ! 1 AU
      v_frag  = gsd_vfrag
      m_min   = 6.7021d-15
      m_max   = 6.7021d3
      !
      !   Make grid based on radius
      !   NOTE: lower limit MRN distribution = 0.025 micron
      !         at 0.01 micron, brownian motion reaches fragmentation velocity 
      !
      r_min=2.5d-2 * micron
!      r_min=1d-1 * micron
      r_max=1d5 * micron
      logstep=log10(r_max/r_min) / (nm - 1)
c      write(*,'("Grid from rmin= ",f9.2," to rmax= ",f9.2)') 
c     1     a_grid(1)/micron,a_grid(nm)/ micron
c      write(*,'("  r= ",f9.2)')  a_grid(1:nm)/ micron
c      write(*,'("Grid from mmin= ",e9.3," to mmax= ",e9.3)') 
c     1     m_grid(1),m_grid(nm)
c      write(*,'("  m= ",e9.3)')  m_grid(1:nm)
      do i=1,nm
         a_grid(i)=r_min*10d0**(logstep*(i-1))
         m_grid(i)=4d0/3d0*pi*rho_s*a_grid(i)**3d0
      enddo
c      write(*,'("Grid from rmin= ",f9.2," to rmax= ",f9.2)') 
c     1     a_grid(1)/micron,a_grid(nm)/ micron
c      write(*,'("  r= ",f9.2)')  a_grid(1:nm)/ micron      
c      write(*,'("Grid from mmin= ",e9.3," to mmax= ",e9.3)') 
c     1     m_grid(1),m_grid(nm)
c      write(*,'("  m= ",e9.3)')  m_grid(1:nm)
      !
      if (.not.test2) goto 22                   ! skip test 2
      write(*,*) 'test2 ...'
c      call fit_function(fit,nm,xi,T,alpha,sigma_g,sigma_d,
c     1     rho_s,m_grid,a_grid,m_star,R,v_frag)
      call fit_function(fit,a_01,a_12,a_l,a_p,a_r,a_sett,
     1     nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
      write(*,*) 'test2 finished'
 22   continue
      !
      !  Check if size grids comply
      !
      ngrains_gsd=gsd_rstep*int(log10(gsd_rmax/gsd_rmin))
      if (ngrains_gsd .ne. ngrains) then
         write(*,'("ngrains=",i02,", ngrains_gsd= ",i02)') ngrains_gsd,ngrains
         write(*,*) "grid doesn't match nr of particles -> stop 78987"
         stop 78987
      endif
      !
      !  Make and check cell interfaces grain size grid for MCMax
      !
      r_int(1)=gsd_rmin
      do ii=1,ngrains
         r_int(ii+1)=gsd_rmin*10**(1d0*ii/gsd_rstep)
c         write(*,'("grain",i02,": from ",3f10.2)') 
c     1        ii,r_int(ii)/micron,Grain(ii)%rv/micron,r_int(ii+1)/ micron
         !
         !  Check
         !
         if (r_int(ii).ge.Grain(ii)%rv.or.r_int(ii+1).le.Grain(ii)%rv) then 
            write(*,'("grain",i02," doesnt fit in the grid")') ii
            write(*,'("check that rgrain",i02," satisfies",
     1        f10.2," < ",f10.2," < ",f10.2)')
     2        ii,r_int(ii)/micron,Grain(ii)%rv/micron,r_int(ii+1)/ micron
         write(*,*) "grid doesn't match input grain size -> stop 67876"
         stop 67876
         endif
      enddo
      !
      !  Now for real
      !
      ith=D%nTheta-1            ! in midplane
      nr=D%nR-1
c$$$      !
c$$$      !  Check surfacedensity
c$$$      !  NOTE: Normalization is not enough to get surfacedens
c$$$      !
c$$$      mdust=0d0
c$$$      do ir=1,nr
c$$$         do j=1,D%nTheta
c$$$            mdust=mdust+C(ir,j)%mass
c$$$         enddo
c$$$      enddo
c$$$      write(*,*) mdust/Msun,D%Mtot/Msun
c$$$      write(*,*)
      !
      !  Calculate surfacedensity
      !  NOTE: Assuming a powerlaw
      !  TODO: check w/Michiel how to improve
      !
      tmp=2d0-D%denspow ! Two Minus P =!= temp
      surf0= tmp*(D%Mtot)/(2d0*pi*((D%Rout*AU)**(tmp)-(D%Rin*AU)**(tmp)))
      surf(1:nr)= surf0 * (D%R_av(1:nr))**(-1d0*D%denspow)
      !
      !  Normalize to dust mass
      !
      mdust=0d0
      do ir=1,nr
         mdust=mdust+(pi*(D%R(ir+1)**2-D%R(ir)**2)*AU**2) * surf(ir)
      enddo
      surf(1:nr)=surf(1:nr)/mdust*D%Mtot
c$$$      write(*,*) mdust/Msun,D%Mtot/Msun
c$$$      write(*,*)
      !
      !  Loop over all radii
      !
      do ir=1,nr
         !
         !  Always nice to have some dots to stare at
         !
         call tellertje(ir,D%nR-1)
         !
         !  Surfacedensity
         !  NOTE: Don't use C(ir,j)%mass here! Doesn't conserve mass and slope
         !
c$$$         sigma_d=0d0
c$$$         do j=1,D%nTheta
c$$$            sigma_d=sigma_d+C(ir,j)%mass
c$$$         enddo
c$$$         sigma_d=sigma_d / (pi*(D%R(ir+1)**2-D%R(ir)**2)*AU**2)
c$$$         sigma_g=sigma_d * gas2dust
c$$$c         write(*,*) sigma_g, sigma_d, gas2dust
c$$$c         write(*,*) D%R_av(ir) / AU, sigma_d, sigma_d / mdust * D%Mtot
         sigma_d=surf(ir)
         sigma_g=sigma_d * gas2dust
c$$$         write(*,*) D%R_av(ir) / AU, sigma_d, surf(ir)
         !
         !  Other radii-dependent parameters
         !
         R=D%R_av(ir)
         T=C(ir,ith)%T
         !
         !   Give it a shot...
         !
c         write(*,*) 'ir=',ir,' ...'
         fit(1:nm)=0d0
         if (.not.polycrash) then
c            call fit_function(fit,nm,xi,T,alpha,sigma_g,sigma_d,
c     1           rho_s,m_grid,a_grid,m_star,R,v_frag)
            call fit_function(fit,a_01,a_12,a_l,a_p,a_r,a_sett,
     1       nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
c         write(*,*) 'finished at ir=', ir
         endif
         !
         !   What does the fit look like
         !
         if (gsd_diag.ge.0) then
            write(*,'("ir= ",i3,", r= ",F8.3," sum= ",F8.5,
     1        " fit(0)= ",f10.7," fit(nm)= ",f10.7)') 
     2        ir,R / AU, sum(fit(1:nm))/sigma_d,fit(1)/sigma_d,fit(nm)/sigma_d
         endif
         !
         !  Is there a NaN?
         !
c$$$         if (IsNaN(sum(fit(1:nm)))) then
c$$$            write(*,'("Found a NaN at:")')
c$$$            write(*,'("   ir= ",i3,", r= ",F8.3," sum= ",F8.5,
c$$$     1           " fit(0)= ",f10.7," fit(nm)= ",f10.7)') 
c$$$     2           ir,R / AU, sum(fit(1:nm))/sigma_d,fit(1)/sigma_d,fit(nm)/sigma_d
c$$$            write(*,'("   Maybe increase gridsize? (nm) ")')
c$$$            write(*,'("  ")')
c$$$         endif
c         write(*,'("Regime, c=",i3,":",6F9.4)') ir+2, 
c     1        a_01/micron,a_12/micron,a_l/micron,
c     2        a_p/micron,a_r/micron,a_sett/micron
         if (IsNaN(sum(fit(1:nm)))) then
            write(*,'("Found a NaN at ir= ",i3,", r= ",F8.3,
     1              " sigma_d= ",F10.5, " temp= ",F6.1)') ir,R / AU, sigma_d, T
         endif
         !
         !  Diagnose
         !
         if (diagnose .and. ir .eq. idiag) then
            write(*,'("  ")')
            write(*,'("M= ",F6.3)') m_star / Msun
            write(*,'("gas2dust= ",F6.1)') gas2dust
            write(*,'("alpha= ",e10.3)')  alpha 
            write(*,'("  ")')      
            write(*,'("nm= ",I3)')  nm 
            write(*,'("Grid from rmin= ",f9.2," to rmax= ",f9.2)') a_grid(1)/micron,a_grid(nm)/ micron      
            write(*,'("Grid from mmin= ",e9.3," to mmax= ",e9.3)') m_grid(1),m_grid(nm)
            write(*,'("dust density= ",F5.2)') rho_s
            write(*,'("  ")')
            write(*,'("xi= ",F6.1)') xi
            write(*,'("vfrag= ",F6.1)') v_frag
            write(*,'("  ")')
            write(*,'("sigma dust= ",F10.5)') sigma_d
            write(*,'("temp= ",F6.1)') T
            write(*,'("  ")')
            if (polycrash) then
c               call fit_function(fit,nm,xi,T,alpha,sigma_g,sigma_d,
c     1              rho_s,m_grid,a_grid,m_star,R,v_frag)
               call fit_function(fit,a_01,a_12,a_l,a_p,a_r,a_sett,
     1       nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
            endif
c            stop 12321
c            write(*,'(" fit= ",e10.3)') fit
            write(*,'(" frac= ",f10.7)') fit / sigma_d
            write(*,'("  ")')
            write(*,'(" Made it through diagnostics -> stop 12343 ")')
            write(*,'("  ")')
            stop 12343
         endif
         !
         !  Convert fit on a_grid(1:nm) to abundances on rgrain(1:ngrains)
         !
         !  r_int(ngrains+1), a_grid(nm)
         do ii=1,ngrains
            !
            !  Outside grid -> abun=0d0
            !
            abun(ii)=0d0
            !
            !  Add everything in cell
            !
            do i=1,nm
               if (a_grid(i).ge.r_int(ii) .and. a_grid(i).lt.r_int(ii+1)) then
                  abun(ii)=abun(ii)+fit(i) / sigma_d
               endif
               !
               !  If there are no bigger particles, assume that contribution
               !  is in the biggest particle.
               !  NOTE: this might give funny errors in the gsd shape.
               !
               if (ii.eq.ngrains) then
                  if (a_grid(i) .ge. r_int(ii+1)) then
                     abun(ii)=abun(ii)+fit(i) / sigma_d
                  endif
               endif
            enddo
            !
            !  Change abundances in MCMax
            !
            do j=1,D%nTheta-1
               C(ir,j)%w(ii)=abun(ii)
            enddo
         enddo
         !
         !  check if they add up to 1
         !
         if (abs(sum(abun(1:ngrains))-1d0).ge.1d-6.and..not.polycrash) then
            write(*,'("Something wrong with abundance normalization")')
            do ii=1,ngrains
               write(*,'(" abun",i02,"= ",f10.8)') ii,abun(ii)
            enddo
            write(*,'("Sum=",f10.8," =!= 1 -> stop 63636")')sum(abun(1:ngrains))
            stop 63636
         endif
         !
         !  Store fit for plotting
         !
         fit2D(ir,1:nm)=fit(1:nm) / sigma_d
         abun2D(ir,1:ngrains)=abun(1:ngrains)
         regime(ir,1)=a_01
         regime(ir,2)=a_12
         regime(ir,3)=a_l
         regime(ir,4)=a_p
         regime(ir,5)=a_r
         regime(ir,6)=a_sett
      enddo
      !
      !  Write fit/gsd to file
      !
      if (plot) then
         !
         !  Prep dr
         !
         do im=1,nm   ! set dr
            if (im.ne.nm) then
               a_dr(im)=a_grid(im+1)-a_grid(im)
            else
               a_dr(nm)=a_dr(nm-1)
            endif
         enddo
         !
         !  First the analytical fit
         !
         write(ffit,'(a,"gsd_fit.dat")') outdir(1:len_trim(outdir))
         !
         open(unit=66,file=ffit,RECL=2000)
         !
         write(66,'("# Written from module GrainsizeDistribution.f ")')
         write(66,'("# It prints nr, nm ")')
         if (plotallradii) then
         write(66,'("# ",i4," ",i4)') nr,nm
         write(66,'("# Column info: columns [4:nr+4] are R[AU]")')
        write(66,'("# radius[um] dr[um]",6x,"mass [g]",4x,150(F12.6,x))') 
     1        D%R_av(1:nr) / AU
         do im=1,nm
            write(66, '(2(f12.3," "),e12.3," ",150(F12.10,x))') 
     1        a_grid(im)/micron, a_dr(im)/micron,
     2        m_grid(im),fit2D(1:nr,im)
         enddo
         else
            !
            !  Select radii 0.1,1.0,10.0,100.0 (not safe)
            !
            ir1=1
            do ir=2,nr        ! nr+1 doesn't exist
               if (D%R_av(ir-1)/AU.le.0.1)   ir1=ir
               if (D%R_av(ir-1)/AU.le.1.0)   ir2=ir
               if (D%R_av(ir-1)/AU.le.10.0)  ir3=ir
               if (D%R_av(ir-1)/AU.le.100.0) ir4=ir
            enddo
         write(66,'("# ",i4," ",i4)') 4,nm
         write(66,'("# Column info: columns [3:nr+3] are R[AU]")')
        write(66,'("# radius[um] dr[um]",6x,"mass [g]",4x,4(F12.6,x))') 
     1        D%R_av(ir1)/AU,D%R_av(ir2)/AU,D%R_av(ir3)/AU,D%R_av(ir4)/AU
         do im=1,nm
            write(66, '(2(f12.3," "),e12.3," ",4(F12.10,x))') 
     1           a_grid(im) / micron ,a_dr(im)/micron,m_grid(im),
     2           fit2D(ir1,im),fit2D(ir2,im),fit2D(ir3,im),fit2D(ir4,im)
         enddo
         endif
         !
         close(unit=66)
         !
         !  Then the gsd used in mcmax
         !
         write(fgsd,'(a,"gsd.dat")') outdir(1:len_trim(outdir))
         !
         open(unit=66,file=fgsd,RECL=2000)
         !
         write(66,'("# Written from module GrainsizeDistribution.f ")')
         write(66,'("# It prints nr, ngrains ")')
         if (plotallradii) then
            write(66,'("# ",i4," ",i4)') nr,ngrains
            write(66,'("# Column info: columns [4:nr+4] are R[AU]")')
        write(66,'("# radius[um] dr[um]",6x,"mass [g]",4x,150(F12.6,x))') 
     1           D%R_av(1:nr) / AU
            do ii=1,ngrains
               write(66, '(2(f12.3," "),e12.3," ",150(F12.10,x))') 
     1             Grain(ii)%rv / micron,(r_int(ii+1)-r_int(ii))/micron,
     2             Grain(ii)%m,abun2D(1:nr,ii)
            enddo
         else
            !
            !  Select radii 0.1,1.0,10.0,100.0 (not safe)
            !
            ir1=1
            do ir=2,nr        ! nr+1 doesn't exist
               if (D%R_av(ir-1)/AU.le.0.1)   ir1=ir
               if (D%R_av(ir-1)/AU.le.1.0)   ir2=ir
               if (D%R_av(ir-1)/AU.le.10.0)  ir3=ir
               if (D%R_av(ir-1)/AU.le.100.0) ir4=ir
            enddo
c            write(*,*) ir1,ir2,ir3,ir4
c            write(*,*) D%R_av(ir1)/AU,D%R_av(ir2)/AU,D%R_av(ir3)/AU,D%R_av(ir4)/AU
         write(66,'("# ",i4," ",i4)') 4,ngrains
         write(66,'("# Column info: columns [3:nr+3] are R[AU]")')
        write(66,'("# radius[um] dr[um]",6x,"mass [g]",4x,4(F12.6,x))') 
     1        D%R_av(ir1)/AU,D%R_av(ir2)/AU,D%R_av(ir3)/AU,D%R_av(ir4)/AU
         do ii=1,ngrains
            write(66, '(2(f12.3," "),e12.3," ",4(F12.10,x))') 
     1           Grain(ii)%rv / micron ,(r_int(ii+1)-r_int(ii))/micron,
     2           Grain(ii)%m,
     3       abun2D(ir1,ii),abun2D(ir2,ii),abun2D(ir3,ii),abun2D(ir4,ii)
         enddo
         endif
         !
         close(unit=66)
         !
         !  Also store the regime changes
         !
         write(freg,'(a,"gsd_regime.dat")') outdir(1:len_trim(outdir))
         !
         open(unit=66,file=freg,RECL=2000)
         !
         write(66,'("# Written from module GrainsizeDistribution.f ")')
         write(66,'("# It prints nr, nreg ")')
         write(66,'("# ",i4," ",i4)') nr,6
         write(66,'("# radius[AU]",6(7x,a6))') 
     1        "  a_01 [mu]","  a_12 [mu]","   a_l [mu]",
     2        "   a_p [mu]","   a_r [mu]","a_sett [mu]"
         if (plotallradii) then
            do ir=1,nr
               write(66,'(f12.6,6(" ",f12.3))') 
     1              D%R_av(ir)/AU,regime(ir,1:6)/micron
            enddo
         else
            !
            !  Select radii 0.1,1.0,10.0,100.0 (not safe)
            !
            ir1=1
            do ir=2,nr        ! nr+1 doesn't exist
               if (D%R_av(ir-1)/AU.le.0.1)   ir1=ir
               if (D%R_av(ir-1)/AU.le.1.0)   ir2=ir
               if (D%R_av(ir-1)/AU.le.10.0)  ir3=ir
               if (D%R_av(ir-1)/AU.le.100.0) ir4=ir
            enddo
            write(66,'(f12.6,6(" ",f12.3))') 
     1           D%R_av(ir1)/AU,regime(ir1,1:6)/micron
            write(66,'(f12.6,6(" ",f12.3))') 
     1           D%R_av(ir2)/AU,regime(ir2,1:6)/micron
            write(66,'(f12.6,6(" ",f12.3))') 
     1           D%R_av(ir3)/AU,regime(ir3,1:6)/micron
            write(66,'(f12.6,6(" ",f12.3))') 
     1           D%R_av(ir4)/AU,regime(ir4,1:6)/micron
         endif
         !
         close(unit=66)
      endif
      !
      if (oneiter) then
         write(*,'("  ")')
         write(*,'(" Made it through all radii -> stop 12345 ")')
         write(*,'("  ")')
         stop 12345
      endif
      !
      !  Made it!
      !
      return
      end
