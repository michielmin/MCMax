c-----------------------------------------------------------------------
c This subroutine calculates the grainsize distribution according to a
c MRN-like powerlaw.
c It reads the grainsize from rgrain (Because Grain()%rv is not allocated),
c and sets the abundance in abun() (which is warg() in Init.f) 
c
c Input keywords:
c (mrn_rmin)    lower bound of the grain size grid.
c (mrn_rmax)    upper bound of the grain size grid.
c (mrn_index)   slope of the size distribution (def: 3.5)
c (mrn_ngrains) number of grains to include in calculation
c
c  Abundances of particles outside the sizegrid (mrn_rmin,mrn_rmax) are 
c  set to zero.
c
c-----------------------------------------------------------------------


      subroutine gsd_MRN(rgrain,abun)
      use Parameters
      implicit none

      doubleprecision    :: abun(1:ngrains),rgrain(1:ngrains),r_int(1:ngrains+1)
      doubleprecision    :: logr_int(1:ngrains+1)
      doubleprecision    :: rmin,rmax,r_index,r_factor,abun_norm
      integer            :: i,j,ii
      logical            :: diagnose

      doubleprecision, PARAMETER :: micron=1d-4 ! confusing! (also MPSet)

      !  Print additional diagnostics?
      diagnose=.false.

	  !  Set grainsize up to mrn_ngrains?
	  !
	  if (mrn_ngrains.eq.0)	then
	  	 mrn_ngrains=ngrains ! set default here to all grains (not done in Init.f)
	  endif
	  if (diagnose) write(*,*) "max grains (mrn / total): ",mrn_ngrains,ngrains
c      do ii=1,ngrains
c      write(*,'("abun",i02,"=",e16.8)') ii,abun(ii)
c	  enddo

      !  Find cell interfaces for grain size grid
      !
      do ii=2,mrn_ngrains  ! inbetween 
         logr_int(ii)=( log(rgrain(ii-1))+log(rgrain(ii)) )/2d0
      enddo
      logr_int(1)=         2*log(rgrain(1))       -logr_int(2)    ! lower bound
      logr_int(mrn_ngrains+1)= 2*log(rgrain(mrn_ngrains)) -logr_int(mrn_ngrains) ! upper bound
      !
      r_int(1:mrn_ngrains+1)= exp(logr_int(1:mrn_ngrains+1))
      !
      if (diagnose) then
         write(*,*)
         write(*,*) "Grid check: r_int(ii) < rgrain(ii) < r_int(ii+1)"
         write(*,*)
         do ii=1,mrn_ngrains 
            write(*,*) r_int(ii)," < ",rgrain(ii)," < ",r_int(ii+1)
         enddo
      endif
      
      !  Calculate mass per grain binsize, i.e. abundance
      !
      !  NOTE: MRN distribution is: f(a)=a^-3.5  , general a^-n
      !        in units of mass:    f(m)=m^-11/6 , general m^-(n+2)/3
      !
      !        For the abundance (or mass per size bin), integrate m f(m) dm
      !
      !        abun = 3/(4-n) m^( (4-n)/3 )  , between m(r_in) and m(r_out)
      !          OR   (check this)
      !        abun = 3/(4-n) r^(4-n)
      !
      !        note that it changes sign at n=4, where abun=log(r)  (ln, not 10log)
      !
      r_index=4d0-mrn_index
      r_factor=3d0/r_index  ! needed for sign change
c      if (diagnose) write(*,'("n= ",f4.1," , rindex= ",f4.1)') mrn_index,r_index
      !
      do ii=1,mrn_ngrains
         rmin=r_int(ii)
         rmax=r_int(ii+1)

         !  check if grainsize is between (mrn_rmin,mrn_rmax)
         if (rmin.ge.mrn_rmax .or. rmax.le.mrn_rmin) then
            abun(ii)=0d0  ! min or max outside grid cell
         else
            if (rmin.le.mrn_rmin) rmin=mrn_rmin ! min in gridcell
            if (rmax.ge.mrn_rmax) rmax=mrn_rmax ! max in gridcell
            if (rmax.le.rmin) stop 88354  ! shouldn't be
                                
            !  calculate abundance
            if (abs(mrn_index-4d0).le.1d-6) then 
               abun(ii)=log(rmax)-log(rmin)
            else
c               abun(ii)=rmin**r_index - rmax**r_index ! not normalized
               abun(ii)=r_factor*(rmax**r_index - rmin**r_index) ! not normalized
            endif
         endif

         ! print abundance
c         if (diagnose) write(*,'("abun",i02,"=",e16.8)') ii,abun(ii)
         
      enddo
      
      !  Calculate total mass & normalize
      !
      rmin=r_int(1)
      rmax=r_int(mrn_ngrains+1)
      if (rmin.ge.mrn_rmax .or. rmax.le.mrn_rmin) stop 16661 ! no grains between mrn_rmin,mrn_rmax
      if (rmin.le.mrn_rmin) rmin=mrn_rmin ! min in gridcell
      if (rmax.ge.mrn_rmax) rmax=mrn_rmax ! max in gridcell
                                
      !  calculate abundance
      if (abs(mrn_index-4d0).le.1d-6) then 
         abun_norm=log(rmax)-log(rmin)
      else
c     abun(ii)=rmin**r_index - rmax**r_index ! not normalized
         abun_norm=r_factor*(rmax**r_index - rmin**r_index) ! not normalized
      endif
      
      !  Normalize abundances
c      if (diagnose) write(*,'("total(abun)",e9.3,"=?=",e9.3," (sum)")') sum(abun(1:mrn_ngrains)),abun_norm
      abun(1:mrn_ngrains)=abun(1:mrn_ngrains) / abun_norm
      
      !  Check normalization
      !
      if (abs(sum(abun(1:mrn_ngrains))-1d0).ge.1d-6) then
         write(*,'("Something wrong with abundance normalization")')
         do ii=1,mrn_ngrains
            write(*,'(" abun",i02,"= ",f10.8)') ii,abun(ii)
         enddo
         write(*,'("Sum=",f10.8," =!= 1 -> stop 63636")') sum(abun(1:mrn_ngrains))
         stop 63636
      endif

	  !  Normalize including grains between mrn_ngrains+1 and ngrains ?
	  if (mrn_ngrains.lt.ngrains) then

         ! print abundance
         if (diagnose) write(*,*)
         do ii=1,ngrains
            if (diagnose) write(*,'("abun",i02,"= ",f10.8)') ii,abun(ii)
         enddo
         if (diagnose) write(*,*)	  
         if (diagnose) write(*,*) sum(abun)

	     abun_norm=1d0 - sum(abun(mrn_ngrains+1:ngrains))
	     abun(1:mrn_ngrains)=abun(1:mrn_ngrains) * abun_norm
	     
	  endif
            
      ! print abundance
      if (diagnose) write(*,*)
      do ii=1,ngrains
         if (diagnose) write(*,'("abun",i02,"= ",f10.8)') ii,abun(ii)
      enddo
      if (diagnose) write(*,*)
      if (diagnose) write(*,*) sum(abun)
      if (diagnose) write(*,*)
         
      !  Made it!
      !
      return
      end
