c	module PAHModule
c	IMPLICIT NONE
c	real*8,allocatable :: BBQHP(:,:)
c	end module PAHModule

	subroutine PAHMCMax(niter)
	use Parameters
	use PAHModule
	IMPLICIT NONE
	integer i,j,ii,niter,itemp,ir
	real*8 tot,tau,temp0,temp1,bplanck
	
	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Computing multi photon PAH emissivity")')
	write(9,'("Computing multi photon PAH emissivity")')

c	open(unit=90,file='LRF.dat')
c	do i=1,nlam
c		write(90,*) lam(i),(C(j,1)%LRF(i),j=1,D%nR-1)
c	enddo
c	close(unit=90)

	temp0=10d0
	temp1=100000d0

	do itemp=1,NTQHP
		tgridqhp(itemp) = temp0 * (temp1/temp0)**
     %                      ((itemp-1.d0)/(NTQHP-1.d0))
	enddo

	if(.not.allocated(BBQHP)) then
		allocate(BBQHP(NTQHP,nlam))
		do itemp=1,NTQHP
			do i=1,nlam
				BBQHP(itemp,i)=bplanck(tgridqhp(itemp),nu(nlam-i+1))
			enddo
		enddo
	endif

	do i=1,D%nR-1
		call tellertje(i,D%nR-1)
		do j=1,D%nTheta-1
			tau=0d0
			if(C(i,j)%Ni.le.0.or.niter.eq.0) then
				do ir=0,i
					do ii=1,ngrains
						tau=tau+C(ir,j)%dens*(D%R(ir+1)-D%R(ir))*AU*C(ir,j)%w(ii)*Grain(ii)%Kpstar(1)
					enddo
				enddo
				if(niter.eq.0) tau=0d0
				C(i,j)%LRF(1:nlam)=D%Fstar(1:nlam)/(4d0*pi*D%R_av(i)**2)*exp(-tau)
			endif
			if(niter.ne.0.or.((tau.lt.10d0.and.tau.gt.1d0).or.j.eq.1)) then
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
						call PAH_Kees(nu,nlam,Grain(ii)%Nc
     &		,Grain(ii)%Kabs(1,1:nlam),Grain(ii)%Ksca(1,1:nlam)
     &		,C(i,j)%LRF(1:nlam),C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)
     &		,tgridqhp,C(i,j)%tdistr(Grain(ii)%qhpnr,1:NTQHP),NTQHP,temp0,temp1)
						call integrate(C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam),tot)
						if(tot.gt.1d-200) then
							C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)/tot
						else
							C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=0d0
						endif
 	 				endif
				enddo
			else if(tau.le.1d0.and.j.ne.1) then
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
						C(i,j)%tdistr(Grain(ii)%qhpnr,1:NTQHP)=C(i,1)%tdistr(Grain(ii)%qhpnr,1:NTQHP)
						C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(i,1)%QHP(Grain(ii)%qhpnr,1:nlam)
					endif
				enddo
			else
				do ii=1,ngrains
					if(Grain(ii)%qhp) then
						C(i,j)%tdistr(Grain(ii)%qhpnr,1:NTQHP)=0d0
						C(i,j)%QHP(Grain(ii)%qhpnr,1:nlam)=C(1,1)%QHP(Grain(ii)%qhpnr,1:nlam)
					endif
				enddo
			endif
		enddo
		do j=1,D%nTheta-1
			do ii=1,nqhp
				C(i,j)%Tqhp(ii)=0d0
				do itemp=1,NTQHP
					C(i,j)%Tqhp(ii)=C(i,j)%Tqhp(ii)+C(i,j)%tdistr(ii,itemp)*tgridqhp(itemp)**4
				enddo
				C(i,j)%Tqhp(ii)=C(i,j)%Tqhp(ii)**0.25d0
			enddo
		enddo
	enddo
	do j=1,D%nTheta-1
		C(0,j)%Tqhp(1:nqhp)=C(1,j)%Tqhp(1:nqhp)
	enddo


	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')


	deallocate(BBQHP)

	return
	end
	
	
	
	subroutine integrate2(spec,dnu,nlam,L)
	IMPLICIT NONE
	integer i,nlam
	real*8 spec(nlam),L,dnu(nlam)
	
	L=0d0
	do i=1,nlam-1
		L=L+(spec(i)+(spec(i+1)-spec(i))/2d0)*dnu(i)
	enddo

	return
	end




	subroutine PAH_Kees(nu,nlam,Ncarb,Kabs,Ksca,LRF,QHPemis,tempgrid2,tdistrib2,ntemp2,temp0,temp1)
#include "configure.h"
c     -----------------------------------------------------------------
c                 PROGRAM FOR COMPUTING PAH EMISSION 
c
c     This program computes PAH emission using the absorption/emission
c     cross-sections computed by the program of Ruud Visser. The 
c     emission is here computed by solving first the temperature
c     distribution function. After that the emission is calculated.
c     -----------------------------------------------------------------
      implicit none
	integer nlam
	real*8 nu(nlam),LRF(nlam),Kabs(nlam),Ksca(nlam),QHPemis(nlam)
	real*8 fracHigh,TQHPMax,Ncarb,tempgrid2(200),tdistrib2(200)
      integer iformat,itemp,inu,inu2,idum,it,ntemp2
      doubleprecision temp0,temp1,rho_pah,mp
      logical fex
      parameter(mp=1.6726d-24)
c     
#include "common_grid.h"
#include "common_dust.h"
#include "common_pah.h"
c
c     Default
c
      rho_pah = 1.d0
c
c     Read the input file
c
     	iformat=1
c      temp0=10d0
c      temp1=100000d0
      ntempgrid=ntemp2
      if(ntempgrid.gt.QUANTUM_NTEMP) then 
          write(*,*) 'ERROR pah_multi: ntempgrid out of range'
          stop
      endif
c
c     Read the frequency grid
c     
	freq_nr=nlam
      if(freq_nr.gt.FRSIZE_FREQ) then
          write(*,*) 'ERROR: frequency.inp has more gridpoints'
          write(*,*) '       than FRSIZE_FREQ. '
          write(*,*) '       Recompile with larger FRSIZE_FREQ'
          stop 13
      endif          
      do inu=1,freq_nr
          freq_nu(inu)=nu(freq_nr-inu+1)
      enddo
c
c     Read the exciting radiation field
c
      do inu=1,freq_nr
          radfield(inu)=LRF(freq_nr-inu+1)
      enddo

c     Read the dust opacities
c
	dust_nr_species=1
	dust_quantum(1) = 1
	dust_mgrain(1) = Ncarb*12.*mp
	dust_nr_temp(1)=1
      dust_setup_nrspecies = dust_nr_species
      dust_temprange_low(1,1)      = 0.d0
      dust_temprange_high(1,1)     = 1.d33
      do inu=1,freq_nr
          dust_kappawgt_abs(inu,1,1,1) = Kabs(freq_nr-inu+1)
          dust_kappawgt_scat(inu,1,1,1) = Ksca(freq_nr-inu+1)
      enddo
c
c     If the frequency weights have still not been specified (for 
c     whatever reason), then we generate them here
c
      call make_dust_freq_weights()
      dust_frwgt_read = 1
c
      dust_done_read = 1
	
c
c     Make a grid of temperatures
c
      do itemp=1,ntempgrid
c          tempgrid(itemp) = temp0 * (temp1/temp0)**
c     %                      ((itemp-1.d0)/(ntempgrid-1.d0))
          tempgrid(itemp) = tempgrid2(itemp)
      enddo
      dtempgrid(2) = 0.5 * ( tempgrid(3) - tempgrid(1) )
      do itemp=1,ntempgrid
          dtempgrid(itemp) = dtempgrid(2) * 
     %          ( tempgrid(itemp) / tempgrid(2) )
      enddo
c
c     Set up the heat content and peak temperatures
c
      call compute_peak_temperatures(1)
c
c     Compute cooling time curve
c
      call pah_cooling(1)
c
c     Solve the temperature distribution
c
      call solve_temp_distrib(1,radfield)
c
c     Compute the emission of PAH molecules
c
      call compute_emission(1)
c
c     Write emissivity
c
c      open(unit=1,file='emissivity.dat')
c      write(1,*) freq_nr
      do inu=1,freq_nr
		QHPemis(freq_nr-inu+1)=emissivity(inu)
c          write(1,*) freq_nu(inu),emissivity(inu)
      enddo
c      close(1)
c
c     Write peak temperatures to file
c
c      open(unit=1,file='peaktempmulti.dat')
c      write(1,*) freq_nr,ntempgrid
c      do itemp=1,ntempgrid
c          do inu=1,freq_nr
c              write(1,*) freq_nu(inu),peak_temp(itemp,inu)
c          enddo
c      enddo
c      close(1)
c
c     Write the cooling curve to file
c
c      open(unit=1,file='cooling.dat')
c      write(1,*) ntempgrid
c      do itemp=1,ntempgrid 
c          write(1,*) tempgrid(itemp),cooltime(itemp)
c      enddo
c      close(1)
c
c     Write the specific heat
c
c      open(unit=1,file='enthalpy.dat')
c      write(1,*) ntempgrid
c      do itemp=1,ntempgrid
c          write(1,*) heatcontent(itemp)
c      enddo
c      close(1)
c
c     Write the temperature distribution
c
c      open(unit=1,file='tempdistr.dat')
c      write(1,*) ntempgrid
      do itemp=1,ntempgrid
		tempgrid2(itemp)=tempgrid(itemp)
		tdistrib2(itemp)=tdistrib(itemp)
c          write(1,*) tempgrid(itemp),tdistrib(itemp),
c     %               tdistrib(itemp)/dtempgrid(itemp)
      enddo
c      close(1)
c
c     Write the matrix
c
c      open(unit=1,file='matrix.dat')
c      write(1,*) ntempgrid
c      do itemp=1,ntempgrid
c          do it=1,ntempgrid
c              write(1,*) mat(itemp,it)
c          enddo
c      enddo
c      close(1)
c
	return
      end


c     -----------------------------------------------------------------
c                 SOLVE FOR THE TEMPERATURE DISTRIBUTION
c     -----------------------------------------------------------------
      subroutine solve_temp_distrib(ispec,Jnu_excite)
      implicit none
c
      doubleprecision Jnu_excite(FRSIZE_FREQ),dummy
      doubleprecision matbk(QUANTUM_NTEMP)
      doubleprecision matwork(QUANTUM_NTEMP,QUANTUM_NTEMP)
      doubleprecision tdtmp(QUANTUM_NTEMP),result(QUANTUM_NTEMP)
      integer ispec,ievap,ntmax,itemp,it
      integer indx(QUANTUM_NTEMP)
      parameter(ntmax=QUANTUM_NTEMP)
c
#include "common_dust.h"
#include "common_grid.h"
#include "common_pah.h"
c
c     Fill the matrix
c
      call fill_matrix(ispec,Jnu_excite)
c
c     This matrix obviously has one zero eigenvalue because all the
c     temperature bins are converted into each other. This cannot be
c     inverted.
c
c     So replace last line with normalization
c     But remember the original line for later in this routine.
c
      do itemp=1,ntempgrid
          matbk(itemp)         = mat(ntempgrid,itemp)
          mat(ntempgrid,itemp) = 1.d0
      enddo
c
c     Put (0,0,...,0,1) as right-hand size
c
      do itemp=1,ntempgrid-1
          tdistrib(itemp) = 0.d0
      enddo
      tdistrib(ntempgrid) = 1.d0
c
c     Precondition this matrix 
c
      do itemp=1,ntempgrid-1
          dummy = 0.d0
          do it=1,ntempgrid
              if(abs(mat(itemp,it)).gt.dummy) dummy=abs(mat(itemp,it))
          enddo
          do it=1,ntempgrid
              mat(itemp,it) = mat(itemp,it) / dummy
          enddo          
      enddo
c
c     Invert this matrix
c
      dummy = 1.d0
      do itemp=1,ntempgrid
          do it=1,ntempgrid
              matwork(itemp,it) = mat(itemp,it)
          enddo
      enddo
      call ludcmp(matwork,ntempgrid,ntmax,indx,dummy)
      call lubksb(matwork,ntempgrid,ntmax,indx,tdistrib)
c
c     Solution has been obtained...
c
c     Check if it is indeed a solution
c
      dummy = 0.d0
      do itemp=1,ntempgrid
          tdtmp(itemp)  = 0.d0
          do it=1,ntempgrid
              tdtmp(itemp) = tdtmp(itemp) + 
     %                       mat(itemp,it) * tdistrib(it)
          enddo
          if(itemp.lt.ntempgrid) then
              dummy = dummy + abs(tdtmp(itemp))
          else
              dummy = dummy + abs(tdtmp(itemp)-1.d0)
          endif
      enddo
      if(dummy.gt.1d-10) then
          write(*,*) 'ERROR in matix inversion: Error = ',dummy
          stop
      endif
c     
c     Check if the solution makes any sense...
c
      dummy = 0.d0
      do itemp=1,ntempgrid
          dummy = dummy + tdistrib(itemp)
      enddo
      if(abs(dummy-1.d0).gt.1d-7) then 
          write(*,*) (tdistrib(itemp),itemp=1,ntempgrid)
          write(*,*) dummy
          stop 95664
      endif
c
c     Some values may be <0 due to numerical errors
c     Put to 1d-90 those values
c
      do itemp=1,ntempgrid
          if(tdistrib(itemp).le.0.d0) then
              tdistrib(itemp) = 1d-90
          endif
      enddo
c
c     Check again if the norm is not too much compromised
c
      dummy = 0.d0
      do itemp=1,ntempgrid
          dummy = dummy + tdistrib(itemp)
      enddo
      if(abs(dummy-1.d0).gt.1d-5) then 
          write(*,*) (tdistrib(itemp),itemp=1,ntempgrid)
          write(*,*) dummy
          stop 95665
      endif
c     
c     The next part of this routine is only for computing
c     the typical life time of the PAHs
c
c     Compute the typical life time of the grains
c     
      if(pahdes_temp(ispec).gt.0.d0) then
          write(*,*) '**********'
c
c         Find the index of the temperature grid 
c
          call hunt(tempgrid,ntempgrid,pahdes_temp(ispec),ievap)
          if(ievap.ge.ntempgrid) stop 90934
          if(ievap.lt.1) then
              write(*,*) pahdes_temp(ispec)
              stop 90935
          endif
c
c         Select all stuff below this temperature
c
          do itemp=1,ievap
              tdtmp(itemp) = tdistrib(itemp)
          enddo
          do itemp=ievap+1,ntempgrid
              tdtmp(itemp) = 0.d0
          enddo
c
c         Normalize
c
          dummy = 0.d0
          do itemp=1,ntempgrid
              dummy = dummy + tdtmp(itemp)
          enddo
          do itemp=1,ntempgrid
              tdtmp(itemp) = tdtmp(itemp) / dummy
          enddo
c
c         Restore last line of the matrix
c
          do itemp=1,ntempgrid
              mat(ntempgrid,itemp) = matbk(itemp)
          enddo
c
c         Do matrix multip
c
          do itemp=1,ntempgrid
              result(itemp) = 0.d0
              do it=1,ntempgrid
                  result(itemp) = result(itemp) + 
     %               mat(itemp,it) * tdtmp(it)
              enddo
          enddo
c
c         Check out how much has exceeded the tevap
c
          dummy = 0.d0
          do itemp=ievap+1,ntempgrid
              dummy = dummy + tdtmp(itemp)
          enddo
          pahdes_time(ispec) = 1.d0 / dummy
      else
          pahdes_time(ispec) = 1d99
      endif
c
      end


c     -----------------------------------------------------------------
c               FILL THE MATRIX, FOR A GIVEN RADIATION FIELD
c     -----------------------------------------------------------------
      subroutine fill_matrix(ispec,Jnu_excite)
      implicit none
c
      doubleprecision kappa(FRSIZE_FREQ),Jnu_excite(FRSIZE_FREQ)
      doubleprecision dnu(FRSIZE_FREQ)
      doubleprecision find_dust_kappa,epspeak,nphot,dum,pi,dt,hh
      integer inu,itemp,itpeak,ispec,it
      parameter(pi=3.1415926535897932385d0)
      parameter(hh=6.6262d-27)
c
#include "common_dust.h"
#include "common_grid.h"
#include "common_pah.h"
c
c     Make dnu
c
      dnu(1)       = 0.5 * ( freq_nu(2) - freq_nu(1) )
      dnu(freq_nr) = 0.5 * ( freq_nu(freq_nr) - freq_nu(freq_nr-1) )
      do inu=2,freq_nr-1
          dnu(inu) = 0.5 * ( freq_nu(inu+1) - freq_nu(inu-1) )
      enddo
c
c     Get the kappa
c
      do inu=1,freq_nr
          kappa(inu) = find_dust_kappa(inu,1,ispec,10d0,1,0)
      enddo
c
c     First clear the matrix
c
      do itemp=1,ntempgrid
          do it=1,ntempgrid
              mat(itemp,it) = 0.d0
          enddo
      enddo
c
c     Fill matrix for the cooling
c     
      do itemp=ntempgrid,2,-1
c
c         Compute the time it takes to cool from itemp to itemp-1
c
c         We do this in an explicit way, i.e. the time scale for
c         cooling from itemp to itemp-1 is computed using the 
c         infrared-emission cooling rate at itemp. In this way
c         we make sure that the shift of the material to lower 
c         temperature over an explicit time step delta t corresponds
c         to the cooling rate at the beginning of the time step
c         multiplied by delta t.
c
          dt = cooltime(itemp-1) - cooltime(itemp)
c
c         Compute the matrix contribution
c
          dum = 1.d0 / dt
c
c         Add this to the matrix, such that mass is conserved
c
          mat(itemp,itemp)   = mat(itemp,itemp)   - dum
          mat(itemp-1,itemp) = mat(itemp-1,itemp) + dum
      enddo
c
c     Fill the matrix for the heating
c
      do itemp=1,ntempgrid-1
c
c         Excite from temperature level itemp to the temperature
c         level given by the input photon energy. So make a loop
c         over all incoming photons
c
          do inu=1,freq_nr
c
c             Find where the energy is dumped
c
              itpeak  = peak_itemp(itemp,inu)
              epspeak = peak_eps(itemp,inu)
c
c             Find out how many photons are absorbed per second
c             by 1 gram of PAHs 
c
              nphot = 4 * pi * Jnu_excite(inu) * dnu(inu) * 
     %                kappa(inu) / ( hh * freq_nu(inu) )
c
c             Convert this into how many photons are absorbed per
c             second by 1 PAH molecule
c
              nphot = nphot * dust_mgrain(ispec)
c
c             This is the excitation rate...
c
c             Put this in the matrix
c
              if(itpeak.lt.ntempgrid) then
c
c                 Normal case
c
                  mat(itemp,itemp)    = mat(itemp,itemp) - nphot
                  mat(itpeak,itemp)   = mat(itpeak,itemp) +
     %                                  (1.d0-epspeak)*nphot
                  mat(itpeak+1,itemp) = mat(itpeak+1,itemp) +
     %                                  epspeak*nphot
              elseif(itpeak.eq.ntempgrid) then
c
c                 Excited to upper bin
c
                  mat(itemp,itemp)    = mat(itemp,itemp) - nphot
                  mat(itpeak,itemp)   = mat(itpeak,itemp) + nphot
              else
c
c                 Internal error...
c
                  stop 93971
              endif
          enddo
      enddo
c
      end


c     -----------------------------------------------------------------
c                           COMPUTE EMISSION
c     -----------------------------------------------------------------
      subroutine compute_emission(ispec)
	use PAHModule
      implicit none
c
      doubleprecision find_dust_kappa,absorb
      doubleprecision kappa(FRSIZE_FREQ),dnu(FRSIZE_FREQ)
      integer ispec,inu,itemp
c
#include "common_dust.h"
#include "common_grid.h"
#include "common_pah.h"
c
c     Make dnu
c
      dnu(1)       = 0.5 * ( freq_nu(2) - freq_nu(1) )
      dnu(freq_nr) = 0.5 * ( freq_nu(freq_nr) - freq_nu(freq_nr-1) )
      do inu=2,freq_nr-1
          dnu(inu) = 0.5 * ( freq_nu(inu+1) - freq_nu(inu-1) )
      enddo
c
c     Get the kappa
c
      do inu=1,freq_nr
          kappa(inu) = find_dust_kappa(inu,1,ispec,10d0,1,0)
      enddo
c
c     Reset emissivity
c
      do inu=1,freq_nr
          emissivity(inu) = 0.d0
      enddo
c
c     Add emissivity of each of the temperature points
c     Except for the itemp=1, because that is the zero-point
c     
      do itemp=2,ntempgrid
          do inu=1,freq_nr
              emissivity(inu) = emissivity(inu) + tdistrib(itemp) *
     %             kappa(inu)*BBQHP(itemp,inu)
          enddo
      enddo
c
c     Find the total energy output
c
      emistot = 0.d0
      do inu=1,freq_nr
          emistot = emistot + emissivity(inu) * dnu(inu)
      enddo
c
c     Check if this is not too far from the input energy
c
      absorb = 0.d0
      do inu=1,freq_nr
          absorb = absorb + kappa(inu)*radfield(inu)*dnu(inu)
      enddo
      if(abs(emistot/absorb-1.d0).gt.1d-5.and.absorb.ne.0d0) then
          write(*,*) 'WARNING: Error = ',abs(emistot/absorb-1.d0),
     %             emistot,absorb
      endif
c
      end


c     -----------------------------------------------------------------
c        COMPUTE PEAK TEMPERATURE FOR EXCITED PAH FROM SOME T-LEVEL 
c
c     This subroutine computes the peak temperature obtained by a
c     single photon with frequency nu, starting from a given temp level 
c     -----------------------------------------------------------------
      subroutine compute_peak_temperatures(ispec)
      implicit none
      doubleprecision hh,cv_pah_gram,energy
      doubleprecision enthalpy_pah_gram
      integer itemp,ispec,inu,itstart
      parameter(hh=6.6262d-27)
c
#include "common_dust.h"
#include "common_grid.h"
#include "common_pah.h"
c
c     Compute a table for the heat content
c
c     Now compute the rest
c     
      do itemp=1,ntempgrid
          heatcontent(itemp) = dust_mgrain(ispec) *
     %              enthalpy_pah_gram(tempgrid(itemp))
      enddo
c
c     First a loop over starting temperatures
c     
      do itstart=1,ntempgrid-1 
c
c         Now do loop over frequency, and compute the peak temperature
c         at every frequency (when a photon of that frequency hits
c         a PAH).
c         
          do inu=1,freq_nr
c         
c             Compute the energy of a single photon
c
              energy = hh * freq_nu(inu)
c
c             Find the index of the heat content corresponding
c             to that energy. This gives the peak temperature
c
              if(energy.lt.heatcontent(ntempgrid-1)-
     %                     heatcontent(itstart)) then
c
c                 Find the index of the temperature
c
                  call hunt(heatcontent,ntempgrid,
     %                      energy+heatcontent(itstart),itemp)
c
c                 Find the actual temperature through linear interpolation
c
                  peak_itemp(itstart,inu) = itemp
                  peak_eps(itstart,inu)   = 
     %                 ( energy + heatcontent(itstart) - 
     %                   heatcontent(itemp) ) / 
     %                 ( heatcontent(itemp+1) - heatcontent(itemp) )
                  if((peak_eps(itstart,inu).lt.0.d0).or.
     %                 (peak_eps(itstart,inu).gt.1.d0)) then
                      stop 82701
                  endif
                  peak_temp(itstart,inu)  = 
     %                 (1.d0-peak_eps(itstart,inu)) * tempgrid(itemp) +
     %                 peak_eps(itstart,inu) * tempgrid(itemp+1)
c
c                 If this temperature exceeds the evaporation temperature
c                 then put peak_itemp to ntempgrid as a signal
c
                  if((peak_temp(itstart,inu).gt.dust_tmax(ispec)).and.
     %                 (dust_tmax(ispec).ne.0.d0)) then
                      peak_itemp(itstart,inu) = ntempgrid
                  endif
              else
c                 write(*,*) 'ERROR: Temperature grid has not high'
c                 write(*,*) '       enough temperature to accomodate'
c                 write(*,*) '       photons of highest frequency...'
c                 stop
c                 
c                 If the temperature is larger than the temperature
c                 grid, then obviously it must have been evaporated.
c                 
                  peak_itemp(itstart,inu) = ntempgrid
              endif
          enddo
c
c         For user-information, store the peak temperature
c
          do inu=1,freq_nr
              peak_temp_store(itstart,inu,ispec) = 
     %               peak_temp(itstart,inu)
          enddo
c
      enddo
c
      end


c     -----------------------------------------------------------------
c                    COMPUTE THE COOLING TIME LINE
c
c     This routine returns the time the PAH reaches each of the 
c     temperatures of the given temperature grid, starting from the
c     highest temperature (i.e. the last temperature in the array),
c     and cooling down to the lowest.
c     For this, the routine uses dust opacity nr ispec.
c     -----------------------------------------------------------------
      subroutine pah_cooling(ispec)
	use PAHModule
      implicit none
c
      integer ispec
      doubleprecision temp,dtime,cv_pah_gram
      doubleprecision kappa(FRSIZE_FREQ),dum(FRSIZE_FREQ)
      doubleprecision enthalpy_pah_gram
      doubleprecision find_dust_kappa,integrate3,pi,mp
      parameter(pi=3.1415926535897932385d0)
      parameter(mp=1.6726d-24)
c     
      integer itemp,inu
c
#include "common_dust.h"
#include "common_grid.h"
#include "common_pah.h"
c
c     Tests
c
      if((ispec.lt.1).or.(ispec.gt.dust_nr_species)) then
          write(*,*) 'ERROR pah_cooling(): ispec out of range'
          stop
      endif
      if(tempgrid(ntempgrid).le.tempgrid(1)) then
          write(*,*) 'Tempeature grid must be in ascending order'
          stop
      endif
c
c     Get the kappa
c
      do inu=1,freq_nr
          kappa(inu) = find_dust_kappa(inu,1,ispec,10d0,1,0)
      enddo
c
c     Time at max temp = 0
c     
      cooltime(ntempgrid) = 0.d0
c
c     Now let's cool!
c
      do itemp=ntempgrid-1,1,-1 
c
c         Use the temperature of itemp+1, which is important for
c         the energy conservation of the multi-photon algorithm...
c
          temp = tempgrid(itemp+1)
c
c         Compute cooling rate per gram of dust
c
          do inu=1,freq_nr
              dum(inu) = BBQHP(itemp+1,inu) * kappa(inu)
          enddo
          qcool(itemp+1) = 4 * pi * integrate3(freq_nr,freq_nu,dum)
c
c         Compute the time it takes to cool...
c
          dtime = ( enthalpy_pah_gram(tempgrid(itemp+1)) - 
     %              enthalpy_pah_gram(tempgrid(itemp)) ) / 
     %            qcool(itemp+1)
c
c         Compute the new time
c
          cooltime(itemp) = cooltime(itemp+1) + dtime 
      enddo
c
c     For user-information, store the cooling time
c
      do itemp=1,ntempgrid
          cooltime_store(itemp,ispec) = cooltime(itemp)
      enddo
c
      end


cc     -----------------------------------------------------------------
cc                SPECIFIC HEAT FOR PAH MOLECULE
cc
cc     This is the specific heat per gram of PAH molecules in
cc     erg/gram/K.
cc     -----------------------------------------------------------------
c      function cv_pah_gram(temp)
c      implicit none
c      doubleprecision cv_pah_gram,temp,enthalpy_pah
c      doubleprecision eps
c      parameter(eps=1d-6)
cc
c      cv_pah_gram = ( enthalpy_pah(temp*(1.d0+eps))-
c     %              enthalpy_pah(temp*(1.d0)) ) / 
c     %              (temp*eps*12*mp)
cc
c      return
c      end


cc     -----------------------------------------------------------------
cc                      ENTHALPY FOR PAH MOLECULES
cc
cc     According to Siebenmorgen & Kruegel (1992) the enthalpy of 
cc     PAH molecules in erg/atom is (formula by Chase et al 1985):
cc     -----------------------------------------------------------------
c      function enthalpy_pah(temp)
c      implicit none
c      doubleprecision temp,enthalpy_pah
cc
cc     The formula of Chase
cc
c      enthalpy_pah = 4.15D-22 * (temp**3.3) / 
c     %               ( 1.D0 + 6.51D-3*temp + 1.5D-6*temp**2 +
c     %                 8.3D-7*(temp**2.3) )
cc
c      return
c      end

c     -----------------------------------------------------------------
c                   ENTHALPY FOR PAH MOLECULES PER GRAM 
c
c     According to Siebenmorgen & Kruegel (1992) the enthalpy of 
c     PAH molecules in erg/atom is (formula by Chase et al 1985):
c     -----------------------------------------------------------------
      function enthalpy_pah_gram(temp)
      implicit none
      doubleprecision temp,enthalpy_pah_gram,mp
      parameter(mp=1.6726d-24)
c
c     The formula of Chase
c
      enthalpy_pah_gram = 4.15D-22 * (temp**3.3) / 
     %                    ( 1.D0 + 6.51D-3*temp + 1.5D-6*temp**2 +
     %                    8.3D-7*(temp**2.3) )
c
c     Normalize to gram^(-1)
c
      enthalpy_pah_gram = enthalpy_pah_gram / (12*mp)      
c
      return
      end


c     --------------------------------------------------------------
c                THE BLACKBODY PLANCK FUNCTION B_nu(T)
c
c     This function computes the Blackbody function 
c
c                    2 h nu^3 / c^2
c        B_nu(T)  = -------------------    [ erg / cm^2 s ster Hz ]
c                   exp(-h nu / kT) - 1
c
c     ARGUMENTS:
c        nu    [Hz]            = Frequency
c        temp  [K]             = Electron temperature
c     --------------------------------------------------------------
      function bplanck(temp,nu)
      implicit none
      doubleprecision temp
      doubleprecision nu
      doubleprecision bplanck
c
      if(temp.eq.0.d0) then 
          bplanck = 0.d0
          return
      endif
c
      bplanck = 1.47455d-47 * nu * nu * nu / 
     %      (exp(4.7989d-11 * nu / temp)-1.d0) + 1.d-290
c
      return
      end

c     --------------------------------------------------------------
c                   A SIMPLE INTEGRATION ROUTINE
c     --------------------------------------------------------------
      function integrate3(n,x,f)
      implicit none
      integer n,i
      doubleprecision x(n),f(n),integrate3,int
      int=0.d0
      do i=2,n
          int=int+0.5d0*(f(i)+f(i-1))*(x(i)-x(i-1))
      enddo
      if(x(n).gt.x(1)) then
          integrate3 = int
      else
          integrate3 = -int
      endif
      return
      end


c     -----------------------------------------------------------------
c                     NUMERICAL RECIPES ROUTINE
c     -----------------------------------------------------------------
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      doubleprecision d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      doubleprecision aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

c     -----------------------------------------------------------------
c                     NUMERICAL RECIPES ROUTINE
c     -----------------------------------------------------------------
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      doubleprecision a(np,np),b(n)
      INTEGER i,ii,j,ll
      doubleprecision sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     --------------------------------------------------------------
c                      MAKE INDEXED FILENAME
c
c     This routine makes filenames of the type dusttemp_4.dat etc.
c     --------------------------------------------------------------
      subroutine make_indexed_filename(base,index,ext,filename)
      implicit none
      character*80 base
      integer index
      character*80 ext
      character*80 filename
      character*12 ch
c     
      if((index.lt.0).or.(index.ge.1000)) then
          write(*,*) 'ERROR in make_indexed_filename()'
          stop 729
      endif
      if(index.lt.10) then
          write(ch,11) index
 11       format(I1)
      elseif(index.lt.100) then
          write(ch,12) index
 12       format(I2)
      elseif(index.lt.1000) then
          write(ch,13) index
 13       format(I3)
      endif
      filename = base(1:len_trim(base))//
     %           ch(1:len_trim(ch))//
     %           ext(1:len_trim(ext))
c
      end



c     --------------------------------------------------------------
c                         FIND DUST OPACITY
c     
c     This function returns the dust opacity of a given species,
c     size, and at a given frequency-index. It takes consistently
c     into account the temperature range in which the dust grains
c     can exist.
c
c     --------------------------------------------------------------
      function find_dust_kappa(inu,isize,ispec,temp,iabs,iscat)
      implicit none
      doubleprecision temp,find_dust_kappa
      integer inu,isize,ispec,iabs,iscat
c
#include "common_grid.h"
#include "common_dust.h"
c
      doubleprecision tmin,tmax,dtmin,dtmax,plindex
      doubleprecision condense,fact,omfact
      integer ntemp,itlo,ithi
c
      condense = 1.d0
      ntemp = dust_nr_temp(ispec)
      tmin  = dust_temprange_low(1,ispec)
      tmax  = dust_temprange_high(ntemp,ispec)
      dtmin = dust_dtmin(ispec)
      dtmax = dust_dtmax(ispec)
      itlo  = 1
      ithi  = 1
      fact  = 1.d0
c
c     Now take care of temperature regimes
c
      if(temp.lt.tmin) then
          if(dtmin.eq.0.d0) then
c             
c             Hard switch-off at low temperature
c         
              condense = 0.d0
          else
c         
c             Smooth switch-off at low temperature
c         
              plindex  = log(1.d0-(dtmin/tmin))
              condense = (tmin/temp)**plindex
              if(condense.gt.1.d0) then
                  stop 9932
              endif
          endif
      elseif(temp.ge.tmax) then
          if((dtmax.eq.0.d0).and.(tmax.gt.0.d0)) then
c         
c             Hard switch-off at high temperature
c         
              condense = 0.d0
          else
c         
c             Smooth switch-off at high temperature
c         
              plindex  = log(1.d0+(dtmax/tmax))
              condense = (tmin/temp)**plindex
              if(condense.gt.1.d0) then
                  stop 9932
              endif              
          endif
      else
          if(ntemp.gt.1) then
              call hunt(dust_temprange_low(1,ispec),ntemp,temp,itlo)
              if(temp.gt.dust_temprange_high(itlo,ispec)) then
c
c                 Consistency check
c
                  if(itlo.ge.ntemp) then
                      write(*,*) 'BUG: Inconsistency in dust opacity'
                      write(*,*) '     temperature range....'
                      stop 13
                  endif
c
c                 We are in the linear interpolation regime
c
                  ithi = itlo + 1
                  fact = ( temp - dust_temprange_high(itlo,ispec) ) 
     %                   / ( dust_temprange_low(ithi,ispec) -
     %                       dust_temprange_high(itlo,ispec) )
              else
c
c                 We are within a range of constant opacity
c
                  ithi = itlo
                  fact = 1.d0
              endif
          endif
      endif
      omfact = 1.d0 - fact
c
c     Now get the opacities
c
      find_dust_kappa = 0.d0
      if(iabs.ne.0) then
          find_dust_kappa = find_dust_kappa +
     %           omfact * dust_kappawgt_abs(inu,itlo,isize,ispec) +
     %           fact * dust_kappawgt_abs(inu,ithi,isize,ispec)
      endif
      if(iscat.ne.0) then
          find_dust_kappa = find_dust_kappa +
     %           omfact * dust_kappawgt_scat(inu,itlo,isize,ispec) +
     %           fact * dust_kappawgt_scat(inu,ithi,isize,ispec)
      endif
c
c     If either iscat or iabs .lt.0 then do not incorporate the
c     evaporation of the grain. This can be handy to ensure that
c     at every temperature (even for the unphysically large ones)
c     a dust temperature can be computed. Simply for consistency.
c
      if((iabs.ge.0).and.(iscat.ge.0)) then
          find_dust_kappa = find_dust_kappa * condense
      endif
c
      return
      end


c     --------------------------------------------------------------
c                COMPUTE THE DUST FREQUENCY WEIGHTS
c
c     This routine computes the integration weights for the integrals
c     over the frequency. In principle this is a very general thing
c     and should not be here in this specialized module. But for now
c     this is the way we do it. Be sure to have generated the frequency
c     grid before calling this routine.
c
c     The frequency weight is such that
c                        __
c         /              \
c         | f(nu) dnu  =  >  f(nu_i) wgt_i
c         /              /_
c
c     We use the simple trapezium rule for the integral. 
c
c     --------------------------------------------------------------
      subroutine make_dust_freq_weights()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
c
      integer ifreq
c
      if(freq_nr.lt.1) then
          write(*,*) 'Zero or negative Number of freqs impossible'
          stop 99
      elseif(freq_nr.eq.1) then
          dust_freq_wgt(1) = 1.d0
cc          write(*,*) 'PROBLEM in dust module:'
cc          write(*,*) '  Attempt to construct frequency integration'
cc          write(*,*) '  weights for case of 1-frequency computation.'
cc          write(*,*) '  This is inconsistent. Are you sure you have '
cc          write(*,*) '  not accidently forgotten to switch off the '
cc          write(*,*) '  dust-temperature solver? This can be done '
cc          write(*,*) '  by setting iradproc_dust to -1 instead of 1.'
cc          stop 13
          return
      elseif((freq_nr.ge.2).and.(freq_nr.le.5)) then
          write(*,*) 'WARNING for dust temperature solver: '
          write(*,*) '  Too few frequency bins for consistently solving'
          write(*,*) '  dust temperatures! Continuuing while hoping '
          write(*,*) '  for the best....'
          dust_warn_few_freqs = 1
      else
          dust_freq_wgt(1) = 0.5d0 * abs( freq_nu(2) - freq_nu(1) )
          dust_freq_wgt(freq_nr) = 0.5d0 * abs( freq_nu(freq_nr) 
     %                              - freq_nu(freq_nr-1) )
          do ifreq=2,freq_nr-1
              dust_freq_wgt(ifreq) = 0.5d0 * 
     %             abs( freq_nu(ifreq+1) - freq_nu(ifreq-1) )
          enddo
          dust_frwgt_read = 1
      endif
c     
      end


