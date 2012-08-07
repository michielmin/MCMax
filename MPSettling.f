c-----------------------------------------------------------------------
c This subroutine computes the scaling factor for the dust scaleheight 
c of different dust species with respect to the gas.
c Based on midplane (mp) temperatures and density, see Dubrulle et al 95
c-----------------------------------------------------------------------


      subroutine MPSettling()
      use Parameters
      implicit none
      integer i,ii,ir,ith,j
      real*8 gamma, omega, omega_r, tauf, tauf_csrho, cs_T, cs, mu
      real*8 rho_g, prefac, omegatau, sh_dtg_w, sh_dtg_d
      real*8, dimension(ngrains,1:D%nR-1) :: scale_w, scale_d
      !
      character*100 outputfile
      !
      doubleprecision, PARAMETER :: GG=6.6720000d-08
      doubleprecision, PARAMETER :: amu=1.66053886d-24
      doubleprecision, PARAMETER :: micron=1d-4
      !
      !  Some constants
      !
c      gamma= 5./3.
      gamma=2.                  ! Hasegawa & Pudritz 2009
      omega_r= sqrt(GG * D%Mstar)
      mu = 2.3
c      cs_T= sqrt(gamma * kb / (mu * amu))  ! gamma here???
      cs_T= sqrt(kb / (mu * amu))  ! gamma here???
      !
      prefac= (1. / (1.+ gamma)) ** (1./4.)
      !
      ith=D%nTheta-1            ! In the midplane
      !
      !  Loop over all grains
      !
      do ii=1,ngrains
         tauf_csrho= Grain(ii)%rho * Grain(ii)%rv
c         write(*,'("species:   ",I0)') ii
c         write(*,'("rgrain=    ",F6.1)') Grain(ii)%rv
c         write(*,'("rhograin=  ",F6.2)') Grain(ii)%rho
         !
         !  Loop over all radii
         !
         do i=1,D%nR-1
            !
            !  Gasdens, use dustdens*gas2dust if not set.
            !
            rho_g=C(i,ith)%gasdens*gas2dust
            if (rho_g.lt.(C(i,ith)%dens / 1d5 )) then
               rho_g=C(i,ith)%dens * gas2dust
            endif
            !
            !  Calculate Omega * Tau
            !
            cs= cs_t * sqrt(C(i,ith)%T)
            tauf= tauf_csrho / (cs * rho_g)
            omega = omega_r / (D%R_av(i) ** (3./2.) )
            omegatau=omega*tauf
            !
            !  Calculate scaling factor according to Dubrulle et al 1995
            !
            sh_dtg_d= prefac * sqrt(alphaturb / omegatau)
            scale_d(ii,i)=sh_dtg_d * (1+sh_dtg_d**2)**(-0.5)
            !
            !  Calculate scaling factor according to Weidenschilling 1977
            !
c            sh_dtg_w=sqrt(alphaturb) / omegatau
c            scale_w(ii,i)=sh_dtg_w * (1+sh_dtg_w**2)**(-0.5)
            !
            !  Pass scaling factor to MCMax
            !
            if(Grain(ii)%shtype.eq.'DISK') then
		         !
		         !  Tell MCMax we are using a scaling factor
		         !
		         Grain(ii)%settle=.true.
		         !
	            Grain(ii)%shscale(i)=scale_d(ii,i)
c            Grain(ii)%shscale(i)=scale_w(ii,i)
			endif
         enddo
      enddo
      !
      !  Write scaling factor(s) to a file
      !
      write(outputfile,'(a,"scalefactor.dat")') outdir(1:len_trim(outdir))
      open(unit=20,file=outputfile)
      do i=1,D%nR-1
c         write(20,*) D%R_av(i) / AU, omegatau(1:ngrains,i)
c         write(20,*) D%R_av(i) / AU, sh_dtg_d(1:ngrains,i)
c         write(20,*) D%R_av(i) / AU, scale_d(4,i), scale_d(10,i)
         write(20,*) D%R_av(i) / AU, scale_d(1:ngrains,i)
c         write(20,*) D%R_av(i) / AU, scale_d(2:ngrains,i),scale_w(2:ngrains,i)
      enddo
      close(unit=20)     
      !
      !  The end
      !
 666  return
      end
