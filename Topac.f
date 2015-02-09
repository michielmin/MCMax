c-----------------------------------------------------------------------
c  This routine calculates the temperature dependent opacities.
c  
c
c-----------------------------------------------------------------------


      subroutine Topac(niter)
      use Parameters
	use InputOutput
      implicit none
      integer i,j,ii,iopac,nopac,niter
      logical InRange,InBetween
      real*8 checksum,lininterpol
	character*500 filename

	if(niter.ne.0.and.NphotUV.gt.0.and.UVdes) call CreateLRF(NphotUV,10000,.true.)
  
      ! loop over all grains
      do ii=1,ngrains 
         if (Grain(ii)%parttype .eq. 5) then
            write(*,'("Calculating opacities for particle ",i2)') ii
            write(9,'("Calculating opacities for particle ",i2)') ii
            nopac=Grain(ii)%nopac

            ! loop over all cells
            do i=1,D%nR-1 
               call tellertje(i,D%nR-1)
               do j=1,D%nTheta-1
                  
                  ! calc opacity using cell temperature
                  C(i,j)%wopac(ii,1:nopac)= 0d0

                  if (topac_interpol) then
                     if (C(i,j)%T.le.Grain(ii)%Topac(1))     C(i,j)%wopac(ii,1)     =1d0
                     if (C(i,j)%T.gt.Grain(ii)%Topac(nopac)) C(i,j)%wopac(ii,nopac) =1d0
                     do iopac=1,nopac-1
                        if (InBetween(C(i,j)%T,iopac,Grain(ii)%Topac,nopac)) then
                           C(i,j)%wopac(ii,iopac)=lininterpol(C(i,j)%T,
     &                        Grain(ii)%Topac(iopac),Grain(ii)%Topac(iopac+1),1d0,0d0)
                           C(i,j)%wopac(ii,iopac+1)=1d0-C(i,j)%wopac(ii,iopac)
                        endif
                     enddo
                  else
                     do iopac=1,nopac
                        if (InRange(C(i,j)%T,iopac,Grain(ii)%Topac,nopac)) then
                           C(i,j)%wopac(ii,iopac)=1d0
                        endif                        
                     enddo
                  endif

                  ! Check one gridcell (remove this)
                  if (0 .eq.1) then
!!                  if (j.eq.10) then
!!                  if (C(i,j)%T.gt.123.and.C(i,j)%T.le.124) then
                     write(*,*)
                     write(*,*) C(i,j)%T
                     write(*,*)
                     write(*,*) Grain(ii)%Topac
                     write(*,*)
                     write(*,*) C(i,j)%wopac(ii,1:nopac)
                     write(*,*)
!!                     stop 666
                  endif

                  ! Check cell
                  checksum=sum(C(i,j)%wopac(ii,1:nopac))
                  if (abs(checksum-1d0).ge.1d-6) then
                     write(*,'("Something wrong with wopac")')
                     do iopac=1,nopac
                        write(*,'(" wopac(",i02,",",i02,")= ",f10.8,f10.8)') 
     &                     ii,iopac,C(i,j)%wopac(ii,iopac),C(i,j)%T
     					print*,i,j,C(i,j)%dens,C(i,j)%w(1:ngrains)
                     enddo
                     write(*,'("Sum=",f10.8," =!= 1 -> stop 63635")')
     &                     checksum
                     stop 63635
                  endif

					if(niter.ne.0.and.NphotUV.gt.0.and.UVdes) call UVdestruction(ii,i,j)
               enddo            ! theta loop
            enddo               ! r loop
         endif                  ! graintype=5
      enddo                     ! graintype loop


	if(niter.ne.0.and.NphotUV.gt.0.and.UVdes) then
		if(outputfits) then
			write(filename,'(a,"/UVdestruction.fits.gz")') outdir(1:len_trim(outdir))
		else
			write(filename,'(a,"/UVdestruction.dat")') outdir(1:len_trim(outdir))
		endif
		call outputstruct(filename,(/'G0     ','NH     '/),2,0)
	endif



      return
      end

      ! This subroutine checks if T is close to Tarr(i)
      logical function InRange(T,i,Tarr,n)
      use Parameters
      implicit none
      integer i,n
      real*8 T,Tarr(1:n)

      InRange=.false.
      if (i.eq.1) then
         if (T .le. (Tarr(i)+Tarr(i+1))/2d0 ) InRange=.true.
      else if (i .eq. n) then
         if (T .gt. (Tarr(i-1)+Tarr(i))/2d0 ) InRange=.true.
      else
         if ( (T .le. (Tarr(i)+Tarr(i+1))/2d0 ) .and.
     &        (T .gt. (Tarr(i-1)+Tarr(i))/2d0 )) InRange=.true.
      endif

      return
      end

      ! This subroutine checks if T is between to Tarr(i) and Tarr(i+1)
      logical function InBetween(T,i,Tarr,n)
      use Parameters
      implicit none
      integer i,n
      real*8 T,Tarr(1:n)

      if (T.gt.Tarr(i).and.T.le.Tarr(i+1)) then 
         InBetween=.true.
      else
         InBetween=.false.
      endif

      return
      end



	subroutine UVdestruction(ii,i,j)
	use Parameters
	IMPLICIT NONE
	real*8 spec(nlam),lam1,lam2,nH,dens0,A,B,T
	integer i,j,l,ii
	real*8 mu,Rgas,mole,dsdt,Pv,G,Omega,GG
	parameter(mu=2.3*1.67262158d-24) !2.3 times the proton mass in gram
	parameter(Rgas=8.314e7) ! erg/mol/K
	parameter(mole=6.022e23)
	parameter(G=6.67300d-8) ! in cm^3/g/s

	lam1=0.0953
	lam2=0.206
	
	spec(1:nlam)=C(i,j)%LRF(1:nlam)
	do l=1,nlam
		if(lam(l).lt.lam1) spec(l)=0d0
		if(lam(l).gt.lam2) spec(l)=0d0
	enddo
	call integrate(spec,GG)

	GG=4d0*pi*GG/5.33d-14/3d10
c	if(GG.lt.1d-6) GG=1d-6

	nH = C(i,j)%gasdens*gas2dust / mu
	
	C(i,j)%G=GG
	C(i,j)%ne=1.5e-4*f_ne*nH

	T=C(i,j)%T
	if(T.lt.2.7d0) T=2.7d0

	A=Grain(ii)%TdesA
	B=Grain(ii)%TdesB
	dens0=10d0**(A/B-1d4/(T*B)-log10(T))+gammaUVdes*GG/sqrt(T)

c------------- determine the time for evaporation --------------
	Pv=(10d0**(A/B-1d4/(T*B)-log10(T))+gammaUVdes*GG/sqrt(T))*Rgas*T/(mu*mole)
	dsdt=(Pv/Grain(ii)%rho(1))*sqrt(mu/(2d0*pi*kb*T))

	Omega=2d0*pi*sqrt(D%R_av(i)**3/(G*D%Mstar))

	if((Grain(ii)%rv/dsdt).lt.Omega.and.C(i,j)%dens.lt.dens0) then
		C(i,j)%wopac(ii,1:Grain(ii)%nopac-1)=0d0
		C(i,j)%wopac(ii,Grain(ii)%nopac)=1d0
	endif

	return
	end

