c-----------------------------------------------------------------------
c  This routine calculates the temperature dependent opacities.
c  
c
c-----------------------------------------------------------------------


      subroutine Topac()
      use Parameters
      implicit none
      integer i,j,ii,iopac,nopac
      logical InRange,InBetween
      real*8 checksum,lininterpol

      ! loop over all grains
      do ii=0,ngrains 
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
                        write(*,'(" wopac(",i02,",",i02,")= ",f10.8)') 
     &                     ii,iopac,C(i,j)%wopac(ii,iopac)
                     enddo
                     write(*,'("Sum=",f10.8," =!= 1 -> stop 63635")')
     &                     checksum
                     stop 63635
                  endif

               enddo            ! theta loop
            enddo               ! r loop
         endif                  ! graintype=5
      enddo                     ! graintype loop

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
