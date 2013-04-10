c     --------------------------------------------------------------
c     This routine computes the value of the viscosity/turbulent
c     mixing strength in every grid cell.
c
c     For a deadzone, there are two criteria for the active layer:
c     - T > deadtemp (ionization of alkali metals)
c     - column < deadcolumn (active column (Xrays/Cosmic rays))
c     
c     --------------------------------------------------------------


      subroutine MakeDeadZone()
      use Parameters
      implicit none

      integer i,j
      real*8 col
      real*8 dcostheta(1:D%nTheta-1)
      character*500 filename

      ! height of cells
      do j=1,D%nTheta-1
         dcostheta(j)=(D%thet(j+1)-D%thet(j)) 
c         write(*,*) j,D%theta_av(j),D%theta_av(j+1),D%thet(j),dcostheta(j)
c         write(*,*) j,D%thet(j),dcostheta(j)
      enddo

      ! Global turbulent mixing strength and active layers
      C(1:D%nR-1,1:D%nTheta-1)%alpha=alphaturb

      ! Set turbulence in each cell for a deadzone
      if (deadzone) then
         do i=1,D%nR-1     

            col=0d0            
            do j=1,D%nTheta-1
               
               ! vertical column
               col=col+ C(i,j)%gasdens*(D%R_av(i)*dcostheta(j))

               if (C(i,j)%T.lt.deadtemp.and.col.gt.deadcolumn) then
                  C(i,j)%alpha=deadalpha
                  if (j.eq.D%nTheta-1) then
                     write(*,'("deadzone at r= ",f8.2," z/r= ",f4.2)') 
     &                    D%R_av(i)/AU, D%theta_av(j)
                  endif
               endif
       
            enddo
            write(*,'("  at r= ",f8.2," column= ",f8.2," gr/cm^2")') 
     &           D%R_av(i)/AU, col
         enddo
        
         !  Write alpha to a file, same format as denstemp
         write(filename,'(a,"deadzone.dat")') outdir(1:len_trim(outdir))
         write(9,'("Writing alpha to: ",a)') filename(1:len_trim(filename))
         call outputstruct(filename,(/'ALPHA  '/),1,0)

      endif

      end
