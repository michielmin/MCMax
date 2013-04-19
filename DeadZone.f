c     --------------------------------------------------------------
c     This routine computes the value of the viscosity and turbulent
c     mixing strength in every grid cell. Both quantities are treated
c     seperately, unless specified by setting the Prandtl number to 1.
c
c     For a deadzone, there are two criteria for the active layer:
c     - T > deadtemp (ionization of alkali metals)
c     - column < deadcolumn (active column (Xrays/Cosmic rays))
c     
c     --------------------------------------------------------------


      subroutine MakeDeadZone()
      use Parameters
      implicit none

      integer i,j,imin,imax,jmin,jmax
      real*8 col,prev_col,Tint,prev_Tint
      real*8 dcostheta(1:D%nTheta-1)
      real*8 talpha(1:D%nR-1,1:D%nTheta-1)
      character*500 filename
      real*8 loginterpol,lininterpol ! function
c      logical bycol,bytemp

      ! deadzone boundaries (radial)
      imin=D%nR
      imax=0
      jmin=D%nTheta
      jmax=0

      ! height of cells
      do j=1,D%nTheta-1
         dcostheta(j)=(D%Theta(j)-D%Theta(j+1)) 
      enddo

      ! Power law viscosity
      do i=1,D%nR-1     
         C(i,1:D%nTheta-1)%alphavis=alphavis*(D%R_av(i)/AU)**alphavispow
      enddo

      ! Set turbulence in each cell for a deadzone
      if (deadzone) then
         do i=1,D%nR-1     

            col=0d0
            Tint=C(i,1)%T

            do j=1,D%nTheta-1
               
               ! vertical column
               prev_col=col
               col=col+ C(i,j)%gasdens*gas2dust*(D%R_av(i)*dcostheta(j))

               ! Temperature at interface
               prev_Tint=Tint
c               print*,cos(D%theta_av(j+1)),D%Theta(j+1),cos(D%theta_av(j))
               if(j.ne.D%nTheta-1) then
                  Tint=loginterpol(D%Theta(j+1),cos(D%theta_av(j+1)),cos(D%theta_av(j)),
     &                             C(i,j+1)%T,C(i,j)%T)
               else
                  Tint=C(i,j)%T ! never used
               endif
c               print*,C(i,j+1)%T,Tint,C(i,j)%T

               ! deadzone or upper boundary
               if (Tint.lt.deadtemp.and.col.gt.deadcolumn) then
                  if(prev_col.gt.deadcolumn) then
                     if(prev_Tint.lt.deadtemp) then ! deadzone
                        C(i,j)%alphaturb=deadalpha
                     else  ! upper deadzone boundary (temperature)
                        C(i,j)%alphaturb=lininterpol(deadtemp,prev_Tint,Tint,
     &                                               deadalpha,alphaturb)
c                        print*,'temp: ',C(i,j)%alphaturb
c                        print*,prev_Tint,C(i,j)%T,Tint
c                        print*
                     endif
                  else ! upper deadzone boundary (column)
                     C(i,j)%alphaturb=lininterpol(deadcolumn,prev_col,col,
     &                                            deadalpha,alphaturb)
c                     print*,'col: ',C(i,j)%alphaturb
c                     print*,prev_col,deadcolumn,col
                  endif

                  imin=min(i,imin)
                  imax=max(i,imax)
                  jmin=min(j,jmin)
                  jmax=max(j,jmax)

               ! lower boundary of deadzone
               elseif (prev_Tint.lt.deadtemp.and.Tint.ge.deadtemp
     &                 .and.col.gt.deadcolumn) then
                  C(i,j)%alphaturb=lininterpol(deadtemp,prev_Tint,Tint,
     &                                         alphaturb,deadalpha)
c                  print*,'temp: ',C(i,j)%alphaturb
c                  print*,prev_Tint,C(i,j)%T,Tint
c                  print*

               else
                  C(i,j)%alphaturb=alphaturb
               endif
       
            enddo ! j=1,D%nTheta-1
c            write(*,'("  at r= ",f8.2," column= ",f8.2," gr/cm^2")') 
c     &           D%R_av(i)/AU, col
         enddo ! i

         ! Smooth array?
         if(.false.) then
            talpha(1:D%nR-1,1:D%nTheta-1)=C(1:D%nR-1,1:D%nTheta-1)%alphaturb
            call LogSmooth2D(talpha,D%nR-1,D%nTheta-1)
            C(1:D%nR-1,1:D%nTheta-1)%alphaturb=talpha(1:D%nR-1,1:D%nTheta-1)
         endif
       
         !  Print dead zone location (if present)
         if(imin.le.imax) then
            write(*,'("Deadzone between:   ",f6.2," and ",f6.2," AU")') 
     &               D%R_av(imin)/AU,D%R_av(imax)/AU
            write(9,'("Deadzone between:   ",f6.2," and ",f6.2," AU")') 
     &               D%R_av(imin)/AU,D%R_av(imax)/AU
         endif
         if(jmin.le.jmax) then
            write(*,'("Aspect ratio between: ",f4.2," and ",f4.2)') 
     &              pi/2d0-D%thet(jmax+1),pi/2d0-D%thet(jmin)
            write(9,'("Aspect ratio between: ",f4.2," and ",f4.2)') 
     &              pi/2d0-D%thet(jmax+1),pi/2d0-D%thet(jmin)
         endif
 
         !  Viscosity folows turbulence
         if (prandtl.gt.0d0) then
            C(1:D%nR-1,1:D%nTheta-1)%alphavis=C(1:D%nR-1,1:D%nTheta-1)%alphaturb / prandtl
         endif

         !  Write alpha to a file, same format as denstemp
         write(filename,'(a,"deadzone.dat")') outdir(1:len_trim(outdir))
         write(9,'("Writing alpha to: ",a)') filename(1:len_trim(filename))
         call outputstruct(filename,(/'ALPHAT ','ALPHAV '/),2,0)

      endif

      end
