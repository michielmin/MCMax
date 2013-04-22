c     --------------------------------------------------------------
c     This routine computes the value of the viscosity and turbulent
c     mixing strength in every grid cell. Both quantities are treated
c     seperately, unless specified by setting the Prandtl number to 1.
c
c     For a deadzone, there are three criteria for the active layer:
c     - T > deadtemp (ionization of alkali metals)
c     - column < deadcolumn (active column (Xrays/Cosmic rays))
c     - radial tau=1 > 1 (exclude first 5 gridcells)
c
c     --------------------------------------------------------------


      subroutine MakeDeadZone(printlocation)
      use Parameters
      implicit none

      real*8 Temp(1:D%nR-1,1:D%nTheta-1)

      integer i,j,imin,imax,jmin,jmax
      real*8 col,col2,prev_col,Tint,prev_Tint
      real*8 dcostheta(1:D%nTheta-1)
      real*8 talpha(1:D%nR-1,1:D%nTheta-1)
      character*500 filename
      real*8 loginterpol,lininterpol ! function
      logical printlocation

      ! Power law viscosity,exponential taper on mdot?
      do i=1,D%nR-1   
         do j=1,D%nTheta-1
            C(i,j)%alphavis=alphavis*(D%R_av(i)/AU)**alphavispow
         enddo
         if(raditer) then
            D%MdotR(i)=D%Mdot*exp(-(D%R_av(i)/(AU*D%Rexp))**2)
         endif
      enddo

      ! the Temperature
      call GetTemp(Temp(1:D%nR-1,1:D%nTheta-1))

      ! deadzone boundaries (radial)
      imin=D%nR
      imax=0
      jmin=D%nTheta
      jmax=0

      ! height of cells
      do j=1,D%nTheta-1
         dcostheta(j)=(D%Theta(j)-D%Theta(j+1)) 
      enddo

      ! Set turbulence in each cell for a deadzone
      if (deadzone) then
         do i=1,D%nR-1     

            col=0d0
            col2=0d0
            Tint=Temp(i,1)

            do j=1,D%nTheta-1
               
               ! vertical column
               prev_col=col
               col=col+ C(i,j)%gasdens*gas2dust*(D%R_av(i)*dcostheta(j))
c               col2=col2+C(i,j)%mass*gas2dust/(pi*(D%R(i+1)**2-D%R(i)**2)*AU**2)
c               print*,col,col2 ! ??

               ! Temperature at interface
               prev_Tint=Tint
c               print*,cos(D%theta_av(j+1)),D%Theta(j+1),cos(D%theta_av(j))
               if(j.ne.D%nTheta-1) then
                  Tint=loginterpol(D%Theta(j+1),cos(D%theta_av(j+1)),cos(D%theta_av(j)),
     &                             Temp(i,j+1),Temp(i,j))
               else
                  Tint=Temp(i,j) ! never used
               endif
c               print*,Temp(i,j+1),Tint,Temp(i,j)

               ! deadzone or upper boundary
               if (i.gt.5.and.Tint.lt.deadtemp.and.col.gt.deadcolumn) then
                  if(prev_col.gt.deadcolumn) then
                     if(prev_Tint.lt.deadtemp) then ! deadzone
                        C(i,j)%alphaturb=deadalpha
                     else  ! upper deadzone boundary (temperature)
                        C(i,j)%alphaturb=lininterpol(deadtemp,prev_Tint,Tint,
     &                                               deadalpha,alphaturb)
c                        print*,'temp: ',C(i,j)%alphaturb
c                        print*,prev_Tint,Temp(i,j),Tint
c                        print*
                     endif
                  else ! upper deadzone boundary (column)
                     C(i,j)%alphaturb=lininterpol(deadcolumn,prev_col,col,
     &                                            deadalpha,alphaturb)
c                     print*,'col: ',C(i,j)%alphaturb
c                     print*,prev_col,deadcolumn,col
                  endif

                  if(j.eq.D%nTheta-1) D%MPdead(i)=.true.

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
c                  print*,prev_Tint,Temp(i,j),Tint
c                  print*
                  if(j.eq.D%nTheta-1) D%MPdead(i)=.true.

               ! active layer
               else
                  C(i,j)%alphaturb=alphaturb
                  if(j.eq.D%nTheta-1) D%MPdead(i)=.false.
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
         if(printlocation.and.imin.le.imax) then
            write(*,'("Deadzone between:   ",f6.2," and ",f6.2," AU")') 
     &               D%R_av(imin)/AU,D%R_av(imax)/AU
            write(9,'("Deadzone between:   ",f6.2," and ",f6.2," AU")') 
     &               D%R_av(imin)/AU,D%R_av(imax)/AU
         endif
         if(printlocation.and.jmin.le.jmax) then
            write(*,'("Aspect ratio between: ",f4.2," and ",f4.2)') 
     &              pi/2d0-D%thet(jmax+1),pi/2d0-D%thet(jmin)
            write(9,'("Aspect ratio between: ",f4.2," and ",f4.2)') 
     &              pi/2d0-D%thet(jmax+1),pi/2d0-D%thet(jmin)
         endif
 
         !  Viscosity folows turbulence
         if (prandtl.gt.0d0) then
            do i=1,D%nR-1
               do j=1,D%nTheta-1
                  C(i,j)%alphavis=C(i,j)%alphaturb / prandtl
               enddo
            enddo
         endif

         !  Write alpha to a file, same format as denstemp
         write(filename,'(a,"deadzone.dat")') outdir(1:len_trim(outdir))
         write(9,'("Writing alpha to: ",a)') filename(1:len_trim(filename))
         call outputstruct(filename,(/'ALPHAT ','ALPHAV '/),2,0)

      endif

      end
