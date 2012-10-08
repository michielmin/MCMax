	subroutine outputstruct(filename,vars,nvars,ipart)
	use Parameters
	IMPLICIT NONE
	integer nvars,ivars,i,j,ii,ipart,l
	character*7 vars(nvars)
	character*500 filename

	if(outputfits) then
		call outputstruct_fits(filename,vars,nvars,ipart)
		return
	endif

	open(unit=20,file=filename,RECL=6000)
	write(20,'("# Format number")')
	write(20,'(i6)') 5 ! including composition and gas density (including gas2dust ratio) and ngrains!
                           ! including ngrains2 for temperature dependent opacities
	do i=1,nvars
		if(vars(i)(1:3).eq.'QHP') then	
			write(20,'("# NR, NT, NLAM")')
			write(20,'(4(i6))') D%nR-1,D%nTheta-1,nlam
			goto 1
		endif
	enddo
	write(20,'("# NR, NT, NGRAINS, NGRAINS2")')
	write(20,'(4(i6))') D%nR-1,D%nTheta-1,ngrains,ngrains2

1	continue
	write(20,'("# Spherical radius grid [cm] (middle of cell)")')
	do i=1,D%nR-1
		write(20,*) D%R_av(i)
	enddo
	write(20,'("# Theta grid [rad, from pole] (middle of cell)")')
	do i=1,D%nTheta-1
		write(20,*) D%theta_av(i)
	enddo

	do ivars=1,nvars
		select case (vars(ivars))
			case ('DENS')
				write(20,'("# Density array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%dens
					enddo
				enddo
			case ('TEMP')
				write(20,'("# Temperature array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%T
					enddo
				enddo
			case ('TEMPMC')
				write(20,'("# Temperature array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%TMC
					enddo
				enddo
			case ('COMP')
				write(20,'("# Composition array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) (C(i,j)%w(ii),ii=1,ngrains)
						if(ngrains2.gt.1) then
							do ii=1,ngrains
								write(20,*) C(i,j)%wopac(ii,1:ngrains2)
							enddo
						endif
					enddo
				enddo
			case ('GASDENS')
				write(20,'("# Gas density array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%gasdens*gas2dust
					enddo
				enddo
			case ('DENS0')
				write(20,'("# Density0 array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%dens0
					enddo
				enddo
			case ('GASTEMP')
				write(20,'("# Temperature array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%Tgas
					enddo
				enddo
			case ('DENSP')
				write(20,'("# Density array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%dens*C(i,j)%w(ipart)
					enddo
				enddo
			case ('TEMPP')
				write(20,'("# Temperature array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%TP(ipart)
					enddo
				enddo
			case ('NPHOT')
				write(20,'("# Number of photons (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%Ni
					enddo
				enddo
			case ('VOLUME')
				write(20,'("# Volume of the cell (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%V
					enddo
				enddo
			case ('dTEMP')
				write(20,'("# Relative dT (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%dT/C(i,j)%T
					enddo
				enddo
			case ('dEJv')
				write(20,'("# Relative dEJv (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%dEJv/C(i,j)%EJv
					enddo
				enddo
			case ('LAM')
				write(20,'("# Wavelength grid")')
				do i=1,nlam
					write(20,*) lam(i)
				enddo
			case ('QHPEJv')
				write(20,'("# QHP emission (for ir=0,nr-1 do for it=0,nt-1 do for ilam=0,nlam-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%EJvQHP(Grain(ipart)%qhpnr)
						do l=1,nlam
							write(20,*) C(i,j)%QHP(Grain(ipart)%qhpnr,l)
						enddo
					enddo
				enddo
			case ('QHPRF')

			case default
				write(9,'("Error in output file specification")')
				write(*,'("Error in output file specification")')
				close(unit=20)
				stop
		end select
	enddo

	close(unit=20)
	
	return
	end


	subroutine outputstruct_fits(filename,vars,nvars,ipart)
	use Parameters
	IMPLICIT NONE
	integer nvars,ivars,i,j,ii,ipart,iopac,l
	character*7 vars(nvars)
	character*500 filename
	real*8,allocatable :: array(:,:,:,:)
	integer status,unit,blocksize,bitpix,naxis,naxes(4)
	integer group,fpixel,nelements
	logical simple,extend,truefalse

	if(filename(len_trim(filename)-4:len_trim(filename)).eq.'.fits') then
		filename=trim(filename)//'.gz'
	endif

	inquire(file=filename,exist=truefalse)
	if(truefalse) then
		write(*,'("FITS file already exists, overwriting")')
		write(9,'("FITS file already exists, overwriting")')
		open(unit=90,file=filename)
		close(unit=90,status='delete')
	endif

	status=0
C	 Get an unused Logical Unit Number to use to create the FITS file
	call ftgiou(unit,status)
C	 create the new empty FITS file
	blocksize=1
	call ftinit(unit,filename,blocksize,status)

	simple=.true.
	extend=.true.
	group=1
	fpixel=1

	bitpix=-64
	naxis=3
	naxes(1)=D%nR-1
	naxes(2)=D%nTheta-1
	naxes(3)=2
	naxes(4)=1
	nelements=naxes(1)*naxes(2)

	! Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Write optional keywords to the header

	call ftpkyd(unit,'Rin',real(D%Rin),8,'[AU]',status)
	call ftpkyd(unit,'Rout',real(D%Rout),8,'[AU]',status)

	call ftpkyj(unit,'nR',D%nR-1,' ',status)
	call ftpkyj(unit,'nTheta',D%nTheta-1,' ',status)
	call ftpkyj(unit,'ngrains',ngrains,' ',status)
	call ftpkyj(unit,'ngrains2',ngrains2,' ',status)
	call ftpkyj(unit,'nlam',nlam,' ',status)

	!  Write the array to the FITS file.

	!------------------------------------------------------------------------------
	! HDU 0: spatial grid
	!------------------------------------------------------------------------------

	allocate(array(D%nR-1,D%nTheta-1,2,1))

	do i=1,D%nR-1
		do j=1,D%nTheta-1
			array(i,j,1,1)=D%R_av(i)/AU
			array(i,j,2,1)=D%theta_av(i)
		enddo
	enddo

	call ftpprd(unit,group,fpixel,nelements,array,status)
	
	deallocate(array)


	do ivars=1,nvars
		! create new hdu
		call ftcrhd(unit, status)
		bitpix=-64
		naxes(1)=D%nR-1
		naxes(2)=D%nTheta-1
		naxes(3)=1
		naxes(4)=1
		naxis=2
		nelements=naxes(1)*naxes(2)
		select case (vars(ivars))
			case ('DENS')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%dens
					enddo
				enddo
			case ('TEMP')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%T
					enddo
				enddo
			case ('TEMPMC')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%TMC
					enddo
				enddo
			case ('COMP')
				naxis=4
				naxes(3)=ngrains
				naxes(4)=ngrains2+1
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						do ii=1,ngrains
							array(i,j,ii,1)=C(i,j)%w(ii)
							do iopac=1,ngrains2
								array(i,j,ii,iopac+1)=C(i,j)%wopac(ii,iopac)
							enddo
						enddo
					enddo
				enddo
			case ('GASDENS')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%gasdens*gas2dust
					enddo
				enddo
			case ('DENS0')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%dens0
					enddo
				enddo
			case ('GASTEMP')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%Tgas
					enddo
				enddo
			case ('DENSP')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%dens*C(i,j)%w(ipart)
					enddo
				enddo
			case ('TEMPP')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%TP(ipart)
					enddo
				enddo
			case ('NPHOT')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%Ni
					enddo
				enddo
			case ('VOLUME')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%V
					enddo
				enddo
			case ('dTEMP')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%dT/C(i,j)%T
					enddo
				enddo
			case ('dEJv')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%dEJv/C(i,j)%EJv
					enddo
				enddo
			case ('LAM')
				naxis=1
				naxes(1)=nlam
				naxes(2)=1
				naxes(3)=1
				naxes(4)=1
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,nlam
					array(i,1,1,1)=lam(i)
				enddo
			case ('QHPEJv')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%EJvQHP(Grain(ipart)%qhpnr)
					enddo
				enddo
			case ('QHPRF')
				naxis=3
				naxes(3)=nlam
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						do l=1,nlam
							array(i,j,l,1)=C(i,j)%QHP(Grain(ipart)%qhpnr,l)
						enddo
					enddo
				enddo
			case default
				write(9,'("Error in output file specification")')
				write(*,'("Error in output file specification")')
				print*,vars(ivars)
				close(unit=20)
				stop
		end select


		!  Write the required header keywords.
		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
		!  Write the array to the FITS file.
		call ftpprd(unit,group,fpixel,nelements,array,status)

		deallocate(array)
	enddo
	
	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	if (status.gt.0) then
	   print*,'error in export to fits file'
	end if


	return
	end




	subroutine readstruct(filename,vars,nvars,ipart,doalloc)
	use Parameters
	IMPLICIT NONE
	integer nvars,ivars,i,j,ii,ipart,l,nr,nt
	character*7 vars(nvars)
	character*500 filename
	logical doalloc,truefalse

c	if(outputfits) then
c		call readstruct_fits(filename,vars,nvars,ipart)
c		return
c	endif

	inquire(file=filename,exist=truefalse)
	if(.not.truefalse) then
		write(*,'("Density file not found")')
		write(9,'("Density file not found")')
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		stop
	endif

	open(unit=20,file=filename,RECL=6000)
	read(20,*) ! comments
	read(20,*) ! format number
	read(20,*) ! comments
	do i=1,nvars
		if(vars(i)(1:3).eq.'QHP') then	
			read(20,*) nr,nt,nlam
			goto 1
		endif
	enddo
	read(20,*,end=1) nr,nt,ngrains,ngrains2

1	continue
	if(doalloc) then
		D%nR=nr+1
		D%nTheta=nt+1
		if(.not.arraysallocated) then
			allocate(C(0:D%nR+1,0:D%nTheta+1))
			allocate(D%theta_av(0:D%nTheta+1))
			allocate(D%Theta(0:D%nTheta+1))
			allocate(D%thet(0:D%nTheta+1))
			allocate(D%SinTheta(0:D%nTheta+1))
			allocate(D%R_av(0:D%nR+1))
			allocate(D%R(0:D%nR+1))
			allocate(shscale(0:D%nR+1))
		endif
		if(nr.ne.D%nR-1.or.nt.ne.D%nTheta-1) then
			write(*,'("File ",a," incompatible with spatial grid")') filename(1:len_trim(filename))
			write(9,'("File ",a," incompatible with spatial grid")') filename(1:len_trim(filename))
			stop
		endif
	endif
	read(20,*) ! comments
	do i=1,D%nR-1
		read(20,*) D%R_av(i)
	enddo
	read(20,*) ! comments
	do i=1,D%nTheta-1
		read(20,*) D%theta_av(i)
	enddo

	do ivars=1,nvars
		select case (vars(ivars))
			case ('DENS')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%dens
					enddo
				enddo
			case ('TEMP')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%T
					enddo
				enddo
			case ('TEMPMC')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%TMC
					enddo
				enddo
			case ('COMP')
				read(20,*,end=3) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*,end=3) (C(i,j)%w(ii),ii=1,ngrains)
						if(ngrains2.gt.1) then
							do ii=1,ngrains
								read(20,*,end=3) C(i,j)%wopac(ii,1:ngrains2)
							enddo
						endif
					enddo
				enddo
			case ('GASDENS')
				read(20,*,end=3) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*,end=3) C(i,j)%gasdens
						C(i,j)%gasdens=C(i,j)%gasdens/gas2dust
					enddo
				enddo
			case ('DENS0')
				read(20,*,end=2) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*,end=2) C(i,j)%dens0
					enddo
				enddo
			case ('GASTEMP')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%Tgas
					enddo
				enddo
			case ('DENSP')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%w(ipart)
					enddo
				enddo
			case ('TEMPP')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%TP(ipart)
					enddo
				enddo
			case ('NPHOT')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%Ni
					enddo
				enddo
			case ('VOLUME')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%V
					enddo
				enddo
			case ('dTEMP')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%dT
						C(i,j)%dT=C(i,j)%dT*C(i,j)%T
					enddo
				enddo
			case ('dEJv')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%dEJv
						C(i,j)%dEJv=C(i,j)%dEJv*C(i,j)%EJv
					enddo
				enddo
			case ('LAM')
				read(20,*) ! comments
				do i=1,nlam
					read(20,*) lam(i)
				enddo
			case ('QHPEJv')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*) C(i,j)%EJvQHP(Grain(ipart)%qhpnr)
						do l=1,nlam
							read(20,*) C(i,j)%QHP(Grain(ipart)%qhpnr,l)
						enddo
					enddo
				enddo
			case ('SKIP')
				read(20,*) ! comments
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						read(20,*)
					enddo
				enddo
			case default
				write(9,'("Error in input file specification")')
				write(*,'("Error in input file specification")')
				close(unit=20)
				stop
		end select
	enddo

	close(unit=20)
	
	return
2	continue

	write(*,'("** Assuming gasdens=dens0 ! **")')
	write(9,'("** Assuming gasdens=dens0 ! **")')
	do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%dens0=C(i,j)%gasdens
		enddo
	enddo
	close(unit=20)

3	continue

	return
	end





