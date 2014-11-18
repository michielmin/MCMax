	module InputOutput
	IMPLICIT NONE

	contains
	
	subroutine outputstruct(filename,vars,nvars,ipart,xx,nxx)
	use Parameters
	IMPLICIT NONE
	integer nvars,ivars,i,j,ii,ipart,l
	integer,optional :: nxx
	real*8,optional :: xx(D%nR-1,D%nTheta-1,*)
	character*7 vars(nvars)
	character*500 filename

	if(outputfits) then
		if(present(xx)) then
			call outputstruct_fits(filename,vars,nvars,ipart,xx,nxx)
		else
			call outputstruct_fits(filename,vars,nvars,ipart)
		endif
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
			case ('ALPHAT')
				write(20,'("# Turbulent mixing strength array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%alphaturb
					enddo
				enddo
			case ('ALPHAV')
				write(20,'("# Viscosity array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%alphavis
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

			case ('TQHP')
				write(20,'("# Average T for QHP (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%Tqhp(Grain(ipart)%qhpnr)
					enddo
				enddo
			case ('OPACITY')

			case ('G0')
				write(20,'("# G0 (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%G
					enddo
				enddo
			case ('NE')
				write(20,'("# Ne (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%ne
					enddo
				enddo
			case ('NH')
				write(20,'("# NH (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						write(20,*) C(i,j)%gasdens*gas2dust / (1.3*1.67262158d-24)
					enddo
				enddo
			case ('ARRAY')
				do l=1,nxx
					write(20,'("# array (for ir=0,nr-1 do for it=0,nt-1 do ...)")')
					do i=1,D%nR-1
						do j=1,D%nTheta-1
							write(20,*) xx(i,j,l)
						enddo
					enddo
				enddo
			case default
				write(9,'("Error in output file specification")')
				write(*,'("Error in output file specification")')
				close(unit=20)
				stop
		end select
	enddo

	close(unit=20)
	
	return
	end subroutine outputstruct


	subroutine outputstruct_fits(filename,vars,nvars,ipart,xx,nxx)
	use Parameters
	IMPLICIT NONE
	integer nvars,ivars,i,j,ii,ipart,iopac,l
	integer,optional :: nxx
	character*7 vars(nvars),hdu
	character*500 filename
	real*8,allocatable :: array(:,:,:,:)
	real*8,optional :: xx(D%nR-1,D%nTheta-1,*)
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
	nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Write optional keywords to the header

	call ftpkyd(unit,'Rin',D%Rin,8,'[AU]',status)
	call ftpkyd(unit,'Rout',D%Rout,8,'[AU]',status)

	call ftpkyj(unit,'nR',D%nR-1,' ',status)
	call ftpkyj(unit,'nTheta',D%nTheta-1,' ',status)
	call ftpkyj(unit,'ngrains',ngrains,' ',status)
	call ftpkyj(unit,'ngrains2',ngrains2,' ',status)
	call ftpkyj(unit,'nlam',nlam,' ',status)

	call ftpkyj(unit,'nHDU',nvars,' ',status)	
	do i=1,nvars
		write(hdu,'("HDU",i2)') i
		call ftpkys(unit,hdu,trim(vars(i)),'',status)
	enddo

	!  Write the array to the FITS file.

	!------------------------------------------------------------------------------
	! HDU 0: spatial grid
	!------------------------------------------------------------------------------

	allocate(array(D%nR-1,D%nTheta-1,2,1))

	do i=1,D%nR-1
		do j=1,D%nTheta-1
			array(i,j,1,1)=D%R_av(i)/AU
			array(i,j,2,1)=D%theta_av(j)
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
			case ('ALPHAT')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%alphaturb
					enddo
				enddo
			case ('ALPHAV')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%alphavis
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
				nelements=naxes(1)*naxes(2)*naxes(3)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						do l=1,nlam
							array(i,j,l,1)=C(i,j)%QHP(Grain(ipart)%qhpnr,l)
						enddo
					enddo
				enddo
			case ('QHPT')
				naxis=4
				naxes(3)=NTQHP
				naxes(4)=2
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						do l=1,NTQHP
							array(i,j,l,1)=tgridqhp(l)
							array(i,j,l,2)=C(i,j)%tdistr(Grain(ipart)%qhpnr,l)
						enddo
					enddo
				enddo
			case ('TQHP')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%Tqhp(Grain(ipart)%qhpnr)
					enddo
				enddo
			case ('G0')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%G
					enddo
				enddo
			case ('NE')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%ne
					enddo
				enddo
			case ('NH')
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						array(i,j,1,1)=C(i,j)%gasdens*gas2dust / (1.3*1.67262158d-24)
					enddo
				enddo
			case ('OPACITY')
				naxis=4
				naxes(3)=nlam
				naxes(4)=3
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						do l=1,nlam
							array(i,j,l,1:3)=0d0
							do ii=1,ngrains
								do iopac=1,Grain(ii)%nopac
									array(i,j,l,1)=array(i,j,l,1)+Grain(ii)%Kext(iopac,l)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
									array(i,j,l,2)=array(i,j,l,2)+Grain(ii)%Kabs(iopac,l)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
									array(i,j,l,3)=array(i,j,l,3)+Grain(ii)%Ksca(iopac,l)*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac)
								enddo
							enddo
						enddo
					enddo
				enddo
			case ('ARRAY')
				if(nxx.gt.1) then
					naxis=3
					naxes(3)=nxx
					nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
				endif
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						do l=1,nxx
							array(i,j,l,1)=xx(i,j,l)
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
	   print*,'error in export to fits file',status
	end if


	return
	end subroutine outputstruct_fits




	subroutine readstruct(filename,vars,nvars,ipart,doalloc)
	use Parameters
	IMPLICIT NONE
	integer nvars,ivars,i,j,ii,ipart,l,nr,nt,ngrains2_tmp,ngrains_tmp,lf
	character*7 vars(nvars)
	character*500 filename
	character*1000 line
	logical doalloc,truefalse

	lf=len_trim(filename)
	if(outputfits.and.(filename(lf-4:lf).eq.'.fits'.or.filename(lf-7:lf).eq.'.fits.gz')) then
		call readstruct_fits(filename,vars,nvars,ipart,doalloc)
		return
	endif

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
	ngrains_tmp=-1
	ngrains2_tmp=-1
	read(20,'(a1000)') line(1:1000)
	read(line,*,end=1) nr,nt,ngrains_tmp,ngrains2_tmp

1	continue
	if(ngrains_tmp.gt.0) ngrains=ngrains_tmp
	if(ngrains2_tmp.gt.0) then
		ngrains2=ngrains2_tmp
	else
		ngrains2_tmp=ngrains2
		ngrains2=1
	endif
	if(doalloc) then
		D%nR=nr+1
		D%nTheta=nt+1
		if(.not.arraysallocated.and..not.allocated(C)) then
			allocate(C(0:D%nR+1,0:D%nTheta+1))
			allocate(D%theta_av(0:D%nTheta+1))
			allocate(D%Theta(0:D%nTheta+1))
			allocate(D%thet(0:D%nTheta+1))
			allocate(D%SinTheta(0:D%nTheta+1))
			allocate(D%R_av(0:D%nR+1))
			allocate(D%R(0:D%nR+1))
			allocate(shscale(0:D%nR+1))
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
			case ('QHPRF')

			case ('TQHP')

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

	ngrains2=ngrains2_tmp

	return
	end subroutine readstruct







	subroutine readstruct_fits(filename,vars,nvars,ipart,doalloc)
	use Parameters
	IMPLICIT NONE
	integer nvars,ivars,i,j,ii,ipart,l,nr,nt,iopac,naxis,nhdu
	character*7 vars(nvars)
	character*500 filename
	logical doalloc,truefalse
	real*8,allocatable :: array(:,:,:,:)
	integer*4 :: status,stat2,stat3,readwrite,unit,blocksize,nfound,group
	integer*4 :: firstpix,nbuffer,npixels,hdunum,hdutype,ix,iz,ilam
	integer*4 :: istat,stat4,tmp_int,stat5,stat6
	real*8  :: nullval
	logical*4 :: anynull
	integer*4, dimension(4) :: naxes
	character*80 comment,errmessage
	character*30 errtext

	! Get an unused Logical Unit Number to use to open the FITS file.
	status=0

	call ftgiou (unit,status)
	! Open file
	readwrite=0
	call ftopen(unit,filename,readwrite,blocksize,status)
	if (status /= 0) then
		write(*,'("Density file not found")')
		write(9,'("Density file not found")')
		print*,trim(filename)
		write(*,'("--------------------------------------------------------")')
		write(9,'("--------------------------------------------------------")')
		stop
	endif
	group=1
	firstpix=1
	nullval=-999


	!------------------------------------------------------------------------
	! HDU0 : grid
	!------------------------------------------------------------------------
	! Check dimensions
	call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

	npixels=naxes(1)*naxes(2)*naxes(3)

	! Read model info

	call ftgkyd(unit,'Rin',D%Rin,comment,status)
	call ftgkyd(unit,'Rout',D%Rout,comment,status)

	call ftgkyj(unit,'nR',D%nR,comment,status)
	call ftgkyj(unit,'nTheta',D%nTheta,comment,status)
	call ftgkyj(unit,'ngrains',ngrains,comment,status)
	call ftgkyj(unit,'ngrains2',ngrains2,comment,status)
	call ftgkyj(unit,'nlam',nlam,comment,status)
	call ftgkyj(unit,'nHDU',nhdu,comment,status)
	if(status.ne.0) then
		nhdu=nvars
		status=0
	endif

	D%nR=D%nR+1
	D%nTheta=D%nTheta+1
	if(doalloc) then
		if(.not.arraysallocated.and..not.allocated(C)) then
			allocate(C(0:D%nR+1,0:D%nTheta+1))
			allocate(D%theta_av(0:D%nTheta+1))
			allocate(D%Theta(0:D%nTheta+1))
			allocate(D%thet(0:D%nTheta+1))
			allocate(D%SinTheta(0:D%nTheta+1))
			allocate(D%R_av(0:D%nR+1))
			allocate(D%R(0:D%nR+1))
			allocate(shscale(0:D%nR+1))
		endif
	endif
 
	! read_image
	allocate(array(D%nR-1,D%nTheta-1,2,1))

	call ftgpvd(unit,group,firstpix,npixels,nullval,array,anynull,status)

	do i=1,D%nR-1
		D%R_av(i)=array(i,1,1,1)*AU
	enddo
	do j=1,D%nTheta-1
		D%theta_av(j)=array(1,j,2,1)
	enddo

	deallocate(array)

	do ivars=1,minval((/nhdu,nvars/))
		!  move to next hdu
		call ftmrhd(unit,1,hdutype,status)
		if(status.ne.0) then
			nhdu=ivars
			status=0
			goto 1
		endif
		select case (vars(ivars))
			case ('COMP')
				naxis=4
			case ('LAM')
				naxis=1
			case ('QHPRF')
				naxis=3
			case default
				naxis=2
		end select

		! Check dimensions
		call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)

		do i=naxis+1,4
			naxes(i)=1
		enddo
		npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

		! read_image
		allocate(array(naxes(1),naxes(2),naxes(3),naxes(4)))

		call ftgpvd(unit,group,firstpix,npixels,nullval,array,anynull,status)
   
		select case (vars(ivars))
			case ('DENS')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%dens=array(i,j,1,1)
					enddo
				enddo
			case ('TEMP')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%T=array(i,j,1,1)
					enddo
				enddo
			case ('TEMPMC')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%TMC=array(i,j,1,1)
					enddo
				enddo
			case ('COMP')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						if(.not.allocated(C(i,j)%w)) allocate(C(i,j)%w(ngrains))
						do ii=1,ngrains
							C(i,j)%w(ii)=array(i,j,ii,1)
							do iopac=1,ngrains2
								C(i,j)%wopac(ii,iopac)=array(i,j,ii,iopac+1)
							enddo
						enddo
					enddo
				enddo
			case ('GASDENS')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%gasdens=array(i,j,1,1)/gas2dust
					enddo
				enddo
			case ('DENS0')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%dens0=array(i,j,1,1)
					enddo
				enddo
			case ('GASTEMP')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%Tgas=array(i,j,1,1)
					enddo
				enddo
			case ('DENSP')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%w(ipart)=array(i,j,1,1)
					enddo
				enddo
			case ('TEMPP')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%TP(ipart)=array(i,j,1,1)
					enddo
				enddo
			case ('NPHOT')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%Ni=array(i,j,1,1)
					enddo
				enddo
			case ('VOLUME')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%V=array(i,j,1,1)
					enddo
				enddo
			case ('dTEMP')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%dT=array(i,j,1,1)*C(i,j)%T
					enddo
				enddo
			case ('dEJv')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%dEJv=array(i,j,1,1)*C(i,j)%EJv
					enddo
				enddo
			case ('LAM')
				do i=1,nlam
					lam(i)=array(i,1,1,1)
				enddo
			case ('QHPEJv')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%EJvQHP(Grain(ipart)%qhpnr)=array(i,j,1,1)
					enddo
				enddo
			case ('QHPRF')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						do l=1,nlam
							C(i,j)%QHP(Grain(ipart)%qhpnr,l)=array(i,j,l,1)
						enddo
					enddo
				enddo
			case ('TQHP')
				do i=1,D%nR-1
					do j=1,D%nTheta-1
						C(i,j)%Tqhp(Grain(ipart)%qhpnr)=array(i,j,1,1)
					enddo
				enddo
			case ('SKIP')
c	just skip this hdu
			case default
				write(9,'("Error in output file specification")')
				write(*,'("Error in output file specification")')
				print*,vars(ivars)
				stop
		end select
		deallocate(array)
	enddo

1	continue
	
	if(nvars.gt.nhdu) then
		do ivars=nhdu+1,nvars
			select case(vars(ivars))
				case ('DENS0')
					write(*,'("** Assuming gasdens=dens0 ! **")')
					write(9,'("** Assuming gasdens=dens0 ! **")')
					do i=1,D%nR-1
						do j=1,D%nTheta-1
							C(i,j)%dens0=C(i,j)%gasdens
						enddo
					enddo
				case default
					stop
			end select
		enddo
	endif


	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	!  Get the text string which describes the error
	if (status > 0) then
	   call ftgerr(status,errtext)
	   print *,'FITSIO Error Status =',status,': ',errtext

	   !  Read and print out all the error messages on the FITSIO stack
	   call ftgmsg(errmessage)
	   do while (errmessage .ne. ' ')
		  print *,errmessage
		  call ftgmsg(errmessage)
	   end do
	endif

	return
	end subroutine readstruct_fits


	end module InputOutput



