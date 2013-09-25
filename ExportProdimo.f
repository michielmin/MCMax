	subroutine DoExportProdimo()
	use Parameters
	IMPLICIT NONE

	  integer status,unit,blocksize,bitpix,naxis,naxes(4),j0
	  integer j,group,fpixel,nelements,ri,zj,l,lambda,iopac,i
	  logical simple,extend,truefalse
	character*500 filename,version
	real*8,allocatable :: grid(:,:,:)
	real,allocatable :: temperature(:,:,:)
	real,allocatable :: spectre(:),J_io(:,:,:),dens(:,:),region_index(:)
	real,allocatable :: opacite(:,:,:,:),N_grains(:,:,:),HRspec(:,:)
	real N,N1,Mdust
	real*8 fUV,tot,pUV
	real*8,allocatable :: spec(:)
	character*2 s

	write(filename,'(a,"forProDiMo.fits.gz")') outdir(1:len_trim(outdir))

	inquire(file=filename,exist=truefalse)
	if(truefalse) then
		write(*,'("FITS file already exists, overwriting")')
		write(9,'("FITS file already exists, overwriting")')
		open(unit=90,file=filename)
		close(unit=90,status='delete')
	endif


	do j0=1,D%nTheta-2
		do i=1,D%nR-1
			if(C(i,j0)%dens.gt.1d-40) goto 1
		enddo
	enddo
1	continue
	if(j0.gt.1) then
		write(*,'("Starting from theta: ",f4.2," (cell nr.",i3,")")') D%theta_av(j0),j0
		write(9,'("Starting from theta: ",f4.2," (cell nr.",i3,")")') D%theta_av(j0),j0
	endif

	allocate(grid(D%nR-1,D%nTheta-j0,2))
	allocate(temperature(D%nR-1,D%nTheta-j0,2))
	allocate(spectre(nlam))
	allocate(J_io(D%nR-1,D%nTheta-j0,nlam))
	allocate(dens(D%nR-1,D%nTheta-j0))
	allocate(opacite(D%nR-1,D%nTheta-j0,2,nlam))
	allocate(N_grains(D%nR-1,D%nTheta-j0,0:3))
	allocate(region_index(D%nR-1))
	allocate(HRspec(nlamHR,2))
	


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

	!------------------------------------------------------------------------------
	! HDU 1 : Grille
	!------------------------------------------------------------------------------
	bitpix=-64
	naxis=3
	naxes(1)=D%nR-1			!n_rad
	naxes(2)=D%nTheta-j0	!nz
	naxes(3)=2
	nelements=naxes(1)*naxes(2)*naxes(3)

	Mdust=0e0
	do ri=1,D%nR-1
	do zj=j0,D%nTheta-1
		Mdust=Mdust+C(ri,zj)%mass
	enddo
	enddo
	Mdust=Mdust/Msun
	
	allocate(spec(nlam))
	do l=1,nlam
		if(lam(l).le.0.25.and.lam(l).ge.0.091) then
			spec(l)=D%Fstar(l)
		else
			spec(l)=0d0
		endif
	enddo
	tot=0d0
	fUV=0d0
	call integrate(D%Fstar,tot)
	call integrate(spec,fUV)
	fUV=fUV/tot
	spec(1)=D%Fstar(1)
	spec(3)=lam(1)
	do l=2,nlam-1
		if(lam(l).ge.0.091.and.lam(l-1).le.0.091) then
			spec(1)=D%Fstar(l)
			spec(3)=lam(l)
		endif
		if(lam(l).le.0.250.and.lam(l+1).ge.0.250) then
			spec(2)=D%Fstar(l)
			spec(4)=lam(l)
		endif
	enddo
	pUV=log10(spec(2)/spec(1))/log10(spec(4)/spec(3))
	deallocate(spec)

c	call RenormalizeLRF()

	! Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Write optional keywords to the header
	call VersionDateTime(version)
	call ftpkys(unit,'MCMax',trim(version),'',status)
c	call ftpkys(unit,'mcfost_id',sha_id,'',status)
	call ftpkyj(unit,'mcfost2prodimo',104,'',status)
c	call ftpkys(unit,'mcfost_model_name',trim(para),'',status)

	call ftpkye(unit,'Teff',real(D%Tstar),-8,'[K]',status)
	call ftpkye(unit,'Rstar',real(D%Rstar/Rsun),-8,'[Rsun]',status)
	call ftpkye(unit,'Mstar',real(D%Mstar/Msun),-8,'[Msun]',status)
	call ftpkye(unit,'fUV',real(fUV),-3,'',status)
	call ftpkye(unit,'slope_UV',real(pUV),-3,'',status)
	call ftpkye(unit,'distance',real(D%distance/parsec),-8,'[pc]',status)
	call ftpkye(unit,'edge',real(D%Rout),-8,'[AU]',status)

	call ftpkye(unit,'disk_dust_mass_tot',real(Mdust),-8,'[Msun]',status)

	call ftpkyj(unit,'n_zones',nzones,'',status)
	call ftpkyj(unit,'n_regions',nzones,'',status)

	if(nzones.le.0) then
		call ftpkye(unit,'disk_dust_mass',real(Mdust),-8,'[Msun]',status)
		call ftpkye(unit,'Rin',real(D%R(1)),-8,'[AU]',status)
		call ftpkye(unit,'Rout',real(D%R(D%nR)),-8,'[AU]',status)
		call ftpkye(unit,'Rref',real(1.0),-8,'[AU]',status)
		call ftpkye(unit,'H0',real(D%sh1au/sqrt(2.0)),-8,'[AU]',status)
		call ftpkye(unit,'beta',real(D%shpow),-8,'',status)
		call ftpkye(unit,'alpha',real(-D%denspow),-8,'',status)
	else
		do i=1,nzones
			if(i.eq.1) then
				s=' '
			else
				write(s,'("_",i1)') i
			endif
			call ftpkye(unit,'disk_dust_mass'//trim(s),real(Zone(i)%Mdust),-8,'[Msun]',status)
			call ftpkye(unit,'Rin'//trim(s),real(Zone(i)%Rin),-8,'[AU]',status)
			call ftpkye(unit,'Rout'//trim(s),real(Zone(i)%Rout),-8,'[AU]',status)
			call ftpkye(unit,'Rref'//trim(s),real(Zone(i)%Rsh),-8,'[AU]',status)
			call ftpkye(unit,'H0'//trim(s),real(Zone(i)%sh/sqrt(2.0)),-8,'[AU]',status)
			call ftpkye(unit,'beta'//trim(s),real(Zone(i)%shpow),-8,'',status)
			call ftpkye(unit,'alpha'//trim(s),real(-Zone(i)%denspow),-8,'',status)

			write(s,'("_",i1)') i
			call ftpkye(unit,'Rmin_region'//trim(s),real(Zone(i)%Rin),-8,'[AU]',status)
			call ftpkye(unit,'Rmax_region'//trim(s),real(Zone(i)%Rout),-8,'[AU]',status)
		enddo
	endif

	call ftpkye(unit,'amin',real(mrn_rmin*1d4),-8,'[micron]',status)
	call ftpkye(unit,'amax',real(mrn_rmax*1d4),-8,'[micron]',status)
	call ftpkye(unit,'aexp',real(mrn_index),-8,'slope of grain size distribution',status)
	call ftpkye(unit,'strat',real(0.0),-8,'stratification exponent',status)
	call ftpkye(unit,'a_settle',real(0.0),-8,'[micron]',status)
	call ftpkye(unit,'rho_grain',real(Grain(1)%rho),-8,'[g.cm^-3]',status)
	call ftpkys(unit,'optical_indices','DHS','',status)

	call ftpkyj(unit,'n_grains',ngrains,' ',status)
	call ftpkyj(unit,'n_rad',D%nR-1,' ',status)
	call ftpkyj(unit,'nz',D%nTheta-j0,' ',status)
	call ftpkyj(unit,'n_rad_in',nspan*nlev,' ',status)

	!  Write the array to the FITS file.
	do ri=1, D%nR-1	!n_rad
		if(ri.eq.1) then
		   grid(ri,:,1) = D%R(1)	!r_grid(ri,:)
		else if (ri.eq.D%nR-1) then
		   grid(ri,:,1) = D%R(D%nR)	!r_grid(ri,:)
		else
		   grid(ri,:,1) = D%R_av(ri)/AU	!r_grid(ri,:)
		endif
	   do zj=j0,D%nTheta-1 !nz
		  grid(ri,D%nTheta-zj,2) = grid(ri,D%nTheta-zj,1)*cos(D%theta_av(zj)) !z_grid(ri,zj)
	   enddo
	enddo
	
	call ftpprd(unit,group,fpixel,nelements,grid,status)

	!------------------------------------------------------------------------------
	! HDU 2: Temperature 
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=2
	naxes(1)=D%nR-1			!n_rad
	naxes(2)=D%nTheta-j0		!nz
	nelements=naxes(1)*naxes(2)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	   do ri=1, D%nR-1	!n_rad
		  do zj=j0, D%nTheta-1	!nz
			 Temperature(ri,D%nTheta-zj,1) = C(ri,zj)%T
		  enddo !j
	   enddo !i

	!  Write the array to the FITS file.
	call ftppre(unit,group,fpixel,nelements,Temperature(:,:,1),status)

	!------------------------------------------------------------------------------
	! HDU 3 : Longueurs d'onde
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=1
	naxes(1)=nlam
	nelements=naxes(1)

	! create new hdu
	call ftcrhd(unit, status)

	! Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Write the array to the FITS file.
	call ftppre(unit,group,fpixel,nelements,real(lam),status)

	!------------------------------------------------------------------------------
	! HDU 4 : Spectre stellaire
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=1
	naxes(1)=nlam
	nelements=naxes(1)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Verif Spectre stellaire
	!write(*,*) sum(spectre_etoiles) / (sigma * etoile(1)%T**4 * 4 * pi * (etoile(1)%r * AU_to_m)**2 )

	! Conversion en lambda.F_lamda
	! Division par 4 * pi * (etoile(1)%r * AU_to_m)**2 pour passer de luminosite a flux
	! Division par pi pour passer en intensite
c	spectre(:) = spectre_etoiles(:) * tab_lambda(:) / tab_delta_lambda(:) &
c		 / (4 * pi * (etoile(1)%r * AU_to_m)**2) / pi

	spectre(:) = D%Fstar(:) * 1d-3 * 2.998e14/lam(:) /(pi*D%Rstar**2)
	 
	!  Write the array to the FITS file.
	call ftppre(unit,group,fpixel,nelements,spectre,status)

	!------------------------------------------------------------------------------
	! HDU 5 : Spectre ISM (input)
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=1
	naxes(1)=nlam
	nelements=naxes(1)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do lambda=1,nlam
c	   wl = tab_lambda(lambda) * 1e-6
		if(use_irf) then
		   spectre(lambda) = IRF(lambda) * 1d-3 * 2.998e14/lam(lambda) /(pi*(D%R(D%nR)*AU)**2)  
		else
			spectre(lambda) = 0d0
		endif
	enddo
 
	!  Write the array to the FITS file.
	call ftppre(unit,group,fpixel,nelements,spectre,status)

	!------------------------------------------------------------------------------
	! HDU 6 : Champ de radiation en W.m-2 (lambda.F_lambda)  
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=3
	naxes(1)=D%nR-1
	naxes(2)=D%nTheta-j0
	naxes(3)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)

	call CreateLRF(NphotUV,10000)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Inversion de l'ordre des dimensions + passage en simple precision
	do ri=1, D%nR-1
	   do zj=j0,D%nTheta-1
		  J_io(ri,D%nTheta-zj,:) = 1d-3 * C(ri,zj)%LRF(:) * 2.998e14/lam(:)
	   enddo
	enddo

	call ftppre(unit,group,fpixel,nelements,J_io,status)

	!------------------------------------------------------------------------------
	! HDU 7 : Statistique du champ de radiation (nombre de paquet)
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=3
	naxes(1)=D%nR-1
	naxes(2)=D%nTheta-j0
	naxes(3)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	
	! Inversion de l'ordre des dimensions 
	do zj=j0,D%nTheta-1
	   do ri=1, D%nR-1
		  J_io(ri,D%nTheta-zj,:) = C(ri,zj)%nLRF(:) 
	   enddo
	enddo

	call ftppre(unit,group,fpixel,nelements,J_io,status)

	!------------------------------------------------------------------------------
	! HDU 8 : Champ de radiation ISM en W.m-2 (lambda.F_lambda)  
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=3
	naxes(1)=D%nR-1
	naxes(2)=D%nTheta-j0
	naxes(3)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)


	! Inversion de l'ordre des dimensions + passage en simple precision
	do ri=1, D%nR-1
	   do zj=j0,D%nTheta-1
		  J_io(ri,D%nTheta-zj,:) = 0d0!C(:,ri,zj)%LRF(:) 
	   enddo
	enddo

	call ftppre(unit,group,fpixel,nelements,J_io,status)

	!------------------------------------------------------------------------------
	! HDU 9 : Statistique du champ de radiation ISM (nombre de paquet)
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=3
	naxes(1)=D%nR-1
	naxes(2)=D%nTheta-j0
	naxes(3)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
   
	! Inversion de l'ordre des dimensions et somaation
	do zj=j0,D%nTheta-1
	   do ri=1, D%nR-1
		  J_io(ri,D%nTheta-zj,:) = 0
	   enddo
	enddo

	call ftppre(unit,group,fpixel,nelements,J_io,status)


	!------------------------------------------------------------------------------
	! HDU 10 : Densite de gaz pour un rapport de masse de 100 / poussiere
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=2
	naxes(1)=D%nR-1
	naxes(2)=D%nTheta-j0
	nelements=naxes(1)*naxes(2)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	!  Write the array to the FITS file.
	do ri=1, D%nR-1
	   do zj=j0,D%nTheta-1
		  dens(ri,D%nTheta-zj) =  C(ri,zj)%gasdens*gas2dust
	   enddo
	enddo

	call ftppre(unit,group,fpixel,nelements,dens,status)

	!------------------------------------------------------------------------------
	! HDU 11 : Opacites
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=4
	naxes(1)=D%nR-1
	naxes(2)=D%nTheta-j0
	naxes(3)=2
	naxes(4)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do zj=1,D%nTheta-j0
	   do ri=1,D%nR-1
			 do lambda=1,nlam
			  opacite(ri,zj,1,lambda) = 0d0
			  opacite(ri,zj,2,lambda) = 0d0
				do l=1,ngrains
					do iopac=1,Grain(l)%nopac
					   opacite(ri,zj,1,lambda) = opacite(ri,zj,1,lambda) + C(ri,D%nTheta-zj)%w(l)*C(ri,D%nTheta-zj)%wopac(l,iopac)
     &													*C(ri,D%nTheta-zj)%dens*Grain(l)%Kext(iopac,lambda)
					   opacite(ri,zj,2,lambda) = opacite(ri,zj,2,lambda) + C(ri,D%nTheta-zj)%w(l)*C(ri,D%nTheta-zj)%wopac(l,iopac)
     &													*C(ri,D%nTheta-zj)%dens*Grain(l)%Kabs(iopac,lambda)
					enddo
				enddo ! k
			 enddo ! lambda
	   enddo ! ri
	enddo !zj	
	opacite=opacite*AU

	call ftppre(unit,group,fpixel,nelements,opacite,status)

	!------------------------------------------------------------------------------
	! HDU 12 : Moments de la distribution en tailles des grains
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=3
	naxes(1)=D%nR-1
	naxes(2)=D%nTheta-j0
	naxes(3)=4
	nelements=naxes(1)*naxes(2)*naxes(3)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do zj=1,D%nTheta-j0
	   do ri=1,D%nR-1
		  ! Nbre total de grain : le da est deja dans densite_pouss
			N=0d0
		  do l=1,ngrains
		  	N = N + 1d6*C(ri,D%nTheta-zj)%dens*C(ri,D%nTheta-zj)%w(l)
     &						/(4d0*pi*(Grain(l)%dust_moment3*1d-12)*Grain(l)%rho/3d0)
		  enddo
		  N_grains(ri,zj,0) = N
		  if (N.gt.1d-40) then
			N1=0d0
			do l=1,ngrains
			  	N1 = N1 + 1d6*Grain(l)%dust_moment1*C(ri,D%nTheta-zj)%dens*C(ri,D%nTheta-zj)%w(l)
     &						/(4d0*pi*(Grain(l)%dust_moment3*1d-12)*Grain(l)%rho/3d0)
			enddo
			N_grains(ri,zj,1) = N1 / N

			N1=0d0
			do l=1,ngrains
			  	N1 = N1 + 1d6*Grain(l)%dust_moment2*C(ri,D%nTheta-zj)%dens*C(ri,D%nTheta-zj)%w(l)
     &						/(4d0*pi*(Grain(l)%dust_moment3*1d-12)*Grain(l)%rho/3d0)
			enddo
			N_grains(ri,zj,2) = N1 / N

			N1=0d0
			do l=1,ngrains
			  	N1 = N1 + 1d6*Grain(l)%dust_moment3*C(ri,D%nTheta-zj)%dens*C(ri,D%nTheta-zj)%w(l)
     &						/(4d0*pi*(Grain(l)%dust_moment3*1d-12)*Grain(l)%rho/3d0)
			enddo
			N_grains(ri,zj,3) = N1 / N
		  else
		  	N=0.0
			N_grains(ri,zj,1) = 0.0
			N_grains(ri,zj,2) = 0.0
			N_grains(ri,zj,3) = 0.0
		  endif
	   enddo
	enddo

	call ftppre(unit,group,fpixel,nelements,N_grains,status)


	!------------------------------------------------------------------------------
	! HDU 13 : Region index
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=1
	naxes(1)=D%nR-1
	nelements=naxes(1)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	
	if(nzones.le.1) then
		region_index=1
	else
		region_index=0
		do ri=1, D%nR-1
			do j=1,nzones
				if((D%R_av(ri)/AU).gt.Zone(j)%Rin.and.(D%R_av(ri)/AU).lt.Zone(j)%Rout) then
					region_index(ri)=j
				endif
			enddo
		enddo
	endif

	call ftppre(unit,group,fpixel,nelements,region_index,status)

	!------------------------------------------------------------------------------
	! HDU 14 : High resolution stellar spectrum
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=2
	naxes(1)=nlamHR
	naxes(2)=2
	nelements=naxes(1)*naxes(2)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	
	do lambda=1,nlamHR
		HRspec(lambda,1)=lamHR(lambda)
		HRspec(lambda,2)=FstarHR(lambda) * 1d-3 * 2.998e14/lamHR(lambda) /(pi*D%Rstar**2)
	enddo

	call ftppre(unit,group,fpixel,nelements,HRspec,status)

	!------------------------------------------------------------------------------


	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	if (status.gt.0) then
	   print*,'error in ExportProDiMo'
	   stop
	end if


	deallocate(grid)
	deallocate(temperature)
	deallocate(spectre)
	deallocate(J_io)
	deallocate(dens)
	deallocate(opacite)
	deallocate(N_grains)
	deallocate(region_index)
	deallocate(HRspec)


	return
	end



	subroutine RenormalizeLRF()
	use Parameters
	IMPLICIT NONE
	integer i,j,iT1,iT2,ii,iopac,ilam
	real*8 kp,epsT1,epsT2
	real*8,allocatable :: Kabs(:),spec(:)
	
	allocate(Kabs(nlam))
	allocate(spec(nlam))

	do i=1,D%nR-1
	do j=1,D%nTheta-1
		iT1=C(i,j)%T/dT
		iT2=iT1+1
		epsT1=1d0-(C(i,j)%T-real(iT1)*dT)
		epsT2=1d0-epsT1	
		kp=0d0
		Kabs(1:nlam)=0d0
		do ii=1,ngrains
			if(.not.Grain(ii)%qhp) then
				do iopac=1,Grain(ii)%nopac
					kp=kp+((epsT1*Grain(ii)%Kp(iopac,iT1)+epsT2*Grain(ii)%Kp(iopac,iT2))
     &						*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac))
					do ilam=1,nlam
						Kabs(ilam)=Kabs(ilam)+(Grain(ii)%Kabs(iopac,ilam)
     &						*C(i,j)%w(ii)*C(i,j)%wopac(ii,iopac))
					enddo
				enddo
			endif
		enddo
		C(i,j)%EJv=kp
		spec(1:nlam)=C(i,j)%LRF(1:nlam)*Kabs(1:nlam)
		call integrate(spec,kp)
		C(i,j)%LRF(1:nlam)=C(i,j)%LRF(1:nlam)*C(i,j)%EJv/kp
	enddo
	enddo

	deallocate(Kabs)
	deallocate(spec)
	return
	end
	

	subroutine DoRunProDiMo()
	use Parameters
	IMPLICIT NONE
	character*500 command
	
	command='cp ' // trim(ProDiModir) // '/* ' // trim(outdir)
	
	call system(command)
	
	command='cd ' // trim(outdir) // '; prodimo'
	
	call system(command)
	
	return
	end
	


