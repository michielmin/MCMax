	subroutine DoExportProdimo()
	use Parameters
	IMPLICIT NONE

	  integer status,unit,blocksize,bitpix,naxis,naxes(2)
	  integer j,group,fpixel,nelements,ri,zj,l,lambda,iopac
	  logical simple,extend,truefalse
	character*500 filename,version
	real*8,allocatable :: grid(:,:,:)
	real,allocatable :: temperature(:,:,:)
	real,allocatable :: spectre(:),J_io(:,:,:),dens(:,:)
	real,allocatable :: opacite(:,:,:,:),N_grains(:,:,:)
	real N,N1,Mdust
	real*8 fUV,tot
	real*8,allocatable :: spec(:)

	allocate(grid(D%nR-1,D%nTheta-1,2))
	allocate(temperature(D%nR-1,D%nTheta-1,2))
	allocate(spectre(nlam))
	allocate(J_io(D%nR-1,D%nTheta-1,nlam))
	allocate(dens(D%nR-1,D%nTheta-1))
	allocate(opacite(D%nR-1,D%nTheta-1,2,nlam))
	allocate(N_grains(D%nR-1,D%nTheta-1,0:3))
	
	write(filename,'(a,"forProDiMo.fits.gz")') outdir(1:len_trim(outdir))

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

	!------------------------------------------------------------------------------
	! HDU 1 : Grille
	!------------------------------------------------------------------------------
	bitpix=-64
	naxis=3
	naxes(1)=D%nR-1			!n_rad
	naxes(2)=D%nTheta-1		!nz
	naxes(3)=2
	nelements=naxes(1)*naxes(2)*naxes(3)

	Mdust=0e0
	do ri=1,D%nR-1
	do zj=1,D%nTheta-1
		Mdust=MDust+C(ri,zj)%mass
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
	print*,fUV
	fUV=fUV/tot
	deallocate(spec)

	! Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Write optional keywords to the header
	call VersionDateTime(version)
	call ftpkys(unit,'MCMax',trim(version),'',status)
c	call ftpkys(unit,'mcfost_id',sha_id,'',status)
c	call ftpkyj(unit,'mcfost2prodimo',mcfost2ProDiMo_version,'',status)
c	call ftpkys(unit,'mcfost_model_name',trim(para),'',status)

	call ftpkye(unit,'Teff',real(D%Tstar),-8,'[K]',status)
	call ftpkye(unit,'Rstar',real(D%Rstar/Rsun),-8,'[Rsun]',status)
	call ftpkye(unit,'Mstar',real(D%Mstar/Msun),-8,'[Msun]',status)
	call ftpkye(unit,'fUV',real(fUV),-3,'',status)
	call ftpkye(unit,'slope_UV',1.0+2.0,-3,'',status)
	call ftpkye(unit,'distance',real(D%distance/parsec),-8,'[pc]',status)
	
	call ftpkye(unit,'disk_dust_mass',real(Mdust),-8,'[Msun]',status)
	call ftpkye(unit,'Rin',real(D%R(1)),-8,'[AU]',status)
	call ftpkye(unit,'Rout',real(D%R(D%nR)),-8,'[AU]',status)
	call ftpkye(unit,'Rref',real(1.0),-8,'[AU]',status)
	call ftpkye(unit,'H0',real(D%sh1au),-8,'[AU]',status)
	call ftpkye(unit,'edge',real(D%Rout),-8,'[AU]',status)
	call ftpkye(unit,'beta',real(D%shpow),-8,'',status)
	call ftpkye(unit,'alpha',real(D%denspow),-8,'',status)

	call ftpkye(unit,'amin',real(mrn_rmin),-8,'[micron]',status)
	call ftpkye(unit,'amax',real(mrn_rmax),-8,'[micron]',status)
	call ftpkye(unit,'aexp',real(mrn_index),-8,'slope of grain size distribution',status)
	call ftpkye(unit,'strat',real(1.0),-8,'stratification exponent',status)
	call ftpkye(unit,'a_settle',real(1.0),-8,'[micron]',status)
	call ftpkye(unit,'rho_grain',real(Grain(1)%rho),-8,'[g.cm^-3]',status)
	call ftpkys(unit,'optical_indices','DHS','',status)

	call ftpkyj(unit,'n_grains',ngrains,' ',status)
	call ftpkyj(unit,'n_rad',D%nR-1,' ',status)
	call ftpkyj(unit,'nz',D%nTheta-1,' ',status)
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
	   do zj=1,D%nTheta-1 !nz
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
	naxes(2)=D%nTheta-1		!nz
	nelements=naxes(1)*naxes(2)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	   do ri=1, D%nR-1	!n_rad
		  do zj=1, D%nTheta-1	!nz
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
	   spectre(lambda) = 0d0	!(chi_ISM * 1.71 * Wdil * Blambda(wl,T_ISM_stars) + Blambda(wl,TCmb)) * wl  
	enddo
 
	!  Write the array to the FITS file.
	call ftppre(unit,group,fpixel,nelements,spectre,status)

	!------------------------------------------------------------------------------
	! HDU 6 : Champ de radiation en W.m-2 (lambda.F_lambda)  
	!------------------------------------------------------------------------------
	bitpix=-32
	naxis=3
	naxes(1)=D%nR-1
	naxes(2)=D%nTheta-1
	naxes(3)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)

	call CreateLRF(100000,10000)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Inversion de l'ordre des dimensions + passage en simple precision
	do ri=1, D%nR-1
	   do zj=1,D%nTheta-1
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
	naxes(2)=D%nTheta-1
	naxes(3)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	
	! Inversion de l'ordre des dimensions 
	do zj=1,D%nTheta-1
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
	naxes(2)=D%nTheta-1
	naxes(3)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)


	! Inversion de l'ordre des dimensions + passage en simple precision
	do ri=1, D%nR-1
	   do zj=1,D%nTheta-1
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
	naxes(2)=D%nTheta-1
	naxes(3)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
   
	! Inversion de l'ordre des dimensions et somaation
	do zj=1,D%nTheta-1
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
	naxes(2)=D%nTheta-1
	nelements=naxes(1)*naxes(2)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	!  Write the array to the FITS file.
	do ri=1, D%nR-1
	   do zj=1,D%nTheta-1
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
	naxes(2)=D%nTheta-1
	naxes(3)=2
	naxes(4)=nlam
	nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do zj=1,D%nTheta-1
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
	naxes(2)=D%nTheta-1
	naxes(3)=4
	nelements=naxes(1)*naxes(2)*naxes(3)

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do zj=1,D%nTheta-1
	   do ri=1,D%nR-1
		  ! Nbre total de grain : le da est deja dans densite_pouss
			N=0d0
		  do l=1,ngrains
		  	N = N + 1d6*C(ri,D%nTheta-zj)%dens*C(ri,D%nTheta-zj)%w(l)/(4d0*pi*(Grain(l)%rv*1d-4)**3*Grain(l)%rho/3d0)
		  enddo
		  N_grains(ri,zj,0) = N
		  if (N.gt.1d-40) then
			N1=0d0
			do l=1,ngrains
				N1 = N1 + 1d6*Grain(l)%rv*C(ri,D%nTheta-zj)%dens*C(ri,D%nTheta-zj)%w(l)/(4d0*pi*(Grain(l)%rv*1d-4)**3*Grain(l)%rho/3d0)
			enddo
			N_grains(ri,zj,1) = N1 / N

			N1=0d0
			do l=1,ngrains
				N1 = N1 + 1d6*Grain(l)%rv**2*C(ri,D%nTheta-zj)%dens*C(ri,D%nTheta-zj)%w(l)/(4d0*pi*(Grain(l)%rv*1d-4)**3*Grain(l)%rho/3d0)
			enddo
			N_grains(ri,zj,2) = N1 / N

			N1=0d0
			do l=1,ngrains
				N1 = N1 + 1d6*Grain(l)%rv**3*C(ri,D%nTheta-zj)%dens*C(ri,D%nTheta-zj)%w(l)/(4d0*pi*(Grain(l)%rv*1d-4)**3*Grain(l)%rho/3d0)
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


	return
	end

