	subroutine ImportProdimo(filename)
	use Parameters
	IMPLICIT NONE
	character*500 filename,line,key,value
	character*50 dummy
	real*8 x,z,ran2,theta,r,tot
	real*8,allocatable :: xProDiMo(:,:),zProDiMo(:,:),Tgas(:,:),rhodust(:,:),rhogas(:,:)
	integer i,j,nint,k,ncell,iint,ix,iz
	integer Nx_ProDiMo,Nz_ProDiMo,Nsp_ProDiMo,Nh_ProDiMo,Nc_ProDiMo,Nlam_ProDiMo
	real*8 dminx1,dminx2,zz1,zz2,dminz1,dminz2,Mdisk,Mgas
	integer izmin,ixmin
	
	nint=50
	
	Mdisk=-1d0
	open(unit=80,file=filename,RECL=6000)
	do k=1,22
		read(80,'(a500)') line
		key=line(1:index(line,'=')-1)
		value=line(index(line,'=')+1:len_trim(line))
		do i=1,len_trim(key)-4
			if(key(i:i+4).eq.'Mdisk') then
				read(value,*) Mdisk
			endif
		enddo
	enddo
	read(value,*) Nx_ProDiMo,Nz_ProDiMo,Nsp_ProDiMo,Nh_ProDiMo,Nc_ProDiMo,Nlam_ProDiMo
	read(80,*)
	read(80,*)
	allocate(xProDiMo(Nx_ProDiMo,Nz_ProDiMo))
	allocate(zProDiMo(Nx_ProDiMo,Nz_ProDiMo))
	allocate(rhodust(Nx_ProDiMo,Nz_ProDiMo))
	allocate(rhogas(Nx_ProDiMo,Nz_ProDiMo))
	allocate(Tgas(Nx_ProDiMo,Nz_ProDiMo))
	ncell=Nx_ProDiMo*Nz_ProDiMo
	do k=1,ncell
		read(80,*) ix,iz,xProDiMo(ix,iz),zProDiMo(ix,iz),dummy,dummy,dummy,dummy,dummy,Tgas(ix,iz),dummy,dummy,rhogas(ix,iz)
c     &  ,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
c     &  ,(dummy,i=1,Nh_ProDiMo),(dummy,i=1,Nc_ProDiMo),dummy,dummy,(dummy,i=1,Nsp_ProDiMo),(dummy,i=1,Nlam_ProDiMo)
c     &  ,dummy,dummy,dummy,dummy,rhodust(ix,iz)
		rhodust(ix,iz)=1d0/gas2dust
		rhodust(ix,iz)=rhogas(ix,iz)*rhodust(ix,iz)
		rhodust(ix,iz)=log10(rhodust(ix,iz))
		rhogas(ix,iz)=log10(rhogas(ix,iz))
	enddo
	close(unit=80)
	
	Mgas=0d0
	do i=1,D%nR-1
	call tellertje(i,D%nR-1)
	do j=1,D%nTheta-1
		C(i,j)%dens=1d-50
		do iint=1,nint
			r=D%R(i)+(D%R(i+1)-D%R(i))*ran2(idum)
			theta=D%thet(j)+(D%thet(j+1)-D%thet(j))*ran2(idum)
			x=r*sin(theta)
			z=r*cos(theta)
			iz=1
			if(x.lt.xProDiMo(1,iz)) ixmin=1
			if(x.gt.xProDiMo(Nx_ProDiMo,iz)) ixmin=Nx_ProDiMo-1
			do ix=1,Nx_ProDiMo-1
				if(x.ge.xProDiMo(ix,iz).and.x.le.xProDiMo(ix+1,iz)) then
					ixmin=ix
				endif
			enddo
			izmin=Nz_ProDiMo-1
			do iz=1,Nz_ProDiMo-1
				dminx1=abs(x-xProDiMo(ixmin,iz))
				dminx2=abs(x-xProDiMo(ixmin+1,iz))
				zz1=(zProDiMo(ixmin,iz)*dminx2+zProDiMo(ixmin+1,iz)*dminx1)/(dminx1+dminx2)
				zz2=(zProDiMo(ixmin,iz+1)*dminx2+zProDiMo(ixmin+1,iz+1)*dminx1)/(dminx1+dminx2)
				if(z.ge.zz1.and.z.le.zz2) then
					izmin=iz
				endif
			enddo
			dminx1=abs(x-xProDiMo(ixmin,izmin))
			dminx2=abs(x-xProDiMo(ixmin+1,izmin))
			zz1=(zProDiMo(ixmin,iz)*dminx2+zProDiMo(ixmin+1,iz)*dminx1)/(dminx1+dminx2)
			zz2=(zProDiMo(ixmin,iz+1)*dminx2+zProDiMo(ixmin+1,iz+1)*dminx1)/(dminx1+dminx2)
			dminz1=abs(z-zz1)
			dminz2=abs(z-zz2)

			C(i,j)%dens=C(i,j)%dens+dminx2*dminz2*rhodust(ixmin,izmin)/((dminx1+dminx2)*(dminz1+dminz2))
			C(i,j)%dens=C(i,j)%dens+dminx2*dminz1*rhodust(ixmin,izmin+1)/((dminx1+dminx2)*(dminz1+dminz2))
			C(i,j)%dens=C(i,j)%dens+dminx1*dminz2*rhodust(ixmin+1,izmin)/((dminx1+dminx2)*(dminz1+dminz2))
			C(i,j)%dens=C(i,j)%dens+dminx1*dminz1*rhodust(ixmin+1,izmin+1)/((dminx1+dminx2)*(dminz1+dminz2))

			C(i,j)%gasdens=C(i,j)%gasdens+dminx2*dminz2*rhogas(ixmin,izmin)/((dminx1+dminx2)*(dminz1+dminz2))
			C(i,j)%gasdens=C(i,j)%gasdens+dminx2*dminz1*rhogas(ixmin,izmin+1)/((dminx1+dminx2)*(dminz1+dminz2))
			C(i,j)%gasdens=C(i,j)%gasdens+dminx1*dminz2*rhogas(ixmin+1,izmin)/((dminx1+dminx2)*(dminz1+dminz2))
			C(i,j)%gasdens=C(i,j)%gasdens+dminx1*dminz1*rhogas(ixmin+1,izmin+1)/((dminx1+dminx2)*(dminz1+dminz2))

			C(i,j)%Tgas=C(i,j)%Tgas+dminx2*dminz2*Tgas(ixmin,izmin)/((dminx1+dminx2)*(dminz1+dminz2))
			C(i,j)%Tgas=C(i,j)%Tgas+dminx2*dminz1*Tgas(ixmin,izmin+1)/((dminx1+dminx2)*(dminz1+dminz2))
			C(i,j)%Tgas=C(i,j)%Tgas+dminx1*dminz2*Tgas(ixmin+1,izmin)/((dminx1+dminx2)*(dminz1+dminz2))
			C(i,j)%Tgas=C(i,j)%Tgas+dminx1*dminz1*Tgas(ixmin+1,izmin+1)/((dminx1+dminx2)*(dminz1+dminz2))
		enddo
		C(i,j)%dens=10d0**(C(i,j)%dens/real(nint))
		C(i,j)%gasdens=10d0**(C(i,j)%gasdens/real(nint))/gas2dust
		C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
		Mgas=Mgas+C(i,j)%gasdens*C(i,j)%V*gas2dust
		C(i,j)%Tgas=C(i,j)%Tgas/real(nint)
	enddo
	enddo
	Mgas=Mgas/Msun
	print*,Mgas,Mdisk

	if(Mdisk.gt.0d0) then
		do i=1,D%nR-1
		do j=1,D%nTheta-1
			C(i,j)%dens=C(i,j)%dens*Mdisk/Mgas
			C(i,j)%gasdens=C(i,j)%gasdens*Mdisk/Mgas
		enddo
		enddo
	endif
	
	deallocate(xProDiMo)
	deallocate(zProDiMo)
	deallocate(rhodust)
	deallocate(rhogas)
	deallocate(Tgas)


	return
	end
	
	
	subroutine PopulationsProdimo(datafile,i_low,i_up,n_low,n_up,nm,mol_mass,poptype)
	use Parameters
	IMPLICIT NONE
	character*500 datafile,line,key,value
	character*50 dummy
	character*10 poptype
	real*8 x,z,ran2,theta,r,tot,n_up(0:D%nR-1,D%nTheta-1),n_low(0:D%nR-1,D%nTheta-1)
	real*8 nm(0:D%nR-1,D%nTheta-1),mol_mass,mp
	real*8,allocatable :: x1(:),x2(:),z1(:),z2(:),up(:),low(:),nm_in(:)
	integer i,j,nint,k,ncell,iint,i_low,i_up
	parameter(mp=1.67262158d-24) ! the proton mass in gram
	real*8 nl(i_up),u0,l0
	
	nint=10

	poptype='ProDiMo'
	
	open(unit=80,file=datafile,RECL=6000)
	do k=1,15
		read(80,'(a500)') line
		key=line(1:index(line,'=')-1)
		value=line(index(line,'=')+1:len_trim(line))
		if(key.eq.'ncell') read(value,*) ncell
	enddo
	allocate(x1(ncell))
	allocate(x2(ncell))
	allocate(z1(ncell))
	allocate(z2(ncell))
	allocate(up(ncell))
	allocate(low(ncell))
	allocate(nm_in(ncell))

	do k=1,ncell
		read(80,*,end=2) j,x1(k),x2(k),z1(k),z2(k),dummy,dummy,nm_in(k),(dummy,i=1,5)
     &		,(nl(i),i=1,i_up)
		low(k)=nl(i_low)
		up(k)=nl(i_up)
	enddo
	goto 3

2	continue
	write(*,'("Upper level not found. Using LTE")')
	write(9,'("Upper level not found. Using LTE")')
	poptype='LTE'
	close(unit=80)
	open(unit=80,file=datafile,RECL=6000)
	do k=1,15
		read(80,'(a500)') 
	enddo
	do k=1,ncell
		read(80,*) j,x1(k),x2(k),z1(k),z2(k),dummy,dummy,nm_in(k)
	enddo

3	continue
	close(unit=80)
	
	do i=1,D%nR-1
	call tellertje(i,D%nR-1)
	do j=1,D%nTheta-1
		if(poptype.ne.'LTE') then
			n_up(i,j)=0d0
			n_low(i,j)=0d0
		endif
		nm(i,j)=0d0
		do iint=1,nint
			r=D%R(i)+(D%R(i+1)-D%R(i))*ran2(idum)
			theta=D%thet(j)+(D%thet(j+1)-D%thet(j))*ran2(idum)
			x=r*sin(theta)*AU*0.01
			z=r*cos(theta)*AU*0.01
			do k=1,ncell
				if(x.ge.x1(k).and.x.le.x2(k).and.z.ge.z1(k).and.z.le.z2(k)) then
					if(poptype.ne.'LTE') then
						n_up(i,j)=n_up(i,j)+(up(k))
						n_low(i,j)=n_low(i,j)+(low(k))
					endif
					nm(i,j)=nm(i,j)+nm_in(k)
					goto 1
				endif
			enddo
1			continue
		enddo
		if(poptype.ne.'LTE') then
			n_up(i,j)=n_up(i,j)/real(nint)
			n_low(i,j)=n_low(i,j)/real(nint)
		endif
		nm(i,j)=mp*mol_mass*nm(i,j)/real(nint)/(C(i,j)%gasdens*gas2dust)
	enddo
	enddo

	deallocate(x1)
	deallocate(x2)
	deallocate(z1)
	deallocate(z2)
	deallocate(up)
	deallocate(low)
	deallocate(nm_in)

	return
	end


