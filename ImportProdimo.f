	subroutine ImportProdimo(filename)
	use Parameters
	IMPLICIT NONE
	character*500 filename,line,key,value
	real*8 x,z,ran2,theta,r,tot,dummy
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
		rhodust(ix,iz)=1d-2
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
		C(i,j)%gasdens=10d0**(C(i,j)%gasdens/real(nint))
		C(i,j)%V=(4d0*pi/3d0)*(D%R(i+1)**3-D%R(i)**3)*
     &			(D%Theta(j)-D%Theta(j+1))*AU**3
		Mgas=Mgas+C(i,j)%gasdens*C(i,j)%V
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
	