	subroutine CreateLRF(Nphot,NphotStar,justUV)
	use Parameters
	IMPLICIT NONE
	integer Nphot,i,j,k,ilam,NphotStar
	logical old_forcefirst,justUV
	real*8 old_pfstop
	
	makeangledependence=.true.
	storescatt=.false.
	useobspol=.false.
	old_forcefirst=forcefirst
	old_pfstop=pfstop
	
	if(scattering.and.Nphot.ne.0) then
		do i=0,D%nR
		do j=1,D%nTheta-1
			C(i,j)%nrg=1
			if(allocated(C(i,j)%thetrg)) deallocate(C(i,j)%thetrg)
			allocate(C(i,j)%thetrg(C(i,j)%nrg+1))
			do k=1,C(i,j)%nrg+1
				C(i,j)%thetrg(k)=cos(D%thet(j)+real(k-1)*(D%thet(j+1)-D%thet(j))/real(C(i,j)%nrg))**2
			enddo
			if(allocated(C(i,j)%scattfield)) deallocate(C(i,j)%scattfield)
			allocate(C(i,j)%scattfield(C(i,j)%nrg,0:nlam,2))
			if(scat_how.eq.2) then
				if(allocated(C(i,j)%scattQ)) deallocate(C(i,j)%scattQ)
				if(allocated(C(i,j)%scattU)) deallocate(C(i,j)%scattU)
				if(allocated(C(i,j)%scattV)) deallocate(C(i,j)%scattV)
				allocate(C(i,j)%scattQ(C(i,j)%nrg,0:nlam,2))
				allocate(C(i,j)%scattU(C(i,j)%nrg,0:nlam,2))
				allocate(C(i,j)%scattV(C(i,j)%nrg,0:nlam,2))
			endif
		enddo
		enddo
	else
		do i=0,D%nR
		do j=1,D%nTheta-1
			C(i,j)%nrg=1
			if(allocated(C(i,j)%thetrg)) deallocate(C(i,j)%thetrg)
			allocate(C(i,j)%thetrg(C(i,j)%nrg+1))
			do k=1,C(i,j)%nrg+1
				C(i,j)%thetrg(k)=cos(D%thet(j)+real(k-1)*(D%thet(j+1)-D%thet(j))/real(C(i,j)%nrg))**2
			enddo
		enddo
		enddo
	endif

	do ilam=1,nlam
		do i=1,D%nR-1
			do j=1,D%nTheta-1
				C(i,j)%scattfield(1,0,1:2)=0d0
				C(i,j)%Ni=0
				if(scat_how.eq.2) then
					C(i,j)%scattQ(1,0,1:2)=0d0
					C(i,j)%scattU(1,0,1:2)=0d0
					C(i,j)%scattV(1,0,1:2)=0d0
				endif
			enddo
		enddo
		if(lam(ilam).lt.maxlamUV) then
			nexits=0
			forcefirst=.true.
			pfstop=0.25
			call TraceMono(lam(ilam),Nphot,45d0,NphotStar)
		else if(.not.justUV) then
			nexits=0
			forcefirst=.false.
			pfstop=0.25
			call TraceMono(lam(ilam),Nphot/10,45d0,NphotStar)
		endif
		do i=1,D%nR-1
			do j=1,D%nTheta-1
				if((C(i,j)%Ni+C(i,j)%nLRF(ilam)).gt.0) then
					C(i,j)%scattfield(1,ilam,1)=(real(C(i,j)%Ni)*C(i,j)%scattfield(1,0,1)*2d0/C(i,j)%V
     &							+real(C(i,j)%nLRF(ilam))*C(i,j)%LRF(ilam))/real(C(i,j)%Ni+C(i,j)%nLRF(ilam))
				else
					C(i,j)%scattfield(1,ilam,1)=C(i,j)%LRF(ilam)
     			endif
				C(i,j)%LRF(ilam)=C(i,j)%scattfield(1,ilam,1)
				C(i,j)%nLRF(ilam)=C(i,j)%Ni+C(i,j)%nLRF(ilam)
			enddo
		enddo
	enddo

	if(scattering.and.Nphot.ne.0) then
		do i=0,D%nR
		do j=1,D%nTheta-1
			deallocate(C(i,j)%scattfield)
			if(scat_how.eq.2) then
				deallocate(C(i,j)%scattQ)
				deallocate(C(i,j)%scattU)
				deallocate(C(i,j)%scattV)
			endif
		enddo
		enddo
	endif
	do i=0,D%nR
	do j=1,D%nTheta-1
		if(allocated(C(i,j)%thetrg)) deallocate(C(i,j)%thetrg)
	enddo
	enddo

	forcefirst=old_forcefirst
	pfstop=old_pfstop

	return
	end
	


