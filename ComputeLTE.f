	subroutine ComputeLTE(datafile,iTr,i_low,i_up,Aul,Bul,Blu,nu0,n_low,n_up,mol_mass)
	use Parameters
	IMPLICIT NONE
	character*500 datafile
	character*10 name
	integer,allocatable :: u(:),l(:)
	integer i,j,nE,nTr,iTr,k,i_low,i_up
	real*8,allocatable :: g(:),E(:),freq(:),EinstA(:)
	real*8 mol_mass,UT,T,Aul,Bul,Blu,n_up(0:D%nR-1,D%nTheta-1),n_low(0:D%nR-1,D%nTheta-1),h,nu0
	parameter(h=6.626068e-27) ! cm^2 g/s
	
	open(unit=80,file=datafile,RECL=6000)
	read(80,*)
	read(80,*) name
	read(80,*)
	read(80,*) mol_mass
	read(80,*)
	read(80,*) nE
	read(80,*)
	allocate(E(nE))
	allocate(g(nE))
	do i=1,nE
		read(80,*) j,E(i),g(i)
	enddo
	read(80,*)
	read(80,*) nTr
	read(80,*)
	allocate(u(nTr))
	allocate(l(nTr))
	allocate(EinstA(nTr))
	allocate(freq(nTr))
	do i=1,nTr
		read(80,*) j,u(i),l(i),EinstA(i),freq(i),E(u(i))
	enddo
	close(unit=80)

	i_low=l(iTr)
	i_up=u(iTr)
	
	Aul=EinstA(iTr)
	nu0=freq(iTr)*1d9

	Bul=Aul*2d0*h*nu0**3/((2.9979d10)**2)
	Blu=Bul*g(i_up)/g(i_low)

	do i=1,D%nR-1
	do j=1,D%nTheta-1
		if(useTgas) then
			T=C(i,j)%Tgas
		else
			T=C(i,j)%T
		endif
		UT=0d0
		do k=1,nE
			UT=UT+g(k)*exp(-E(k)/T)
		enddo
		n_up(i,j)=g(i_up)*exp(-E(i_up)/T)/UT
		n_low(i,j)=g(i_low)*exp(-E(i_low)/T)/UT
	enddo
	enddo

	deallocate(E)
	deallocate(g)
	deallocate(u)
	deallocate(l)
	deallocate(EinstA)
	deallocate(freq)

	
	return
	end
	

	