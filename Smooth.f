	subroutine do_smooth()
	use Parameters
	IMPLICIT NONE
	real*8,allocatable :: T1(:),T2(:),dT1(:)
	integer i,j
	
	allocate(T1(D%nR-1))
	allocate(dT1(D%nR-1))
	allocate(T2(D%nR-1))
	do j=1,D%nTheta-1
		do i=1,D%nR-1
			T1(i)=C(i,j)%T
			dT1(i)=C(i,j)%dT
		enddo
		call Smooth(T1,dT1,T2,D%nR-1)
		do i=1,D%nR-1
			if(T2(i).gt.dT) C(i,j)%T=T2(i)
		enddo
	enddo
	deallocate(T1)
	deallocate(dT1)
	deallocate(T2)

	return
	end


	subroutine Smooth(y,dy,z,n)
	IMPLICIT NONE
	integer m,n,MAXM,MAXITER
	parameter(MAXM=1000,MAXITER=100)
	real*8 y(n),dy(n),x(n),a(MAXM,n),b(MAXM),min,z(n),maxsigma
	real*8 rnorm,w(MAXM*3),beta(MAXITER),error(MAXITER),eps
	integer info,i,j,niter,ii
	m=2*n-1

	eps=1d-2
	beta(1)=0.1d0

	maxsigma=1d0

	niter=1
1	continue
	do i=1,n
		do j=1,n
			a(i,j)=0d0
			a(i+n,j)=0d0
		enddo
		a(i,i)=1d0/dy(i)
		b(i)=y(i)/dy(i)
		if(i.lt.n.and.i.gt.1) then
			a(i+n,i-1)=-beta(niter)
			a(i+n,i)=2d0*beta(niter)
			a(i+n,i+1)=-beta(niter)
			b(i+n)=0d0
		endif
	enddo

	call dgels('N',m,n,1,a,MAXM,b,MAXM,w,MAXM*3,info)

	error(niter)=0d0
	do i=1,n
		error(niter)=error(niter)+((y(i)-b(i))/dy(i))**2
	enddo
	error(niter)=error(niter)/real(n-1)
	
	if(abs(error(niter)-maxsigma).gt.eps) then
		niter=niter+1
		if(niter.gt.MAXITER) goto 2
		beta(niter)=beta(niter-1)/sqrt(error(niter-1))
		goto 1
	endif

	goto 3

2	min=0d0
	ii=MAXITER
	do i=1,MAXITER
		if(error(i).gt.min.and.error(i).lt.1d0) then
			min=error(i)
			ii=i
		endif
	enddo
	
	do i=1,n
		do j=1,n
			a(i,j)=0d0
			a(i+n,j)=0d0
		enddo
		a(i,i)=1d0/dy(i)
		b(i)=y(i)/dy(i)
		if(i.lt.n.and.i.gt.1) then
			a(i+n,i-1)=-beta(ii)
			a(i+n,i)=2d0*beta(ii)
			a(i+n,i+1)=-beta(ii)
			b(i+n)=0d0
		endif
	enddo
	call dgels('N',m,n,1,a,MAXM,b,MAXM,w,MAXM*3,info)
	
3	continue
	do i=1,n
		z(i)=b(i)
	enddo

	return
	end

c ---------------------------------------------------------------------
c A general 2D smoothing routine, using kernel:
c
c  (1,2,1)
c  (2,4,2) / 16
c  (1,2,1)
c
c ---------------------------------------------------------------------

	subroutine LogSmooth2D(array,nx,ny)
	implicit none
	integer i,j,nx,ny
	real*8 array(1:nx,1:ny),temp(0:nx+1,0:ny+1)

	! Needs to be positive
	if(minval(array).le.0d0) then
	   print*,"LogSmooth2D: no log smoothing with negative or zero values!!"
	   stop
	endif

	!  Copy array into larger array
c$$$	temp(1:nx  ,1:ny)=log(array(1:nx  ,1:ny)) ! core array
c$$$	temp(0     ,1:ny)=log(array(1     ,1:ny)) ! left
c$$$	temp(nx+1  ,1:ny)=log(array(nx    ,1:ny)) ! right
c$$$	temp(0:nx+1,0   )=log( temp(0:nx+1,1   )) ! top+corner
c$$$	temp(0:nx+1,ny+1)=log( temp(0:nx+1,ny  )) ! bottom+corner

	! *sigh...
	do i=1,nx
	   do j=1,ny
	      temp(i,j)=log(array(i,j))	! core array
	   enddo
	enddo

	! array boundaries
	temp(0     ,1:ny)= temp(1     ,1:ny) ! left
	temp(nx+1  ,1:ny)= temp(nx    ,1:ny) ! right
	temp(0:nx+1,0   )= temp(0:nx+1,1   ) ! top+corner
	temp(0:nx+1,ny+1)= temp(0:nx+1,ny  ) ! bottom+corner
	
	!  Smooth with the kernel
	do i=1,nx
	   do j=1,ny
	      array(i,j)=    temp(i-1,j-1)+2d0*temp(i,j-1)+    temp(i+1,j-1)+
     &	                 2d0*temp(i-1,j  )+4d0*temp(i,j  )+2d0*temp(i+1,j  )+
     &	                     temp(i-1,j+1)+2d0*temp(i,j+1)+    temp(i+1,j+1)
	      array(i,j)=exp(array(i,j)/16d0)
	   enddo
	enddo

	end subroutine
