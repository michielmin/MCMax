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

