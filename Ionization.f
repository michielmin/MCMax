	subroutine PAHionization(ii,i,j)
	use Parameters
	IMPLICIT NONE
	real*8 G,spec(nlam),lam1,lam2,mu,gamma0,ne
	integer i,j,l,ii
	parameter(mu=1.3*1.67262158d-24) !1.3 times the proton mass in gram

	lam1=0.0953
	lam2=0.206
	
	spec(1:nlam)=C(i,j)%LRF(1:nlam)
	do l=1,nlam
		if(lam(l).lt.lam1) spec(l)=0d0
		if(lam(l).gt.lam2) spec(l)=0d0
	enddo
	call integrate(spec,G)

	G=4d0*pi*G/5.33d-14/3d10
	if(G.lt.1d-6) G=1d-6
	
	ne = C(i,j)%gasdens*gas2dust / mu
	
	gamma0 = 3.5e-6 * sqrt(Grain(ii)%Nc) * G * sqrt(C(i,j)%T) / ne

	C(i,j)%wopac(ii,1)=1d0/(1d0+gamma0)
	C(i,j)%wopac(ii,2)=1d0-C(i,j)%wopac(ii,1)
		
	return
	end

