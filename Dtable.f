	module DiffTable
	IMPLICIT NONE
	integer NALB,NGSYM,i,j,iopac
	parameter(NALB=9,NGSYM=51)
	real*8 a(NALB),g(NGSYM),D(NALB,NGSYM)
 
	data (a(i),i=1,9)  / 0.0000, 0.2000, 0.4000, 0.6000, 0.8000, 0.9000, 0.9500, 0.9900, 1.0000/
	data (g(i),i=1,51) /-1.0000,-0.9750,-0.9500,-0.9250,-0.9000,-0.8750,-0.8500,-0.8250,-0.8000,-0.7750,
     & -0.7500,-0.7000,-0.6500,-0.6000,-0.5500,-0.5000,-0.4500,-0.4000,-0.3500,-0.3000,
     & -0.2500,-0.2000,-0.1500,-0.1000,-0.0500, 0.0000, 0.0500, 0.1000, 0.1500, 0.2000,
     &  0.2500, 0.3000, 0.3500, 0.4000, 0.4500, 0.5000, 0.5500, 0.6000, 0.6500, 0.7000,
     &  0.7500, 0.7750, 0.8000, 0.8250, 0.8500, 0.8750, 0.9000, 0.9250, 0.9500, 0.9750, 1.0000/
	data ((D(i,j),i=1,9),j=1,51) /
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0003, 0.9995, 0.3333,
     &  1.0000, 0.9958, 0.9308, 0.8687, 0.7585, 0.6441, 0.5365, 0.3851, 0.3333,
     &  1.0000, 0.9917, 0.8977, 0.8106, 0.6712, 0.5492, 0.4569, 0.3602, 0.3333,
     &  1.0000, 0.9875, 0.8736, 0.7696, 0.6165, 0.4990, 0.4222, 0.3517, 0.3333,
     &  1.0000, 0.9451, 0.8543, 0.7380, 0.5781, 0.4680, 0.4029, 0.3474, 0.3333,
     &  1.0000, 0.9411, 0.8386, 0.7123, 0.5494, 0.4471, 0.3907, 0.3446, 0.3333,
     &  1.0000, 0.9371, 0.8252, 0.6909, 0.5272, 0.4320, 0.3823, 0.3430, 0.3333,
     &  1.0000, 0.9331, 0.8135, 0.6726, 0.5095, 0.4208, 0.3764, 0.3416, 0.3333,
     &  1.0000, 0.9280, 0.8034, 0.6568, 0.4951, 0.4120, 0.3717, 0.3408, 0.3333,
     &  1.0000, 0.9240, 0.7941, 0.6429, 0.4832, 0.4050, 0.3680, 0.3400, 0.3333,
     &  1.0000, 0.9200, 0.7856, 0.6307, 0.4731, 0.3993, 0.3651, 0.3394, 0.3333,
     &  1.0000, 0.9120, 0.7708, 0.6097, 0.4571, 0.3906, 0.3607, 0.3386, 0.3333,
     &  1.0000, 0.9040, 0.7572, 0.5921, 0.4449, 0.3842, 0.3575, 0.3380, 0.3333,
     &  1.0000, 0.8960, 0.7444, 0.5772, 0.4354, 0.3795, 0.3552, 0.3376, 0.3333,
     &  1.0000, 0.8880, 0.7320, 0.5640, 0.4278, 0.3758, 0.3534, 0.3373, 0.3333,
     &  1.0000, 0.8800, 0.7200, 0.5523, 0.4215, 0.3729, 0.3520, 0.3371, 0.3333,
     &  1.0000, 0.8720, 0.7081, 0.5417, 0.4164, 0.3705, 0.3509, 0.3367, 0.3333,
     &  1.0000, 0.8640, 0.6964, 0.5322, 0.4122, 0.3686, 0.3500, 0.3364, 0.3333,
     &  1.0000, 0.8560, 0.6850, 0.5235, 0.4085, 0.3670, 0.3494, 0.3363, 0.3333,
     &  1.0000, 0.8480, 0.6738, 0.5158, 0.4055, 0.3657, 0.3488, 0.3362, 0.3333,
     &  1.0000, 0.8400, 0.6628, 0.5088, 0.4031, 0.3647, 0.3483, 0.3363, 0.3333,
     &  1.0000, 0.8320, 0.6524, 0.5027, 0.4010, 0.3638, 0.3478, 0.3361, 0.3333,
     &  1.0000, 0.8240, 0.6425, 0.4973, 0.3993, 0.3632, 0.3476, 0.3359, 0.3333,
     &  1.0000, 0.8160, 0.6333, 0.4928, 0.3980, 0.3627, 0.3474, 0.3362, 0.3333,
     &  1.0000, 0.8080, 0.6249, 0.4889, 0.3969, 0.3624, 0.3472, 0.3361, 0.3333,
     &  1.0000, 0.8002, 0.6174, 0.4858, 0.3963, 0.3623, 0.3472, 0.3361, 0.3333,
     &  1.0000, 0.7922, 0.6109, 0.4834, 0.3959, 0.3622, 0.3471, 0.3360, 0.3333,
     &  1.0000, 0.7846, 0.6054, 0.4818, 0.3958, 0.3623, 0.3472, 0.3362, 0.3333,
     &  1.0000, 0.7766, 0.6010, 0.4808, 0.3961, 0.3625, 0.3474, 0.3360, 0.3333,
     &  1.0000, 0.7700, 0.5978, 0.4805, 0.3967, 0.3630, 0.3477, 0.3360, 0.3333,
     &  1.0000, 0.7644, 0.5958, 0.4811, 0.3976, 0.3634, 0.3480, 0.3362, 0.3333,
     &  1.0000, 0.7597, 0.5949, 0.4823, 0.3988, 0.3642, 0.3485, 0.3362, 0.3333,
     &  1.0000, 0.7559, 0.5953, 0.4843, 0.4004, 0.3652, 0.3489, 0.3363, 0.3333,
     &  1.0000, 0.7537, 0.5970, 0.4872, 0.4026, 0.3665, 0.3496, 0.3364, 0.3333,
     &  1.0000, 0.7526, 0.5999, 0.4912, 0.4053, 0.3680, 0.3504, 0.3369, 0.3333,
     &  1.0000, 0.7532, 0.6045, 0.4962, 0.4088, 0.3699, 0.3514, 0.3371, 0.3333,
     &  1.0000, 0.7448, 0.6107, 0.5028, 0.4131, 0.3723, 0.3526, 0.3373, 0.3333,
     &  1.0000, 0.7596, 0.6187, 0.5109, 0.4183, 0.3753, 0.3541, 0.3374, 0.3333,
     &  1.0000, 0.7510, 0.6291, 0.5211, 0.4253, 0.3792, 0.3562, 0.3380, 0.3333,
     &  1.0000, 0.7747, 0.6423, 0.5341, 0.4342, 0.3845, 0.3590, 0.3387, 0.3333,
     &  1.0000, 0.7657, 0.6590, 0.5510, 0.4463, 0.3916, 0.3630, 0.3394, 0.3333,
     &  1.0000, 0.7612, 0.6693, 0.5612, 0.4537, 0.3962, 0.3657, 0.3403, 0.3333,
     &  1.0000, 0.8018, 0.6810, 0.5734, 0.4629, 0.4020, 0.3687, 0.3410, 0.3333,
     &  1.0000, 0.7970, 0.6945, 0.5874, 0.4737, 0.4090, 0.3728, 0.3420, 0.3333,
     &  1.0000, 0.7922, 0.7102, 0.6043, 0.4873, 0.4180, 0.3783, 0.3428, 0.3333,
     &  1.0000, 0.7875, 0.7288, 0.6249, 0.5047, 0.4296, 0.3852, 0.3446, 0.3333,
     &  1.0000, 0.8504, 0.7513, 0.6506, 0.5276, 0.4460, 0.3955, 0.3476, 0.3333,
     &  1.0000, 0.8549, 0.7794, 0.6838, 0.5594, 0.4704, 0.4114, 0.3509, 0.3333,
     &  1.0000, 0.8497, 0.8165, 0.7297, 0.6070, 0.5101, 0.4388, 0.3592, 0.3333,
     &  1.0000, 0.8444, 0.8697, 0.8011, 0.6912, 0.5899, 0.5021, 0.3810, 0.3333,
     &  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.3333/

	end module DiffTable



	real*8 function DiffCoeff(a0,g0)
	use DiffTable
	IMPLICIT NONE
	real*8 a0,g0,epsa,epsg
	
	do i=1,NAlB-1
		if(a(i).le.a0.and.a(i+1).gt.a0) then
			epsa=(a(i+1)-a0)/(a(i+1)-a(i))
			goto 1
		endif
	enddo
	i=NALB-1
	epsa=0d0
1	do j=1,NGSYM-1
		if(g(j).le.g0.and.g(j+1).gt.g0) then
			epsg=(g(j+1)-g0)/(g(j+1)-g(j))
			goto 2
		endif
	enddo
	j=NGSYM-1
	epsg=0d0
2	continue
	DiffCoeff=epsa*(epsg*D(i,j)+(1d0-epsg)*D(i,j+1))+(1d0-epsa)*(epsg*D(i+1,j)+(1d0-epsg)*D(i+1,j+1))
	return
	end
	

	subroutine DiffCoeffCell(ci,cj,iT)
	use Parameters
	IMPLICIT NONE
	integer i,j,iT,ci,cj,ii,iopac
	real*8 Ksca(nlam),Kext(nlam),dBB(nlam),spec(nlam),int1,int2
	real*8 wfunc(nlam),Kabs(nlam),g(nlam),Kqhp(nlam)

	if(iT.lt.1) iT=1
	if(iT.gt.TMAX-1) iT=TMAX-1

	if(ngrains.eq.1.and.iTD(iT).eq.iT) then
		C(ci,cj)%KDabs=KDabs(iT)
		C(ci,cj)%KDext=KDext(iT)
		C(ci,cj)%iTD=iTD(iT)
		return
	endif

	if(C(ci,cj)%opacity_set) then
		Ksca(1:nlam)=C(ci,cj)%KscaTot(1:nlam)
		Kabs(1:nlam)=C(ci,cj)%KabsTot(1:nlam)
		Kext(1:nlam)=Kabs(1:nlam)+Ksca(1:nlam)
	else
		Kext(1:nlam)=0d0
		Ksca(1:nlam)=0d0
		do ii=1,ngrains
		   do iopac=1,Grain(ii)%nopac
			Kext(1:nlam)=Kext(1:nlam)+
     &	        Grain(ii)%Kext(iopac,1:nlam)*C(ci,cj)%w(ii)*C(ci,cj)%wopac(ii,iopac)
			Ksca(1:nlam)=Ksca(1:nlam)+
     &          Grain(ii)%Ksca(iopac,1:nlam)*C(ci,cj)%w(ii)*C(ci,cj)%wopac(ii,iopac)
		   enddo
		enddo
		Kabs(1:nlam)=Kext(1:nlam)-Ksca(1:nlam)
	endif

	Kqhp(1:nlam)=0d0
	g(1:nlam)=0d0
	do ii=1,ngrains
		do iopac=1,Grain(ii)%nopac
			g(1:nlam)=g(1:nlam)+
     &          Grain(ii)%g(iopac,1:nlam)*Grain(ii)%Ksca(iopac,1:nlam)*C(ci,cj)%w(ii)*C(ci,cj)%wopac(ii,iopac)
			if(Grain(ii)%qhp) then
				Kqhp(1:nlam)=Kqhp(1:nlam)+Grain(ii)%Kabs(iopac,1:nlam)*C(ci,cj)%w(ii)*C(ci,cj)%wopac(ii,iopac)
			endif
		enddo
	enddo
	if(scattering) then
		g(1:nlam)=g(1:nlam)/Ksca(1:nlam)
	else
		g(1:nlam)=0d0
	endif

	dBB(1:nlam)=BB(1:nlam,iT+1)-BB(1:nlam,iT)

	spec(1:nlam)=dBB(1:nlam)/(Kext(1:nlam)-g(1:nlam)**2*Ksca(1:nlam))
	call integrate(spec,int1)
	call integrate(dBB,int2)
	C(ci,cj)%KDext=int2/int1

	spec(1:nlam)=dBB(1:nlam)*(Kabs(1:nlam)-Kqhp(1:nlam))
	call integrate(spec,int1)
	C(ci,cj)%KDabs=int1/int2
	if(use_qhp) then
		spec(1:nlam)=dBB(1:nlam)*Kqhp(1:nlam)
		call integrate(spec,int1)
		C(ci,cj)%KDQHP=int1/int2
	else
		C(ci,cj)%KDQHP=0d0
	endif

	C(ci,cj)%iTD=iT

	if(ngrains.eq.1) then
		KDext(iT)=C(ci,cj)%KDext
		KDabs(iT)=C(ci,cj)%KDabs
		iTD(iT)=C(ci,cj)%iTD
	endif

	return
	end





	subroutine DiffCoeffCellMin2009(ci,cj,iT)
	use Parameters
	IMPLICIT NONE
	integer i,j,iT,ci,cj,ii,iopac
	real*8 Ksca(nlam),Kext(nlam),dBB(nlam),spec(nlam),int1,int2
	real*8 wfunc(nlam),Kabs(nlam),g(nlam),Kqhp(nlam)

	if(iT.lt.1) iT=1
	if(iT.gt.TMAX-1) iT=TMAX-1

	if(ngrains.eq.1.and.iTD(iT).eq.iT) then
		C(ci,cj)%KDabs=KDabs(iT)
		C(ci,cj)%KDext=KDext(iT)
		C(ci,cj)%iTD=iTD(iT)
		return
	endif

	Kext(1:nlam)=0d0
	Ksca(1:nlam)=0d0
	Kqhp(1:nlam)=0d0
	g(1:nlam)=0d0
	do ii=1,ngrains
	   do iopac=1,Grain(ii)%nopac
		Kext(1:nlam)=Kext(1:nlam)+
     &	        Grain(ii)%Kext(iopac,1:nlam)*C(ci,cj)%w(ii)*C(ci,cj)%wopac(ii,iopac)
		Ksca(1:nlam)=Ksca(1:nlam)+
     &          Grain(ii)%Ksca(iopac,1:nlam)*C(ci,cj)%w(ii)*C(ci,cj)%wopac(ii,iopac)
		g(1:nlam)=g(1:nlam)+
     &          Grain(ii)%g(iopac,1:nlam)*Grain(ii)%Ksca(iopac,1:nlam)*C(ci,cj)%w(ii)*C(ci,cj)%wopac(ii,iopac)
		if(Grain(ii)%qhp) then
			Kqhp(1:nlam)=Kqhp(1:nlam)+Grain(ii)%Kabs(iopac,1:nlam)*C(ci,cj)%w(ii)*C(ci,cj)%wopac(ii,iopac)
		endif
	   enddo
	enddo
	Kabs(1:nlam)=Kext(1:nlam)-Ksca(1:nlam)
	if(scattering) then
		g(1:nlam)=g(1:nlam)/Ksca(1:nlam)
	else
		g(1:nlam)=0d0
	endif

	dBB(1:nlam)=BB(1:nlam,iT+1)-BB(1:nlam,iT)

	spec(1:nlam)=dBB(1:nlam)*Kabs(1:nlam)
	call integrate(spec,int1)
	dBB(1:nlam)=dBB(1:nlam)/int1
	wfunc(1:nlam)=dBB(1:nlam)*Kabs(1:nlam)

	do j=1,10
		spec(1:nlam)=Kabs(1:nlam)*wfunc(1:nlam)/Kext(1:nlam)
		call integrate(spec,int1)
		wfunc(1:nlam)=Kabs(1:nlam)*dBB(1:nlam)*int1+Ksca(1:nlam)*wfunc(1:nlam)/Kext(1:nlam)
	enddo

	spec(1:nlam)=wfunc(1:nlam)/(Kext(1:nlam)-g(1:nlam)**2*Ksca(1:nlam))
	call integrate(spec,int1)
	spec(1:nlam)=wfunc(1:nlam)/((Kext(1:nlam)-g(1:nlam)**2*Ksca(1:nlam))**2)
	call integrate(spec,int2)
	C(ci,cj)%KDext=int1/int2

	spec(1:nlam)=wfunc(1:nlam)*(Kabs(1:nlam)-Kqhp(1:nlam))/Kext(1:nlam)
	call integrate(spec,int1)
	spec(1:nlam)=wfunc(1:nlam)/Kext(1:nlam)
	call integrate(spec,int2)
	C(ci,cj)%KDabs=int1/int2
	spec(1:nlam)=wfunc(1:nlam)*Kqhp(1:nlam)/Kext(1:nlam)
	call integrate(spec,int1)
	C(ci,cj)%KDQHP=int1/int2

	C(ci,cj)%iTD=iT

	if(ngrains.eq.1) then
		KDext(iT)=C(ci,cj)%KDext
		KDabs(iT)=C(ci,cj)%KDabs
		iTD(iT)=C(ci,cj)%iTD
	endif

	return
	end



