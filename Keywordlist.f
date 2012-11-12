	subroutine Keywordlist(key)
	IMPLICIT NONE
	character*25 keyword(300)
	character*500 key
	integer i,Nnormal,Nkey

	if(key(1:1).eq.'*'.or.key.eq.' ') return
	
	i=1
	keyword(i)='startype'
	i=i+1
	keyword(i)='starfile'
	i=i+1
	keyword(i)='tstar'
	i=i+1
	keyword(i)='rstar'
	i=i+1
	keyword(i)='mstar'
	i=i+1
	keyword(i)='lstar'
	i=i+1
	keyword(i)='distance'
	i=i+1
	keyword(i)='av'
	i=i+1
	keyword(i)='posangle'
	i=i+1
	keyword(i)='incangle'
	i=i+1
	keyword(i)='tstar2'
	i=i+1
	keyword(i)='rstar2'
	i=i+1
	keyword(i)='rin'
	i=i+1
	keyword(i)='rout'
	i=i+1
	keyword(i)='nrad'
	i=i+1
	keyword(i)='ntheta'
	i=i+1
	keyword(i)='nspan'
	i=i+1
	keyword(i)='nlev'
	i=i+1
	keyword(i)='ntspan'
	i=i+1
	keyword(i)='ntlev'
	i=i+1
	keyword(i)='gridrefine'
	i=i+1
	keyword(i)='denstype'
	i=i+1
	keyword(i)='radfile'
	i=i+1
	keyword(i)='densfile'
	i=i+1
	keyword(i)='denspow'
	i=i+1
	keyword(i)='denspow2'
	i=i+1
	keyword(i)='rpow2'
	i=i+1
	keyword(i)='rpow3'
	i=i+1
	keyword(i)='rexp'
	i=i+1
	keyword(i)='shpow'
	i=i+1
	keyword(i)='sh1au'
	i=i+1
	keyword(i)='mdust'
	i=i+1
	keyword(i)='tau550'
	i=i+1
	keyword(i)='reprocess'
	i=i+1
	keyword(i)='mdot'
	i=i+1
	keyword(i)='mu0max'
	i=i+1
	keyword(i)='vexp'
	i=i+1
	keyword(i)='vexp1'
	i=i+1
	keyword(i)='vexp2'
	i=i+1
	keyword(i)='scaledisk'
	i=i+1
	keyword(i)='tmax'
	i=i+1
	keyword(i)='forcefirst'
	i=i+1
	keyword(i)='iter'
	i=i+1
	keyword(i)='raditer'
	i=i+1
	keyword(i)='epsiter'
	i=i+1
	keyword(i)='maxiter'
	i=i+1
	keyword(i)='niter0'
	i=i+1
	keyword(i)='startiter'
	i=i+1
	keyword(i)='lamgrid'
	i=i+1
	keyword(i)='lam1'
	i=i+1
	keyword(i)='lam2'
	i=i+1
	keyword(i)='nlam'
	i=i+1
	keyword(i)='zlam1'
	i=i+1
	keyword(i)='zlam2'
	i=i+1
	keyword(i)='nzlam'
	i=i+1
	keyword(i)='scattype'
	i=i+1
	keyword(i)='tcontact'
	i=i+1
	keyword(i)='tdesiter'
	i=i+1
	keyword(i)='1d'
	i=i+1
	keyword(i)='halo'
	i=i+1
	keyword(i)='fweight'
	i=i+1
	keyword(i)='storescatt'
	i=i+1
	keyword(i)='maxinteract'
	i=i+1
	keyword(i)='randomwalk'
	i=i+1
	keyword(i)='etrace'
	i=i+1
	keyword(i)='factrw'
	i=i+1
	keyword(i)='nphotdiffuse'
	i=i+1
	keyword(i)='dtdiffuse'
	i=i+1
	keyword(i)='fld'
	i=i+1
	keyword(i)='nphotfinal'
	i=i+1
	keyword(i)='nphotfirst'
	i=i+1
	keyword(i)='nfirst'
	i=i+1
	keyword(i)='idum'
	i=i+1
	keyword(i)='obstmc'
	i=i+1
	keyword(i)='tracestar'
	i=i+1
	keyword(i)='dimstar'
	i=i+1
	keyword(i)='traceemis'
	i=i+1
	keyword(i)='tracescat'
	i=i+1
	keyword(i)='tracegas'
	i=i+1
	keyword(i)='radpress'
	i=i+1
	keyword(i)='fluxcontr'
	i=i+1
	keyword(i)='topac_interpol'
	i=i+1
	keyword(i)='ngrains'
	i=i+1
	keyword(i)='compositionfile'
	i=i+1
	keyword(i)='tdesqhp'
	i=i+1
	keyword(i)='scalesh'
	i=i+1
	keyword(i)='shscale'
	i=i+1
	keyword(i)='shpow'
	i=i+1
	keyword(i)='shscalefile'
	i=i+1
	keyword(i)='psettr0'
	i=i+1
	keyword(i)='psettpow'
	i=i+1
	keyword(i)='psettphi0'
	i=i+1
	keyword(i)='rimscale'
	i=i+1
	keyword(i)='rimwidth'
	i=i+1
	keyword(i)='nruns'
	i=i+1
	keyword(i)='viscous'
	i=i+1
	keyword(i)='fastviscous'
	i=i+1
	keyword(i)='alphaviscous'
	i=i+1
	keyword(i)='getalpha'
	i=i+1
	keyword(i)='powviscous'
	i=i+1
	keyword(i)='tgas'
	i=i+1
	keyword(i)='g2d_heat'
	i=i+1
	keyword(i)='convection'
	i=i+1
	keyword(i)='forcediff'
	i=i+1
	keyword(i)='multiwav'
	i=i+1
	keyword(i)='nbw'
	i=i+1
	keyword(i)='nspike'
	i=i+1
	keyword(i)='scset'
	i=i+1
	keyword(i)='mpset'
	i=i+1
	keyword(i)='mpstr'
	i=i+1
	keyword(i)='scseteq'
	i=i+1
	keyword(i)='scsetsave'
	i=i+1
	keyword(i)='lifetime'
	i=i+1
	keyword(i)='alphaturb'
	i=i+1
	keyword(i)='gas2dust'
	i=i+1
	keyword(i)='qturb'
	i=i+1
	keyword(i)='thinparticle'
	i=i+1
	keyword(i)='mrn'
	i=i+1
	keyword(i)='mrn_index'
	i=i+1
	keyword(i)='mrn_rmin'
	i=i+1
	keyword(i)='mrn_rmax'
	i=i+1
	keyword(i)='mrn_ngrains'
	i=i+1
	keyword(i)='gsd'
	i=i+1
	keyword(i)='gsd_full'
	i=i+1
	keyword(i)='gsd_plot'
	i=i+1
	keyword(i)='gsd_rmin'
	i=i+1
	keyword(i)='gsd_rmax'
	i=i+1
	keyword(i)='gsd_xi'
	i=i+1
	keyword(i)='gsd_vfrag'
	i=i+1
	keyword(i)='gsd_diag'
	i=i+1
	keyword(i)='deadzone'
	i=i+1
	keyword(i)='deadheight'
	i=i+1
	keyword(i)='deadalpha'
	i=i+1
	keyword(i)='hotgasminrad'
	i=i+1
	keyword(i)='hotgasmaxrad'
	i=i+1
	keyword(i)='hotgast'
	i=i+1
	keyword(i)='meixa'
	i=i+1
	keyword(i)='meixb'
	i=i+1
	keyword(i)='meixc'
	i=i+1
	keyword(i)='meixd'
	i=i+1
	keyword(i)='meixe'
	i=i+1
	keyword(i)='meixf'
	i=i+1
	keyword(i)='meixg'
	i=i+1
	keyword(i)='meixrsw'
	i=i+1
	keyword(i)='meixrin'
	i=i+1
	keyword(i)='timeshift'
	i=i+1
	keyword(i)='irf'
	i=i+1
	keyword(i)='t_irf'
	i=i+1
	keyword(i)='f_irf'
	i=i+1
	keyword(i)='exportprodimo'
	i=i+1
	keyword(i)='tsmooth'
	i=i+1
	keyword(i)='outputfits'
	
	Nnormal=i
c from here on keywords with numbers
	i=i+1
	keyword(i)='tdesa'
	i=i+1
	keyword(i)='tdesb'
	i=i+1
	keyword(i)='forcegf'
	i=i+1
	keyword(i)='tdesfast'
	i=i+1
	keyword(i)='trace'
	i=i+1
	keyword(i)='abun'
	i=i+1
	keyword(i)='mat'
	i=i+1
	keyword(i)='part'
	i=i+1
	keyword(i)='fitspart'
	i=i+1
	keyword(i)='opac'
	i=i+1
	keyword(i)='aggmix'
	i=i+1
	keyword(i)='topac'
	i=i+1
	keyword(i)='powmix'
	i=i+1
	keyword(i)='radmix'
	i=i+1
	keyword(i)='tmix'
	i=i+1
	keyword(i)='qhp'
	i=i+1
	keyword(i)='minrad'
	i=i+1
	keyword(i)='maxrad'
	i=i+1
	keyword(i)='shaperad'
	i=i+1
	keyword(i)='roundtype'
	i=i+1
	keyword(i)='roundwidth'
	i=i+1
	keyword(i)='roundpow'
	i=i+1
	keyword(i)='roundpeak'
	i=i+1
	keyword(i)='settle'
	i=i+1
	keyword(i)='settlefile'
	i=i+1
	keyword(i)='shtype'
	i=i+1
	keyword(i)='asym'
	i=i+1
	keyword(i)='asymb'
	i=i+1
	keyword(i)='wasymb'
	i=i+1
	keyword(i)='pmax'
	i=i+1
	keyword(i)='gradslope'
	i=i+1
	keyword(i)='gradinner'
	i=i+1
	keyword(i)='gradclose'
	i=i+1
	keyword(i)='gradfar'
	i=i+1
	keyword(i)='gradrad'
	i=i+1
	keyword(i)='startgap'
	i=i+1
	keyword(i)='endgap'
	i=i+1
	keyword(i)='gap'
	i=i+1
	keyword(i)='shapegap'
	i=i+1
	keyword(i)='roundtypegap'
	i=i+1
	keyword(i)='roundpowgap'
	i=i+1
	keyword(i)='rfix'
	i=i+1
	keyword(i)='coupledabun'
	i=i+1
	keyword(i)='coupledfrac'
	i=i+1
	keyword(i)='frac'
	i=i+1
	keyword(i)='tau1_lam'
	i=i+1
	keyword(i)='rgrain'
	i=i+1
	keyword(i)='rhograin'
	i=i+1
	keyword(i)='zone'

	Nkey=i
	
	do i=1,Nkey
		if(trim(keyword(i)).eq.key(1:len_trim(keyword(i)))) return
	enddo
	write(*,'("Keyword not recognized: ",a)') trim(key)
	write(9,'("Keyword not recognized: ",a)') trim(key)

	open(unit=81,file='keywordlist.txt',RECL=1000)
	do i=1,Nkey
		write(81,*) trim(keyword(i))
	enddo
	close(unit=81)

	write(*,'("Keyword list written to file")')
	write(9,'("Keyword list written to file")')

	stop
	
	return
	end
