MODULE FIT_MODULE

doubleprecision :: mu,m_p,k_b,pi            ! some mutual constants (see below)

doubleprecision :: poly_coeff(1:6)          ! used for polynomial solving
doubleprecision :: polyn

parameter(pi       = 3.141593)              ! PI
parameter(k_b      = 1.380658d-16)          ! Boltzmann constant in erg/K
parameter(mu       = 2.3d0)                 ! mean molecular mass in proton masses
parameter(m_p      = 1.6726231d-24)         ! proton mass in g
parameter(sig_h2   = 2d-15)                 ! cross section of H2
parameter(Grav     = 6.67259d-8)            ! gravitational constant in cm^3 g^-1 s^-2

PRIVATE

PUBLIC :: fit_function                      ! fit_function is the public function
PUBLIC :: poly_coeff                        ! poly_coeff need to be shared for solving the polynomial

CONTAINS
! __________________________________________________________________________________________________________________________________
! This subroutine is a wrapper function which calls the proper fitting
! function according to the input xi
!
! USAGE: fit_function(fit,a_01,a_12,a_l,a_p,a_r,a_sett,nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
!
! INPUT:    nm      = number of mass grid points   []
!           xi      = fragmentation power law      []
!           T       = mid-plane temperature        [K]
!           alpha   = turbulence alpha             []
!           sigma_g = gas surface density          [g cm^-2]
!           sigma_d = dust surface density         [g cm^-2]
!           rho_s   = solid density of the grains  [g cm^-3]
!           m_grid  = mass grid                    [g]
!           a_grid  = size grid                    [cm]
!           m_star  = stellar mass                 [g]
!           R       = mid-plane distance to star   [cm]
!           v_frag  = fragmentation velocity       [cm s^-1]
!
! OUTPUT:   fit     = fit-distribution in g cm^-2 in each mass bin
!           a_01    = first regime boundary
!           a_12    = second regime boundary
!           a_l     = left wing of peak
!           a_p     = center of peak
!           a_r     = right wing of peak
!           a_sett  = where settling starts to play a role
! __________________________________________________________________________________________________________________________________
subroutine fit_function(fit,a_01,a_12,a_l,a_p,a_r,a_sett,nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
    implicit none
    integer,         intent(in)  :: nm
    doubleprecision, intent(in)  :: xi,T,alpha,sigma_g,sigma_d,rho_s,m_star,R,v_frag
    doubleprecision, intent(in)  :: m_grid(1:nm),a_grid(1:nm)
    doubleprecision, intent(out) :: fit(1:nm),a_01,a_12,a_l,a_p,a_r,a_sett
    integer                      :: i_01,i_12,i_l,i_p,i_r,i_sett,i_inc,i,ir(1:4)
    doubleprecision              :: nu,L,N,sig
    doubleprecision              :: alpha_0,alpha_1,alpha_2,slope_0,slope_1,slope_2,dv,jumpfactor,inc_factor
    doubleprecision              :: alpha_0_s,alpha_1_s,alpha_2_s,slope_0_s,slope_1_s,slope_2_s
    doubleprecision              :: slopes_s(1:3),slopes(1:3)
    doubleprecision              :: bump(1:nm)
    !
    ! some general things
    !
    if (a_grid(1)>=1d-5) then
        write(*,*) 'NOTE:    recommended a_min is 0.025 micron.'
    endif
    !
    ! call the function
    !
    if (abs(xi-1.8d0)<0.05d0) then
        call fit_function18(fit,a_01,a_12,a_l,a_p,a_r,a_sett,nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
    else
        write(*,*) 'ERROR:   xi-values much different that 1.8 are not implemented yet!'
        stop 34242
    endif
end subroutine fit_function
! ==================================================================================================================================


! __________________________________________________________________________________________________________________________________
! This subroutine provides a fit function to the grain size distribution of the
! given parameters BUT ONLY FOR A XI OF (about) 1.8!!
!
! USAGE: fit_function18(fit,a_01,a_12,a_l,a_p,a_r,a_sett,nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
!
! INPUT:    nm      = number of mass grid points   []
!           xi      = fragmentation power law      []
!           T       = mid-plane temperature        [K]
!           alpha   = turbulence alpha             []
!           sigma_g = gas surface density          [g cm^-2]
!           sigma_d = dust surface density         [g cm^-2]
!           rho_s   = solid density of the grains  [g cm^-3]
!           m_grid  = mass grid                    [g]
!           a_grid  = size grid                    [cm]
!           m_star  = stellar mass                 [g]
!           R       = mid-plane distance to star   [cm]
!           v_frag  = fragmentation velocity       [cm s^-1]
!
! OUTPUT:   fit     = fit-distribution in g cm^-2 in each mass bin
!           a_01    = first regime boundary
!           a_12    = second regime boundary
!           a_l     = left wing of peak
!           a_p     = center of peak
!           a_r     = right wing of peak
!           a_sett  = where settling starts to play a role
! __________________________________________________________________________________________________________________________________
subroutine fit_function18(fit,a_01,a_12,a_l,a_p,a_r,a_sett,nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
    implicit none
    integer,         intent(in)  :: nm
    doubleprecision, intent(in)  :: xi,T,alpha,sigma_g,sigma_d,rho_s,m_star,R,v_frag
    doubleprecision, intent(in)  :: m_grid(1:nm),a_grid(1:nm)
    doubleprecision, intent(out) :: fit(1:nm),a_01,a_12,a_l,a_p,a_r,a_sett
    integer                      :: i_01,i_12,i_l,i_p,i_r,i_sett,i_inc,i,ir(1:4)
    doubleprecision              :: nu,L,N,sig
    doubleprecision              :: alpha_0,alpha_1,alpha_2,slope_0,slope_1,slope_2,dv,jumpfactor,inc_factor
    doubleprecision              :: alpha_0_s,alpha_1_s,alpha_2_s,slope_0_s,slope_1_s,slope_2_s
    doubleprecision              :: slopes_s(1:3),slopes(1:3)
    doubleprecision              :: bump(1:nm)
    !
    ! get the sizes of the boundaries and the according indices
    !
    call get_boundaries(nm,T,alpha,sigma_g,rho_s,m_grid,a_grid,m_star,R,v_frag,a_01,a_12,a_l,a_p,a_r,a_sett)
    i_01   = first(nm,a_grid>=a_01)
    i_12   = first(nm,a_grid>=a_12)
    i_l    = first(nm,a_grid>=a_l)
    i_p    = first(nm,a_grid>=a_p)
    i_r    = first(nm,a_grid>=a_r)
    i_sett = first(nm,a_grid>=a_sett)
    !
    ! some consistency checks
    !
    if ((i_01==0).or.(i_12==0).or.(i_l==0).or.(i_p==0).or.(i_r==0)) then
        write(*,*) 'ERROR:   indices are not correct!'
        write(*,*) 'i_01   = ',i_01
        write(*,*) 'i_12   = ',i_12
        write(*,*) 'i_l    = ',i_l
        write(*,*) 'i_p    = ',i_p
        write(*,*) 'i_r    = ',i_r
        write(*,*) 'i_sett = ',i_sett
        stop 18072
    endif
    if (log10(maxval((/a_l,a_p,a_r/)))-log10(a_12)<0.1d0) then
        write(*,*) 'WARNING: maximum grain size is very small,'
        write(*,*) '         fit results are not reliable!'
    endif
    !
    ! define the slopes without settling
    !
    nu = 1d0/6d0
    alpha_0 = (nu+xi+1d0)/2d0
    slope_0 = (6d0-3d0*alpha_0)

    nu = 1d0
    alpha_1 = (nu+xi+1d0)/2d0
    slope_1 = (6d0-3d0*alpha_1)

    nu = 5d0/6d0
    alpha_2 = (nu+xi+1d0)/2d0
    slope_2 = (6d0-3d0*alpha_2)
    !
    ! define the slopes with settling
    !
    nu        = 2d0/6d0
    alpha_0_s = (nu+xi+1d0)/2d0
    slope_0_s = 6d0-3d0*alpha_0_s

    nu        = 7d0/6d0
    alpha_1_s = (nu+xi+1d0)/2d0
    slope_1_s = 6d0-3d0*alpha_1_s

    nu        = 1d0
    alpha_2_s = (nu+xi+1d0)/2d0
    slope_2_s = 6d0-3d0*alpha_2_s
    !
    ! plug them together to an array each
    !
    slopes   = (/ slope_0,  slope_1,  slope_2   /)
    slopes_s = (/ slope_0_s,slope_1_s,slope_2_s /)
    !
    ! construct the fit function
    !
    ! EXPERIMENTAL: the jump factor is something between 1 and 2.5 it
    ! seems.
    !
    dv = sqrt(3d0/2d0*alpha)*sqrt(3d0)*sqrt(k_b*T/mu/m_p)*sqrt(1d0/(sigma_g*2d0*1.6d0))*(8d0*mu*m_p*sigma_g/alpha/sig_h2)**0.25d0
    !
    ! make the powerlaws with settling
    !
    ir = (/1,i_01,i_12,i_p/)
    fit    = 0d0
    fit(1) = 1
    do i = 1,3
        if (i==3) then
            ! old one:
            !jumpfactor = max(1,1.2+log(dv/30));
            ! this one works good for V_FRAG = 1 m/s
            !jumpfactor = min(2.5,max(1.1,1+dv/40));
            ! lets try to scale it with V_FRAG
            jumpfactor = min(2.5d0,max(1.1d0,1d0 + dv/50d0 * 1d0/(v_frag/100d0)))
        else
            jumpfactor = 1
        endif
        if (i_sett<=ir(i)) then
            fit(ir(i):ir(i+1))  = fit(ir(i))/jumpfactor*(a_grid(ir(i):ir(i+1)) /a_grid(ir(i))) **slopes_s(i)
        elseif ((ir(i)<=i_sett).and.(i_sett<=ir(i+1))) then
            fit(ir(i):i_sett)   = fit(ir(i))/jumpfactor*(a_grid(ir(i):i_sett)  /a_grid(ir(i))) **slopes(i)
            fit(i_sett:ir(i+1)) = fit(i_sett)          *(a_grid(i_sett:ir(i+1))/a_grid(i_sett))**slopes_s(i)
        elseif (i_sett>ir(i+1)) then
            fit(ir(i):ir(i+1))  = fit(ir(i))/jumpfactor*(a_grid(ir(i):ir(i+1)) /a_grid(ir(i))) **slopes(i)
        else
            write(*,*) 'ERROR:   something went wrong with the regimes!'
            stop 19723
        endif
    enddo
    !
    ! cut off values larger than i_right
    !
    if (i_r>i_p) then
        fit( (max(i_p,max(i_r,i_l))+1) :nm) = 0d0
    endif
    !
    ! smooth the fit
    !
    call smooth1D(nm,3,fit)
    !
    ! EXPERIMENTAL: SLOPE-INCREASE AT UPPER END OF DISTRIBUTION
    ! we will do this if a_12 and a_p are too close to each other
    !
    !
    if (log10(a_p/a_12)<0.7) then
        i_inc = last(nm,a_grid<0.2d0*a_grid(i_p))
        inc_factor = 1d0
        if (i_inc==0) then
            write(*,*) 'WARNING: i_inc == 0 => i_p seems to be quite small, consider a smaller m_min'
            i_inc = 1
        endif
        fit(i_inc:i_p) = fit(i_inc:i_p)+inc_factor*fit(i_inc:i_p)*(1-(a_grid(i_inc:i_p)-a_grid(i_p))/(a_grid(i_inc)-a_grid(i_p)))
    endif
    !
    ! produce a bump
    !
    ! GAUSSIAN:
    L    = fit(i_l)
    N    = 2d0*fit(i_l)
    sig  = max( abs(a_grid(i_p-1)-a_grid(i_p)) , min(abs(a_grid(i_r)-a_grid(i_p)),abs(a_grid(i_l)-a_grid(i_p))))/sqrt(2d0*log(N/L))
    bump = N*exp(-(a_grid-a_grid(i_p))**2d0/(2d0*sig**2d0))
    !
    ! add the bump
    !
    do i = i_12,nm
       fit(i) = max(fit(i),bump(i))
    enddo
    !
    ! normalize the thing
    !
    fit = sigma_d*fit/sum(fit)
    if (any(isnan(fit))) write(*,*) 'NaN occured in renormalization'
end subroutine fit_function18
! ==================================================================================================================================

! __________________________________________________________________________________________________________________________________
! A new function that gives the velocities according to Ormel and Cuzzi (2007)
!
! INPUT:   tau_1       =   stopping time 1
!          tau_2       =   stopping time 2
!          t0          =   large eddy turnover time
!          v0          =   large eddy velocity
!          ts          =   small eddy turnover time
!          vs          =   small eddy velocity
!          Reynolds    =   Reynolds number
!
! RETURNS: v_rel_ormel =   relative velocity SQUARED
!
! __________________________________________________________________________________________________________________________________
doubleprecision function v_rel_ormel(tau_1, tau_2,t0,v0,ts,vs,Reynolds)
    implicit none
    doubleprecision, intent(in)            :: tau_1, tau_2
    doubleprecision, intent(in)            :: t0,v0,ts,vs,Reynolds
    doubleprecision                        :: St1, St2, tau_mx, tau_mn, Vg2
    doubleprecision                        :: c0,c1,c2,c3, y_star, ya, eps
    doubleprecision                        :: hulp1, hulp2, coll_term
    !
    ! sort tau's 1--> correspond to the max. now
    ! (a bit confusing perhaps)
    !
    if (tau_1 .GE. tau_2) then
        tau_mx = tau_1
        tau_mn = tau_2
        St1 = tau_mx/t0
        St2 = tau_mn/t0
    else
        tau_mx = tau_2
        tau_mn = tau_1
        St1 = tau_mx/t0
        St2 = tau_mn/t0
    endif
    !
    ! note the square
    !
    Vg2 = 1.5 *v0**2.0
    !
    ! approximate solution for St*=y*St1; valid for St1 << 1.
    !
    ya = 1.6
    if (tau_mx .LT. 0.2*ts) then
       !
       ! very small regime
       !
       v_rel_ormel = 1.5 *(vs/ts *(tau_mx - tau_mn))**2.0
    elseif (tau_mx .LT. ts/ya) then
       v_rel_ormel = Vg2 *(St1-St2)/(St1+St2)*(St1**2.0/(St1+Reynolds**(-0.5)) - St2**2.0/(St2+Reynolds**(-0.5)))
    elseif (tau_mx .LT. 5.0*ts) then
        !
        ! Eq. 17 of OC07. The second term with St_i**2.0 is negligible (assuming !Re>>1)
        ! hulp1 = Eq. 17
        ! hulp2 = Eq. 18
        !
        hulp1 = ( (St1-St2)/(St1+St2) * (St1**2.0/(St1+ya*St1) - St2**2.0/(St2+ya*St1)) )!note the -sign
        hulp2 = 2.0*(ya*St1-Reynolds**(-0.5)) + St1**2.0/(ya*St1+St1) - St1**2.0/(St1+Reynolds**(-0.5)) +&
        St2**2.0/(ya*St1+St2) - St2**2.0/(St2+Reynolds**(-0.5))
        v_rel_ormel = Vg2 *(hulp1 + hulp2)
    elseif (tau_mx .LT. t0/5.0) then
        !
        ! stopping time ratio
        !
        eps=St2/St1
        !
        ! Full intermediate regime
        !
        v_rel_ormel = Vg2 *( St1*(2.0*ya - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+ya) + eps**3.0/(ya+eps) )) )
    elseif (tau_mx .LT. t0) then
        !
        ! now y* lies between 1.6 (St1 << 1) and 1.0 (St1>=1).
        ! The fit below fits ystar to less than 1%
        !
        c3 =-0.29847604
        c2 = 0.32938936
        c1 =-0.63119577
        c0 = 1.6015125
        y_star = c0 + c1*St1 + c2*St1**2.0 + c3*St1**3.0
        !
        !stopping time ratio
        !
        eps=St2/St1
        !
        ! we can then employ the same formula as before
        !
        v_rel_ormel = Vg2 *( St1*(2.0*y_star - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+y_star) + eps**3.0/(y_star+eps) )) )
    else
        !
        ! heavy particle limit
        !
        v_rel_ormel = Vg2 *( 1.0/(1.0+St1) + 1.0/(1.0+St2) )
    endif
end function v_rel_ormel
! ==================================================================================================================================

! __________________________________________________________________________________________________________________________________
! This function calculates the boundaries of the different regimes in the
! grain size distribution.
!
! INPUT:
!   nm          = # of size/mass points
!   T           = temperature               [K]
!   alpha       = turbulence parameter      []
!   sigma_g     = gas surface density       [g/cm^2]
!   rho_s       = dust grain volume density [g/cm^3]
!   m_grid      = mass array                [g]
!   a_grid      = grain size array          [cm]
!   m_star      = stellar mass              [g]
!   R           = radial distance to star   [cm]
!   v_frag      = fragmentation velocity    [cm/s]
!
! OUTPUT:
! position (in grain size) of the transisions:
!   a_01    = first transition between brownian motion and turbulence
!   a_12    = linear to intermediate regime of turbulende
!   a_left  = left side of cratering peak
!   a_peak  = position of cratering peak
!   a_right = right side of cratering peak
! __________________________________________________________________________________________________________________________________
subroutine get_boundaries(nm,T,alpha,sigma_g,rho_s,m_grid,a_grid,m_star,R,v_frag,a_01,a_12,a_left,a_peak,a_right,a_sett)
    implicit none
    integer,         intent(in)  :: nm
    doubleprecision, intent(in)  :: T,alpha,sigma_g,rho_s,m_star,R,v_frag
    doubleprecision, intent(in)  :: m_grid(1:nm),a_grid(1:nm)
    doubleprecision, intent(out) :: a_01,a_12,a_left,a_peak,a_right,a_sett
    !
    ! used for polynomial solving
    !
    doubleprecision,external  :: polyn
    integer                   :: ierr
    !
    ! other
    !
    doubleprecision :: cs,Re,a1,A,ya,omega,ro,tn,ts,vn,vs,tau_1,tau_2,xL,xR,yL,yR
    doubleprecision :: dv_BM(1:nm,1:nm),dv_TM(1:nm,1:nm),dv_ii(1:nm),dv_i1(1:nm)
    integer         :: i,j,i_peak,i_left,i_right
    !
    ! initialize output
    !
    a_01    = 0d0
    a_12    = 0d0
    a_left  = 0d0
    a_peak  = 0d0
    a_right = 0d0
    a_sett  = 0d0
    !
    ! sound speed, reynolds number and grainsize(1)
    !
    cs     = sqrt(k_b*T/mu/m_p)
    Re     = alpha*sig_h2*sigma_g/(2d0*mu*m_p)
    a1     = a_grid(1)
    !
    ! get the size where particles start to settle
    !
    a_sett = 2d0*alpha*sigma_g/(pi*rho_s)
    !
    ! get the transition from BM to turbulence
    !
    A          = 0.1d0 * 32d0*k_b*T*sigma_g**2d0/(pi**4d0*rho_s**3d0*sqrt(Re)*alpha*cs**2d0)
    poly_coeff = (/1d0,-2d0*a1,a1**2d0,0d0,0d0,-A/)
    !
    ! solve polynomial
    !
    a_01 = zbrent(polyn,a_grid(1),a_grid(nm),1d-10,ierr)
    !
    ! if ZBRENT fails, then we give up
    !
    if (ierr/=0) then
        write(*, *) 'ERROR:   Failure by ZBRENT, exiting'
        stop 92378
    endif
    !
    ! get the turbulence bump
    !
    ya     = 1.6d0
    a_12   = 1d0/(ya*pi*rho_s)*sqrt(8d0*mu*m_p*sigma_g/(alpha*sig_h2))
    !
    ! get the cratering bump positions
    !
    !
    ! calculate relative velocities due to brownian motion
    !
    do i=1,nm
        do j=1,nm
            dv_BM(i,j) = sqrt(8d0*k_b*T*(m_grid(i)+m_grid(j))/pi/m_grid(i)/m_grid(j))
        enddo
    enddo
    !
    ! calculate relative velocities due to turbulence
    !
    do i=1,nm
        do j=1,nm
            omega = sqrt(Grav*m_star/R**3d0)
            ro    = sigma_g*omega/(sqrt(2*pi)*cs)
            tn    = 1d0/omega
            ts    = tn*Re**(-0.5d0)
            vn    = sqrt(alpha)*cs
            vs    = vn*Re**(-0.25d0)

            tau_1 = rho_s*a_grid(i)/sigma_g/omega*pi/2d0
            tau_2 = rho_s*a_grid(j)/sigma_g/omega*pi/2d0

            dv_TM(i,j) = v_rel_ormel(tau_1,tau_2,tn,vn,ts,vs,Re)
        enddo
    enddo
    !
    ! put the velocities together,
    ! note that the turbulent ones are already squared
    !
    do i = 1,nm
        dv_ii(i) = sqrt(dv_BM(i,i)**2d0 + dv_TM(i,i))
        dv_i1(i) = sqrt(dv_BM(i,1)**2d0 + dv_TM(i,1))
    enddo
    !
    ! the position of the peak
    !
    if (maxval(dv_ii)<V_FRAG) then
        write(*,*) 'ERROR:   particles will not fragment!'
        write(*,*) '         Either increase turbulent velocities or decrease'
        write(*,*) '         the fragmentation velocity'
        stop 38332
    endif
    i_peak = first(nm,(dv_ii-0.8d0*v_frag >=0d0))
    if (i_peak<2) then
        i_peak = first(nm,(dv_ii-0.8d0*v_frag>=0d0).and.(a_grid>1d-4))
        if (i_peak<2) then
            write (*,*) 'ERROR:   cannot find i_peak'
            stop 12783
        else
            a_peak = a_grid(i_peak)
        endif
    else
        xL = a_grid(i_peak-1)
        xR = a_grid(i_peak)
        yL = dv_ii(i_peak-1)-0.8d0*v_frag
        yR = dv_ii(i_peak)  -0.8d0*v_frag
        a_peak = xL-yL*(xR-xL)/(yR-yL)
    endif
    !
    ! the position of the left wing
    !
    i_left = first(nm,dv_i1-0.8d0*v_frag>=0)
    if (i_left<2) then
        i_left    = first(nm,(dv_i1-0.8d0*v_frag>=0d0).and.(a_grid>1d-4))
        if (i_left<2) then
            write(*,*) 'ERROR:   cannot find i_left'
            stop 23212
        else
            a_left = a_grid(i_left)
        endif
    else
        xL = a_grid(i_left-1)
        xR = a_grid(i_left)
        yL = dv_i1(i_left-1)-0.8d0*v_frag
        yR = dv_i1(i_left)  -0.8d0*v_frag
        a_left = xL-yL*(xR-xL)/(yR-yL)
    endif
    !
    ! the position of the right wing
    !
    i_right = first(nm,dv_ii-v_frag>=0d0)
    if (i_right<2) then
        i_right = first(nm,(dv_ii-v_frag>=0d0).and.(a_grid>1d-4))
        if (i_right<2) then
            write(*,*) 'ERROR:   cannot find i_right'
            stop 27431
        else
            a_right = a_grid(i_right)
        endif
    else
        xL = a_grid(i_right-1)
        xR = a_grid(i_right)
        yL = dv_ii(i_right-1)-v_frag
        yR = dv_ii(i_right)  -v_frag
        a_right = xL-yL*(xR-xL)/(yR-yL)
    endif

end subroutine get_boundaries
! ==================================================================================================================================

! __________________________________________________________________________________________________________________________________
!                   FUNCTION: FIND LAST OF MASK
!
! finds the last element where mask is true
! example:
!
! A = (/0,1,2,3,4,5,6,7/)
!
! last(size(A),A<2) = 2
! __________________________________________________________________________________________________________________________________
integer function last(n,mask)
implicit none
integer,intent(in) :: n
logical,intent(in) :: mask(1:n)
integer :: i

last = 0
do i=n,1,-1
    if (mask(i)) then
        last=i
        return
    endif
enddo

end function last
! ==================================================================================================================================

! __________________________________________________________________________________________________________________________________
!                   FUNCTION: FIND FIRST OF MASK
!
! finds the first element where mask is true
! example:
!
! A = (/0,1,2,3,4,5,6,7/)
!
! first(size(A),A>2) = 4
! __________________________________________________________________________________________________________________________________
integer function first(n,mask)
implicit none
integer,intent(in) :: n
logical,intent(in) :: mask(1:n)
integer :: i

first = 0
do i=1,n
    if (mask(i)) then
        first=i
        return
    endif
enddo

end function first
! ==================================================================================================================================



! __________________________________________________________________________________________________________________________________
! This function smoothes the given 1D array of size n by a running average over
! 3 points. The averaging is performed n_smooth times.
!
! USAGE:   smooth1D(n,n_smooth,arr_in)
!
! INPUT:   n        = size of in/output array (1D)
!          n_smooth = how much the array is to be smoothed
!          arr_in   = array which is to be smoothed
!
! OUTPUT:  the input array is smoothed
! __________________________________________________________________________________________________________________________________
subroutine smooth1D(n,n_smooth,arr_in)
implicit none
integer,         intent(in)    :: n,n_smooth
doubleprecision, intent(inout) :: arr_in(1:n)
doubleprecision                :: arr_dummy(1:n)
integer                        :: i_smooth,i
!
! some checks
!
if (n_smooth<1) then
    write(*,*) 'ERROR:   n_smooth has to be larger than 0!'
    stop 72311
endif
if (n_smooth>=0.5d0*n) then
    write(*,*) 'WARNING: n_smooth is too large, this way you will average out'
    write(*,*) '         any trends in your array!'
endif
!
! average n_smooth times
!
arr_dummy = arr_in
do i_smooth = 1,n_smooth
    do i = 1,n
        arr_in(i) = ( arr_dummy(max(1,i-1))+arr_dummy(i)+arr_dummy(min(i+1,n)) )/3d0
    enddo
    arr_dummy = arr_in
enddo

end subroutine smooth1D
! ==================================================================================================================================

! __________________________________________________________________________________________________________________________________
! Brent's Algorithm for root finding
!
! USAGE:     zbrent(func,x1,x2,tol,ierr)
!
! INPUT:     func   =    the double precision function to be worked with
!            x1,x2  =    the interval around the root
!            tol    =    the tolerance around the root
!            ierr   =    0 means no error, 1 means error
!
! OUTPUT:    zbrent =    the root of function "func" with a precision of "tol"
!
! __________________________________________________________________________________________________________________________________
doubleprecision function zbrent(func,x1,x2,tol,ierr)
    implicit none
    external func
    doubleprecision             :: func
    doubleprecision, intent(in) :: tol,x1,x2
    integer, intent(out)        :: ierr
    integer                     :: iter
    doubleprecision             :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    integer,parameter         :: ITMAX=100
    doubleprecision,parameter :: EPSS=3.d-8
    ierr = 0
    a    = x1
    b    = x2
    fa   = func(a)
    fb   = func(b)

    if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
        write(*,*) 'root must be bracketed for zbrent'
        ierr = 1
    endif
    c  = b
    fc = fb
    do iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
            c  = a
            fc = fa
            d  = b-a
            e  = d
        endif
        if(abs(fc).lt.abs(fb)) then
            a  = b
            b  = c
            c  = a
            fa = fb
            fb = fc
            fc = fa
        endif
        tol1 = 2d0*EPSS*abs(b)+0.5d0*tol
        xm   = 0.5d0*(c-b)
        if(abs(xm)<=tol1 .or. fb==0.)then
            zbrent = b
            return
        endif
        if (abs(e)>=tol1 .and. abs(fa)>abs(fb)) then
            s = fb/fa
            if(a==c) then
                p = 2.d0*xm*s
                q = 1.d0-s
            else
                q = fa/fc
                r = fb/fc
                p = s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
                q = (q-1.)*(r-1.)*(s-1.)
            endif
            if(p>0.) q = -q
            p = abs(p)
            if(2.d0 * p < min(3.*xm*q-abs(tol1*q),abs(e*q))) then
                e = d
                d = p/q
            else
                d = xm
                e = d
            endif
        else
            d = xm
            e = d
        endif
        a  = b
        fa = fb
        if (abs(d) > tol1) then
            b = b+d
        else
            b = b + sign(tol1,xm)
        endif
        fb = func(b)
    enddo
    write(*,*) 'WARNING: zbrent exceeding maximum iterations'
    zbrent=b
    ierr = 2
    return
end function zbrent
! ==================================================================================================================================

END MODULE FIT_MODULE
! ==================================================================================================================================

! __________________________________________________________________________________________________________________________________
! polyn(x) is a function which calculates the 5th order polynomial with the coefficients given by the array poly_coeff which is
! supplied (and set within) fit_module
!
! USAGE:  polyn(x)
!
! INPUT:  x = value at which the polynomial is evaluated
!
! OUTPUT: returns the value of the 5th order poynomial evaluated at x. Coefficients are supplied by fit_module
!
! __________________________________________________________________________________________________________________________________
doubleprecision function polyn(x)
    use fit_module, ONLY: poly_coeff
    implicit none
    doubleprecision, intent(in) :: x
    polyn = poly_coeff(1)*x**5+poly_coeff(2)*x**4+poly_coeff(3)*x**3+poly_coeff(4)*x**2+poly_coeff(5)*x+poly_coeff(6)
end function polyn
! ==================================================================================================================================
