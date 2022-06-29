!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Kazem Rezazadeh                                   June 19, 2022
! kazem.rezazadeh@ipm.ir
! https://github.com/krezazadeh
!
! This program calculates the present-day fractional energy density
! of the induced gravitational waves (Omega_GW0) in the setup of
! alpha-attractor inflation with a tiny-bump. This program takes
! the file "Ps.txt" as input. The file "Ps.txt" contains the results
! of the model for the scalar power spectrum and it can be generated
! by either "Ps_alpha_attractor.f90" or "fPBH_alpha_attractor.f90".
! The output results of the code are stored in the file "OmegaGW.txt".
!
! To cite this code, please cite the following paper:
! https://arxiv.org/abs/2110.01482
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program OmegaGW0_alpha_attractor

    implicit none

    ! real(8), parameter :: pi = acos(-1.0d0)

    real(8), parameter :: MP = 1.0d0
    real(8), parameter :: Mpc = 3.8082082148306456d56/MP

    integer, parameter :: nPs = 1000

    real(8), dimension(nPs) :: array_k, array_Ps
    common /arrays_Ps/ array_k, array_Ps

    integer, parameter :: n_array_Si = 100000

    real(8), dimension(n_array_Si) :: array_x_Si, array_Si
    common /arrays_Si/ array_x_Si, array_Si

    integer, parameter :: n_array_Ci = 100000

    real(8), dimension(n_array_Ci) :: array_x_Ci, array_Ci
    common /arrays_Ci/ array_x_Ci, array_Ci

    integer :: i, j, l
    integer, parameter :: ni = 1000
    integer, parameter :: nj = 500
    integer, parameter :: nl = 100
    real(8) :: kmin, kmax
    real(8) :: ftoHzinitial, ftoHzfinal
    real(8) :: result
    real(8) :: ftoHzi, ki
    real(8) :: vmini, vmaxi
    real(8) :: vj, vjp1, dvj
    real(8) :: uminj, umaxj
    real(8) :: ul, ulp1, dul
    real(8) :: OmegaGW0i
    real(8) :: kftoHz, ftoHzk
    real(8) :: fkuv

    integer :: n = 100000

    real(8) :: x, xmin, xmax, dx
    real(8) :: Si, Ci
    real(8) :: Si2, Ci2

    real(8) :: Tc, Ts, Tst, IRD2bx2

    ! write(*, *) Ts(1.0d-15, 1.0d0, 1.0d0)
    ! stop

    ! call fill_array_Si(array_x_Si, array_Si, n_array_Si)
    ! call fill_array_Ci(array_x_Ci, array_Ci, n_array_Ci)

    ! x = 1.0d16
    ! write(*, *) Si(x), Si2(x)
    ! write(*, *) Ci(x), Ci2(x)
    !
    ! stop

    ! open(unit=11, file = 'Si.txt')
    !     xmin = -1.0d4
    !     xmax = 1.d4
    !     dx = (xmax - xmin)/real(n - 1, 8)
    !     do i = 1, n
    !         x = xmin + real(i - 1, 8)*dx
    !         write(11, "(4e25.16)") x, Si(x), Si2(x), Si(x) - Si2(x)
    !     end do
    ! close(11)

    ! open(unit=11, file = 'Ci.txt')
    !     xmin = 1.0d-10
    !     xmax = 1.d4
    !     dx = (xmax - xmin)/real(n - 1, 8)
    !     do i = 1, n
    !         x = xmin + real(i - 1, 8)*dx
    !         write(11, "(4e25.16)") x, Ci(x), Ci2(x), Ci(x) - Ci2(x)
    !     end do
    ! close(11)
    !
    ! stop


    open (unit = 11, file = 'Ps.txt', status = 'old')

    do i = 1, nPs
        read(11, *) (array_k(nPs - i + 1), array_Ps(nPs - i + 1))
    end do

    close(11)

    open(unit=11, file = 'OmegaGW0.txt')

    kmin = array_k(1)
    kmax = array_k(nPs)

    ftoHzinitial = ftoHzk(kmin)
    ftoHzfinal = ftoHzk(kmax)

    do i = 1, ni
        result = 0.0d0
        ftoHzi = 10.0d0**((real(i - 1, 8)*(log10(ftoHzfinal) - &
        & log10(ftoHzinitial)))/real(ni - 1, 8) + &
        & log10(ftoHzinitial))
        ki = kftoHz(ftoHzi)
        vmini = kmin/ki
        vmaxi = kmax/ki
        do j = 1, nj
            vj = 10.0d0**((real(j - 1, 8)*(log10(vmaxi) - log10(vmini)))/ &
            & real(nj - 1, 8) + log10(vmini))
            vjp1 = 10.0d0**((real(j, 8)*(log10(vmaxi) - log10(vmini)))/ &
            & real(nj - 1, 8) + log10(vmini))
            dvj = vjp1 - vj
            uminj = abs(1.0d0 - vj)
            umaxj = abs(1.0d0 + vj)
            do l = 1, nl
                ul = (real(l - 1, 8)*(umaxj - uminj))/real(nl - 1, 8) + uminj
                ulp1 = (real(l, 8)*(umaxj - uminj))/real(nl - 1, 8) + uminj
                dul = ulp1 - ul
                result = result + dvj*dul*fkuv(ki, ul, vj)
            end do ! l
        end do ! j
        OmegaGW0i = abs(result)
        write(11, "(2e25.16)") ftoHzi, OmegaGW0i
    end do ! i

    close(11)

end program OmegaGW0_alpha_attractor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function kftoHz(ftoHz)

    implicit none

    real(8) :: kftoHz
    real(8) :: ftoHz

    kftoHz = 1.6982007371081617d-42*ftoHz

end function kftoHz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function ftoHzk(k)

    implicit none

    real(8) :: ftoHzk
    real(8) :: k

    ftoHzk = 5.888585360661683d41*k

end function ftoHzk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function IRD2bx2(u, v)

    implicit none

    real(8), parameter :: pi = acos(-1.0d0)

    real(8) :: IRD2bx2
    real(8) :: u, v

    real(8) :: Tst, Tc, Theta

! {
! approximation
    IRD2bx2 = (9.0d0*(-3.0d0 + u**2 + v**2)**2* &
    & ((-4.0d0*u*v + (-3.0d0 + u**2 + v**2)* &
    & log(abs((3.0d0 - (u + v)**2)/ &
    & (3.0d0 - (u - v)**2))))**2 + &
    & pi**2*(-3.0d0 + u**2 + v**2)**2* &
    & Theta(-sqrt(3.0d0) + u + v)))/(32.0d0*u**6*v**6)
! }

! {
! exact
    ! IRD2bx2 = (Tst(u, v, 1.0d0)**2/81.0d0 + &
    ! & (Tc(u, v, 1.0d0)/9.0d0 + &
    ! & (3.0d0*pi*(-3.0d0 + u**2 + v**2)**2* &
    ! & Theta(-sqrt(3.0d0) + u + v))/(4.0d0*u**3*v**3))**2)/2.0d0
! }

end function IRD2bx2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Tc(u, v, x)

    implicit none

    real(8) :: Tc
    real(8) :: u, v, x

    real(8) :: Si

    ! u = 1.0d0
    ! v = 1.0d0
    ! x = 1.0d0

    Tc = (-27.0d0*(-3.0d0 + u**2 + v**2)**2* &
    &    (Si((1.0d0 - (u - v)/sqrt(3.0d0))*x) + &
    &    Si((1.0d0 + (u - v)/sqrt(3.0d0))*x) - &
    &    Si((1.0d0 - (u + v)/sqrt(3.0d0))*x) - &
    &    Si((1.0d0 + (u + v)/sqrt(3.0d0))*x)))/ &
    &    (4.0d0*u**3*v**3) - &
    &    (27.0d0*(-48.0d0*u*v*x**2*cos((u*x)/sqrt(3.0d0))* &
    &    cos((v*x)/sqrt(3.0d0))*(x*cos(x) + 3*sin(x)) &
    &    + 24.0d0*x*(-6.0d0 + (3.0d0 - u**2 - v**2)*x**2)* &
    &    cos(x)*sin((u*x)/sqrt(3.0d0))* &
    &    sin((v*x)/sqrt(3.0d0)) + &
    &    24.0d0*(-18.0d0 + (3.0d0 + u**2 + v**2)*x**2)*sin(x)* &
    &    sin((u*x)/sqrt(3.0d0))*sin((v*x)/sqrt(3.0d0)) + &
    &    48.0d0*sqrt(3.0d0)*x**2*cos(x)* &
    &    (v*cos((v*x)/sqrt(3.0d0))* &
    &    sin((u*x)/sqrt(3.0d0)) + &
    &    u*cos((u*x)/sqrt(3.0d0))*sin((v*x)/sqrt(3.0d0))) &
    &    + 8.0d0*sqrt(3.0d0)*x*sin(x)* &
    &    (v*(18.0d0 - (3.0d0 + u**2 - v**2)*x**2)* &
    &    cos((v*x)/sqrt(3.0d0))*sin((u*x)/sqrt(3.0d0)) &
    &    + u*(18.0d0 - (3.0d0 - u**2 + v**2)*x**2)* &
    &    cos((u*x)/sqrt(3.0d0))*sin((v*x)/sqrt(3.0d0)))) &
    &    )/(8.0d0*u**3*v**3*x**4)

    ! write(*, *) Tc
    ! stop

end function Tc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Ts(u, v, x)

    implicit none

    real(8) :: Ts
    real(8) :: u, v, x

    real(8) :: Ci

    ! u = 1.0d0
    ! v = 1.0d0
    ! x = 1.0d0



    Ts = (-27.0d0*(-3.0d0 + u**2 + v**2))/(u**2*v**2) + &
    &    (27.0d0*(-3.0d0 + u**2 + v**2)**2* &
    &    (Ci( abs( (1.0d0 - (u - v)/sqrt(3.0d0))*x ) ) + &
    &    Ci( abs( (1.0d0 + (u - v)/sqrt(3.0d0))*x ) ) - &
    &    Ci( abs( (1.0d0 + (u + v)/sqrt(3.0d0))*x ) ) - &
    &    Ci( abs( x*abs(1.0d0 - (u + v)/sqrt(3.0d0)) ) ) + &
    &    log(abs((3.0d0 - (u + v)**2)/(3.0d0 - (u - v)**2))) &
    &    ))/(4.0d0*u**3*v**3) + &
    &    (27.0d0*(48.0d0*u*v*x**2*cos((u*x)/sqrt(3.0d0))* &
    &    cos((v*x)/sqrt(3.0d0))*(-3.0d0*cos(x) + x*sin(x)) &
    &    + 24.0d0*(-18.0d0 + (3.0d0 + u**2 + v**2)*x**2)* &
    &    cos(x)*sin((u*x)/sqrt(3.0d0))* &
    &    sin((v*x)/sqrt(3.0d0)) + &
    &    24.0d0*x*(6.0d0 - (3.0d0 - u**2 - v**2)*x**2)*sin(x)* &
    &    sin((u*x)/sqrt(3.0d0))*sin((v*x)/sqrt(3.0d0)) - &
    &    48.0d0*sqrt(3.0d0)*x**2*sin(x)* &
    &    (v*cos((v*x)/sqrt(3.0d0))* &
    &    sin((u*x)/sqrt(3.0d0)) + &
    &    u*cos((u*x)/sqrt(3.0d0))*sin((v*x)/sqrt(3.0d0))) &
    &    + 8.0d0*sqrt(3.0d0)*x*cos(x)* &
    &    (v*(18.0d0 - (3.0d0 + u**2 - v**2)*x**2)* &
    &    cos((v*x)/sqrt(3.0d0))*sin((u*x)/sqrt(3.0d0)) &
    &    + u*(18.0d0 - (3.0d0 - u**2 + v**2)*x**2)* &
    &    cos((u*x)/sqrt(3.0d0))*sin((v*x)/sqrt(3.0d0)))) &
    &    )/(8.0d0*u**3*v**3*x**4)

    ! write(*, *) Ts
    ! stop

end function Ts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Si(x)

    implicit none

    real(8) :: Si
    real(8) :: x

! {
    ! integer, parameter :: n = 10000
    ! real(8), dimension(n) :: array_f
    ! real(8) :: ymin, ymax, dy, y
    ! integer :: i
    ! real(8) :: simpson
    !
    ! ymin = 0.0d0
    ! ymax = x
    ! dy = (ymax - ymin)/real(n - 1, 8)
    ! do i = 1, n
    !     y = ymin + real(i - 1, 8)*dy
    !     if (y == 0.0d0) y = 1.0d-16
    !     array_f(i) = sin(y)/y
    ! end do
    !
    ! Si = simpson(array_f, n, dy)
! }

! {
    real(8), parameter :: pi = acos(-1.0d0)
    integer, parameter :: n_array_Si = 100000

    real(8), dimension(n_array_Si) :: array_x_Si, array_Si
    common /arrays_Si/ array_x_Si, array_Si

    real(8) :: linear_interpolation

    if (x < -1.0d3) then
        Si = -pi/2.0d0
        return
    else if (x > 1.0d3) then
        Si = pi/2.0d0
        return
    else
        Si = linear_interpolation(array_x_Si, array_Si, n_array_Si, x)
        return
    end if
! }

end function Si

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Ci(x)

    implicit none

    real(8) :: Ci
    real(8) :: x

! {
    ! integer, parameter :: n = 10000
    ! real(8), dimension(n) :: array_f
    ! real(8) :: ymin, ymax, dy, y
    ! integer :: i
    ! real(8) :: simpson
    !
    ! ymin = x
    ! ymax = 1.0d3
    ! dy = x/real(n - 1, 8)
    ! do i = 1, n
    !     y = ymin + real(i - 1, 8)*dy
    !     if (y == 0.0d0) y = 1.0d-16
    !     array_f(i) = cos(y)/y
    ! end do
    !
    ! Ci = -simpson(array_f, n, dy)
! }

! {
    integer, parameter :: n_array_Ci = 100000

    real(8), dimension(n_array_Ci) :: array_x_Ci, array_Ci
    common /arrays_Ci/ array_x_Ci, array_Ci

    real(8) :: linear_interpolation

    if (x < 1.0d-10) then
        Ci = array_Ci(1)
        return
    else if (x > 1.0d3) then
        Ci = 0.0d0
        return
    else
        Ci = linear_interpolation(array_x_Ci, array_Ci, n_array_Ci, x)
        return
    end if
! }

end function Ci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Tst(u, v, x)

    implicit none

    real(8) :: Tst
    real(8) :: u, v, x

    real(8) :: Ts

    Tst = (27.0d0*(-3.0d0 + u**2 + v**2))/(u**2*v**2) - &
    &   (27.0d0*(-3.0d0 + u**2 + v**2)**2* &
    &   log(abs((3.0d0 - (u + v)**2)/(3.0d0 - (u - v)**2))))/ &
    &   (4.0d0*u**3*v**3) + Ts(u, v, x)

end function Tst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Theta(x)

    implicit none

    real(8) :: Theta
    real(8) :: x

    if (x < 0.0d0) then
        Theta = 0.0d0
    else
        Theta = 1.0d0
    end if

end function Theta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function fkuv(k, u, v)

    implicit none

    real(8) :: fkuv
    real(8) :: k, u, v

    integer, parameter :: nPs = 1000

    real(8), dimension(nPs) :: array_k, array_Ps
    common /arrays_Ps/ array_k, array_Ps

    real(8) :: Omegar0
    real(8) :: gc
    real(8) :: IRD2bx2
    real(8) :: linear_interpolation

    ! arXiv:1907.11896
    ! h = 0.6727
    Omegar0 = 9.17d-5

    ! Omegar0h2 = Omegar0*h**2

    ! arXiv:1912.05927
    gc = 106.75d0
    ! Omegar0h2 = 4.2d-5

    ! arXiv:1912.05927
    fkuv = (0.83d0)*(gc/10.75d0)**(-1.0d0/3.0d0)*Omegar0*(1.0d0/6.0d0)* &
    &   ((4.0d0*v**2 - (1.0d0 - u**2 + v**2)**2)**2* &
    &   IRD2bx2(u, v)* &
    &   linear_interpolation(array_k, array_Ps, nPs, k*u)* &
    &   linear_interpolation(array_k, array_Ps, nPs, k*v) &
    &   )/(16.0d0*u**2*v**2)

end function fkuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function simpson(array_f, n, dx)

    implicit none

    real(8) :: simpson
    integer :: n
    real(8), dimension(n) :: array_f
    real(8) :: dx

    integer :: i
    real(8) :: s

    s = 0.0d0
    do i = 1, int(real(n - 1, 8)/2.0d0)
        s = s + array_f(2*i - 1) + 4.0d0*array_f(2*i) + array_f(2*i + 1)
    end do
    s = s/3.0d0
    if ( mod(n, 2) == 0 ) then
        s = s + (5.0d0*array_f(n) + 8.0d0*array_f(n - 1) - array_f(n - 2))/12.0d0
    end if

    simpson = dx*s

end function simpson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function linear_interpolation(arrayx, arrayy, imax, x)

    implicit none

    real(8) :: linear_interpolation
    integer :: imax
    real(8), dimension(imax) :: arrayx, arrayy
    real(8) :: x

    integer :: leftpoint, rightpoint, midpoint
    integer :: i

    leftpoint = 1
    rightpoint = imax
    do while ((rightpoint - leftpoint) > 1)
       midpoint = (rightpoint + leftpoint)/2
       if (arrayx(midpoint) > x) then
          rightpoint = midpoint
       else
          leftpoint = midpoint
       endif
    enddo
    i = leftpoint

    ! i = 1 + (real(imax - 1, 16)/(arrayx(imax) - arrayx(1)))*(x - arrayx(1))
    ! if(i < 1) then
    !     i = 1
    ! end if
    ! if(i > imax) then
    !     i = imax
    ! end if

    linear_interpolation = arrayy(i) + ((arrayy(i+1) - arrayy(i))/(arrayx(i+1) - arrayx(i)))*(x - arrayx(i))

end function linear_interpolation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fill_array_Si(array_x_Si, array_Si, n_array_Si)

    implicit none

    integer, intent(in) :: n_array_Si
    real(8), dimension(n_array_Si), intent(out) :: array_x_Si, array_Si
    real(8) :: xmin, xmax, dx
    real(8) :: s
    integer :: i
    real(8) :: xim1, xi

    xmin = -1.0d3
    xmax = 1.0d3
    dx = (xmax - xmin)/real(n_array_Si - 1, 8)
    ! s = -1.562225466889056d0
    s = -1.570233121968771d0

    i = 1
    array_x_Si(i) = xmin
    array_Si(i) = s

    do i = 2, n_array_Si
        xim1 = xmin + real(i - 2, 8)*dx
        if (xim1 == 0.0d0) xim1 = 1.0d-16
        xi = xmin + real(i - 1, 8)*dx
        if (xi == 0.0d0) xi = 1.0d-16
        s = s + 0.5d0*dx*(sin(xim1)/xim1 + sin(xi)/xi)
        array_x_Si(i) = xi
        array_Si(i) = s
    end do

end subroutine fill_array_Si

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fill_array_Ci(array_x_Ci, array_Ci, n_array_Ci)

    implicit none

    integer, intent(in) :: n_array_Ci
    real(8), dimension(n_array_Ci), intent(out) :: array_x_Ci, array_Ci
    real(8) :: xmin, xmax, dx
    real(8) :: s
    integer :: i
    real(8) :: xim1, xi

    xmin = 1.0d3
    xmax = 1.0d-10
    dx = (xmax - xmin)/real(n_array_Ci - 1, 8)
    ! s = -0.00514882514261049
    s = 0.0008263155110907115d0

    i = 1
    array_x_Ci(i) = xmin
    array_Ci(i) = s

    do i = 2, n_array_Ci
        xim1 = xmin + real(i - 2, 8)*dx
        xi = xmin + real(i - 1, 8)*dx
        s = s + 0.5d0*dx*(Cos(xim1)/xim1 + Cos(xi)/xi)
        array_x_Ci(n_array_Ci - i + 1) = xi
        array_Ci(n_array_Ci - i + 1) = s
    end do

end subroutine fill_array_Ci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Si2(x)

    implicit none

    real(8) :: Si2
    real(8) :: x

! {
    integer, parameter :: n_array_Si = 100000

    real(8), dimension(n_array_Si) :: array_x_Si, array_Si
    common /arrays_Si/ array_x_Si, array_Si

    real(8) :: linear_interpolation

    if (x < -1.0d3) then
        Si2 = array_Si(1)
    else if (x > 1.0d3) then
        Si2 = array_Si(n_array_Si)
    else
        Si2 = linear_interpolation(array_x_Si, array_Si, n_array_Si, x)
    end if
! }

end function Si2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Ci2(x)

    implicit none

    real(8) :: Ci2
    real(8) :: x

! {
    integer, parameter :: n_array_Ci = 100000

    real(8), dimension(n_array_Ci) :: array_x_Ci, array_Ci
    common /arrays_Ci/ array_x_Ci, array_Ci

    real(8) :: linear_interpolation

    if (x < 1.0d-16) then
        Ci2 = array_Ci(1)
    else if (x > 1.0d3) then
        ! Ci2 = array_Ci(n_array_Ci)
        Ci2 = 0.0d0
    else
        Ci2 = linear_interpolation(array_x_Ci, array_Ci, n_array_Ci, x)
    end if
! }

end function Ci2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
