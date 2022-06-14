!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Kazem Rezazadeh                                   June 14, 2022
! kazem.rezazadeh@ipm.ir
! https://github.com/krezazadeh
!
! This program calculates the power spectrum of the scalar
! perturbations which is noted by P_s, in the setup of
! alpha-attractor inflation with a tiny bump. For this purpose, it
! first computes the background dynamics, and stores the results in
! the file "phi.txt". Then it solves the Mukhanov-Sasaki equation
! by using the Bunch-Davies initial conditions and computes the
! scalar power spectrum P_s. The results of P_s are stored in the
! file "Ps.txt".
!
! To cite this code, please cite the following paper:
! https://arxiv.org/abs/2110.01482
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program Ps_alpha_attractor

    implicit none

    real(8), parameter :: pi = acos(-1.0d0)

    real(8), parameter :: MP = 1.0d0
    real(8), parameter :: kg = 2.3034255543767065d8*MP
    real(8), parameter :: Mpc = 3.8082082148306456d56/MP

    integer, parameter :: nNe = 1000
    integer, parameter :: nk = 1000
    integer, parameter :: nNe2 = 10000
    integer :: i, j

    real(8) :: V0, AA, sigma, phi0
    common /potential_parameters/ V0, AA, sigma, phi0

    real(8) :: Nei, Nef, dNe
    real(8) :: phii, phipi
    real(8) :: Ne, phi, phip

    real(8), dimension(1000) :: array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe
    common /background_arrays/ array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe

    real(8) :: kphi1, kphi2, kphi3, kphi4
    real(8) :: kphip1, kphip2, kphip3, kphip4
    real(8) :: dphidNe, dphipdNe
    real(8) :: d2phipdNe2
    real(8) :: a, z
    real(8) :: H, dHdNe, d2HdNe2
    real(8) :: epsilon_H
    real(8) :: ns_sr, r_sr
    real(8) :: linear_interpolation
    real(8) :: NeBD
    real(8) :: k
    real(8) :: NeEvaluation
    complex(8) :: uki, ukpi
    complex(8) :: ukBD, dukBDdNe
    real(8) :: Ne_i
    complex(8) :: uk_i, ukp_i
    complex(8) :: kuk1, kuk2, kuk3, kuk4
    complex(8) :: kukp1, kukp2, kukp3, kukp4
    complex(8) :: dukdNe, dukpdNe
    real(8) :: fk, dzdNe, d2zdNe2
    real(8) :: k60

    integer :: nPs
    common /nPs/ nPs

    real(8) :: k_j, Ps_j

    real(8), dimension(1000) :: array_k, array_Ps
    common /arrays_Ps/ array_k, array_Ps

    V0 = 1.448d-10
    AA = 2.0448801d-3
    sigma = 2.524999d-2
    phi0 = 4.850001d0

    Nei = 80.0d0
    Nef = 0.0d0
    dNe = (Nef - Nei)/real(nNe - 1, 8)

    phii = 6.41032d0
    phipi = 0.0143d0

    open(unit=11, file='phi.txt')

    i = 1
    Ne = Nei
    phi = phii
    phip = phipi
    array_Ne(i) = Ne
    array_phi(i) = phi
    array_phip(i) = phip
    write(11, "(3e25.16)") Ne, phi, phip

    do i = 2, nNe
        kphi1 = dNe*dphidNe(Ne, phi, phip)
        kphip1 = dNe*dphipdNe(Ne, phi, phip)
        kphi2 = dNe*dphidNe(Ne + dNe/2.0d0, phi + kphi1/2.0d0, phip + kphip1/2.0d0)
        kphip2 = dNe*dphipdNe(Ne + dNe/2.0d0, phi + kphi1/2.0d0, phip + kphip1/2.0d0)
        kphi3 = dNe*dphidNe(Ne + dNe/2.0d0, phi + kphi2/2.0d0, phip + kphip2/2.0d0)
        kphip3 = dNe*dphipdNe(Ne + dNe/2.0d0, phi + kphi2/2.0d0, phip + kphip2/2.0d0)
        kphi4 = dNe*dphidNe(Ne + dNe, phi + kphi3, phip + kphip3)
        kphip4 = dNe*dphipdNe(Ne + dNe, phi + kphi3, phip + kphip3)
        Ne = Ne + dNe
        phi = phi + (1.0d0/6.0d0)*(kphi1 + 2.0d0*kphi2 + 2.0d0*kphi3 + kphi4)
        phip = phip + (1.0d0/6.0d0)*(kphip1 + 2.0d0*kphip2 + 2.0d0*kphip3 + kphip4)
        array_Ne(i) = Ne
        array_phi(i) = phi
        array_phip(i) = phip
        array_H(i) = H(Ne, phi, phip)
        array_dHdNe(i) = dHdNe(Ne, phi, phip)
        array_dphipdNe(i) = dphipdNe(Ne, phi, phip)
        array_d2phipdNe2(i) = d2phipdNe2(Ne, phi, phip)
        ! write(11, "(3e25.16)") Ne, phi, phip
        write(11, "(3e25.16)") Ne, phi, epsilon_H(Ne, phi, phip)
    end do

    close(11)

    phi = linear_interpolation(array_Ne, array_phi, nNe, 60.0d0)
    write(*, "(1a25, 1es25.16)") "phi(Ne = 60) = ", phi

    write(*, "(1a25, 1es25.16)") "ns_sr(Ne = 60) = ", ns_sr(phi)
    write(*, "(1a25, 1es25.16)") "r_sr(Ne = 60) = ", r_sr(phi)


    write(*, "(1a25, 1es25.16)") "phi(Ne = 0) = ", array_phi(nNe)
    write(*, "(1a25, 1es25.16)") "epsilon_H(Ne = 0) = ", array_dHdNe(nNe)/array_H(nNe)

    open(unit=11, file='Ps.txt')

    k60 = a(60.0d0)*linear_interpolation(array_Ne, array_H, nNe, 60.0d0)
    nPs = 0

    do j = 1, nk
        Ne = 60.0d0*real(j - 1, 16)/real(nk - 1, 8)
        NeBD = Ne + 5.0d0
        k = a(Ne)*linear_interpolation(array_Ne, array_H, nNe, Ne)
        if(Ne >= 5.0d0) then
            NeEvaluation = Ne - 5.0d0
        else
            NeEvaluation = 0.0d0
        end if

        Nei = NeBD
        Nef = NeEvaluation
        dNe = (Nef - Nei)/real(nNe2 - 1, 8)

        uki = ukBD(k)
        ukpi = dukBDdNe(k, NeBD)

        i = 1

        Ne_i = Nei
        uk_i = uki
        ukp_i = ukpi

        do i = 2, nNe2
            kuk1 = dNe*dukdNe(Ne_i, uk_i, ukp_i)
            kukp1 = dNe*dukpdNe(k, Ne_i, uk_i, ukp_i)
            kuk2 = dNe*dukdNe(Ne_i + dNe/2.0d0, uk_i + kuk1/2.0d0, ukp_i + kukp1/2.0d0)
            kukp2 = dNe*dukpdNe(k, Ne_i + dNe/2.0d0, uk_i + kuk1/2.0d0, ukp_i + kukp1/2.0d0)
            kuk3 = dNe*dukdNe(Ne_i + dNe/2.0d0, uk_i + kuk2/2.0d0, ukp_i + kukp2/2.0d0)
            kukp3 = dNe*dukpdNe(k, Ne_i + dNe/2.0d0, uk_i + kuk2/2.0d0, ukp_i + kukp2/2.0d0)
            kuk4 = dNe*dukdNe(Ne_i + dNe, uk_i + kuk3, ukp_i + kukp3)
            kukp4 = dNe*dukpdNe(k, Ne_i + dNe, uk_i + kuk3, ukp_i + kukp3)
            Ne_i = Ne_i + dNe
            uk_i = uk_i + (1.0d0/6.0d0)*(kuk1 + 2.0d0*kuk2 + 2.0d0*kuk3 + kuk4)
            ukp_i = ukp_i + (1.0d0/6.0d0)*(kukp1 + 2.0d0*kukp2 + 2.0d0*kukp3 + kukp4)
        end do

        k_j = k*((0.05d0/Mpc)/k60)
        Ps_j = (k**3*abs(uk_i)**2)/(2.0d0*pi**2*abs(z(NeEvaluation))**2)
        if(.not.(isnan(k_j) .or. isnan(Ps_j))) then
            nPs = nPs + 1
            array_k(nPs) = k_j
            array_Ps(nPs) = Ps_j
            write(11, "(2e25.16)") k_j, Ps_j
        end if

    end do

    close(11)

    write(*, *)
    write(*, "(1a25, 1es25.16)") "Ps(Ne = 60) = ", array_Ps(nPs)

    write(*, "(1a25, 1es25.16)") "k_peak/Mpc^-1 = ", array_k(maxloc(array_Ps))*Mpc
    write(*, "(1a25, 1es25.16)") "Ps_peak = ", maxval(array_Ps)

end program Ps_alpha_attractor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function V(phi)

    implicit none

    real(8) :: V
    real(8) :: phi

    real(8) :: V0, AA, sigma, phi0
    common /potential_parameters/ V0, AA, sigma, phi0

    real(8) :: sech

    V = V0*(1.0d0 + AA*sech((phi - phi0)/sigma)**2)*tanh(phi/sqrt(6.0d0))**2

end function V

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dVdphi(phi)

    implicit none

    real(8) :: dVdphi
    real(8) :: phi

    real(8) :: V0, AA, sigma, phi0
    common /potential_parameters/ V0, AA, sigma, phi0

    real(8) :: sech

    dVdphi = sqrt(2.0d0/3.0d0)*V0* &
     &   sech(phi/sqrt(6.0d0))**2* &
     &   (1.0d0 + AA*sech((phi - phi0)/sigma)**2)* &
     &   tanh(phi/sqrt(6.0d0)) - &
     &  (2.0d0*AA*V0*sech((phi - phi0)/sigma)**2* &
     &     tanh(phi/sqrt(6.0d0))**2*tanh((phi - phi0)/sigma) &
     &     )/sigma

end function dVdphi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function d2Vdphi2(phi)

    implicit none

    real(8) :: d2Vdphi2
    real(8) :: phi

    real(8) :: V0, AA, sigma, phi0
    common /potential_parameters/ V0, AA, sigma, phi0

    real(8) :: sech

    d2Vdphi2 = (V0*sech(phi/sqrt(6.0d0))**4* &
     &     (1.0d0 + AA*sech((phi - phi0)/sigma)**2))/3.0d0 - &
     &  (2.0d0*AA*V0*sech((phi - phi0)/sigma)**4* &
     &     tanh(phi/sqrt(6.0d0))**2)/sigma**2 - &
     &  (2.0d0*V0*sech(phi/sqrt(6.0d0))**2* &
     &     (1.0d0 + AA*sech((phi - phi0)/sigma)**2)* &
     &     tanh(phi/sqrt(6.0d0))**2)/3.0d0 - &
     &  (4.0d0*sqrt(2.0d0/3.0d0)*AA*V0* &
     &     sech(phi/sqrt(6.0d0))**2* &
     &     sech((phi - phi0)/sigma)**2* &
     &     tanh(phi/sqrt(6.0d0))*tanh((phi - phi0)/sigma))/ &
     &   sigma + (4.0d0*AA*V0*sech((phi - phi0)/sigma)**2* &
     &     tanh(phi/sqrt(6.0d0))**2* &
     &     tanh((phi - phi0)/sigma)**2)/sigma**2

end function d2Vdphi2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function H(Ne, phi, phip)

    implicit none

    real(8) :: H
    real(8) :: Ne, phi, phip
    real(8) :: V

    H = sqrt(V(phi)/(3.0d0 - phip**2/2.0d0))

end function H

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dHdNe(Ne, phi, phip)

    implicit none

    real(8) :: dHdNe
    real(8) :: Ne, phi, phip
    real(8) :: H

    dHdNe = (H(Ne, phi, phip)*phip**2)/2.0d0

end function dHdNe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function d2HdNe2(Ne, phi, phip)

    implicit none

    real(8) :: d2HdNe2
    real(8) :: Ne, phi, phip
    real(8) :: H, dHdNe, dphipdNe

    d2HdNe2 = (phip*(phip*dHdNe(Ne, phi, phip) + &
    &      2.0d0*H(Ne, phi, phip)*dphipdNe(Ne, phi, phip)))/2.0d0

end function d2HdNe2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function epsilon_H(Ne, phi, phip)

    implicit none

    real(8) :: epsilon_H
    real(8) :: Ne, phi, phip
    real(8) :: H, dHdNe

    epsilon_H = dHdNe(Ne, phi, phip)/H(Ne, phi, phip)

end function epsilon_H

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dphidNe(Ne, phi, phip)

    implicit none

    real(8) :: dphidNe
    real(8) :: Ne, phi, phip

    dphidNe = phip

end function dphidNe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dphipdNe(Ne, phi, phip)

    implicit none

    real(8) :: dphipdNe
    real(8) :: Ne, phi, phip
    real(8) :: dVdphi, H, dHdNe

    dphipdNe = (-dVdphi(phi) + phip*H(Ne, phi, phip)*(-dHdNe(Ne, phi, phip) + &
    & 3.0d0*H(Ne, phi, phip)))/H(Ne, phi, phip)**2

end function dphipdNe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function d2phipdNe2(Ne, phi, phip)

    implicit none

    real(8) :: d2phipdNe2
    real(8) :: Ne, phi, phip
    real(8) :: dVdphi, d2Vdphi2, H, dHdNe, dphipdNe, d2HdNe2

    d2phipdNe2 = (2.0d0*dHdNe(Ne, phi, phip)* &
     &     dVdphi(phi))/H(Ne, phi, phip)**3 + &
     &  3.0d0*dphipdNe(Ne, phi, phip) - &
     &  (phip*d2HdNe2(Ne, phi, phip) + &
     &     dHdNe(Ne, phi, phip)*dphipdNe(Ne, phi, phip))/ &
     &   H(Ne, phi, phip) + (phip* &
     &     (dHdNe(Ne, phi, phip)**2 - &
     &       d2Vdphi2(phi)))/H(Ne, phi, phip)**2

end function d2phipdNe2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function sech(x)

    implicit none

    real(8) :: sech
    real(8) :: x

    sech = 1.0d0/cosh(x)

end function sech

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function csch(x)

    implicit none

    real(8) :: csch
    real(8) :: x

    csch = 1.0d0/sinh(x)

end function csch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function ns_sr(phi)

    implicit none

    real(8) :: ns_sr
    real(8) :: phi

    real(8) :: csch

    ns_sr = ((-11.0d0 + 3.0d0*cosh(sqrt(2.0d0/3.0d0)*phi))*csch(phi/sqrt(6.0d0))**2)/6.0d0

end function ns_sr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function r_sr(phi)

    implicit none

    real(8) :: r_sr
    real(8) :: phi

    real(8) :: csch

    r_sr = (64.0d0*csch(sqrt(2.0d0/3.0d0)*phi)**2)/3.0d0

end function r_sr

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
       if (arrayx(midpoint) < x) then
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

function a(Ne)

    implicit none

    real(8) :: a
    real(8) :: Ne

    a = exp(50.0d0 - Ne)

end function a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dadNe(Ne)

    implicit none

    real(8) :: dadNe
    real(8) :: Ne

    dadNe = -exp(50.0d0 - Ne)

end function dadNe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function d2adNe2(Ne)

    implicit none

    real(8) :: d2adNe2
    real(8) :: Ne

    d2adNe2 = exp(50.0d0 - Ne)

end function d2adNe2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function z(Ne)

    implicit none

    real(8) :: z
    real(8) :: Ne
    real(8) :: a, linear_interpolation

    real(8), dimension(1000) :: array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe
    common /background_arrays/ array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe

    z = -a(Ne)*linear_interpolation(array_Ne, array_phip, 1000, Ne)

end function z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dzdNe(Ne)

    implicit none

    real(8) :: dzdNe
    real(8) :: Ne
    real(8) :: a, dadNe
    real(8) :: linear_interpolation

    real(8), dimension(1000) :: array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe
    common /background_arrays/ array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe

    dzdNe = -linear_interpolation(array_Ne, array_phip, 1000, Ne)*dadNe(Ne) - &
    & a(Ne)*linear_interpolation(array_Ne, array_dphipdNe, 1000, Ne)

end function dzdNe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function d2zdNe2(Ne)

    implicit none

    real(8) :: d2zdNe2
    real(8) :: Ne
    real(8) :: a, dadNe, d2adNe2

    real(8), dimension(1000) :: array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe
    common /background_arrays/ array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe

    real(8) :: linear_interpolation

    d2zdNe2 = -2.0d0*dadNe(Ne)* &
     &   linear_interpolation(array_Ne, array_dphipdNe, 1000, Ne) - &
     &  linear_interpolation(array_Ne, array_phip, 1000, Ne)*d2adNe2(Ne) - &
     &  a(Ne)*linear_interpolation(array_Ne, array_d2phipdNe2, 1000, Ne)

end function d2zdNe2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function fk(Ne)

    implicit none

    real(8) :: fk
    real(8) :: Ne
    real(8) :: a, dadNe, z, dzdNe, d2zdNe2
    real(8) :: linear_interpolation

    real(8), dimension(1000) :: array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe
    common /background_arrays/ array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe

    fk = ((dadNe(Ne)*dzdNe(Ne))/ &
     &     a(Ne) + (linear_interpolation(array_Ne, array_dHdNe, 1000, Ne)* &
     &       dzdNe(Ne))/linear_interpolation(array_Ne, array_H, 1000, Ne) + &
     &    d2zdNe2(Ne))/z(Ne)

end function fk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function ukBD(k)

    implicit none

    complex(8) :: ukBD
    real(8) :: k

    ukBD = 1.0d0/sqrt(2.0d0*k)

end function ukBD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dukBDdNe(k, NeBD)

    implicit none

    complex(8) :: dukBDdNe
    real(8) :: k, NeBD
    real(8) :: a
    real(8) :: linear_interpolation

    real(8), dimension(1000) :: array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe
    common /background_arrays/ array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe

    dukBDdNe = ((0.0d0,1.0d0)*sqrt(k/2.0d0))/(a(NeBD)*linear_interpolation(array_Ne, array_H, 1000, NeBD))

end function dukBDdNe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dukdNe(Ne, uk, ukp)

    implicit none

    complex(8) :: dukdNe
    real(8) :: Ne
    complex(8) :: uk, ukp

    dukdNe = ukp

end function dukdNe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dukpdNe(k, Ne, uk, ukp)

    implicit none

    complex(8) :: dukpdNe
    real(8) :: k, Ne
    complex(8) :: uk, ukp
    real(8) :: a, dadNe, fk
    real(8) :: linear_interpolation
    real(8) :: kaH

    real(8), dimension(1000) :: array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe
    common /background_arrays/ array_Ne, array_phi, array_phip, array_dphipdNe, array_d2phipdNe2, array_H, array_dHdNe

    kaH = k/(a(Ne)*linear_interpolation(array_Ne, array_H, 1000, Ne))

    dukpdNe = -((kaH**2 - fk(Ne))*uk) - &
     &  (dadNe(Ne)*ukp)/ &
     &   a(Ne) - (linear_interpolation(array_Ne, array_dHdNe, 1000, Ne)* &
     &     ukp)/linear_interpolation(array_Ne, array_H, 1000, Ne)

end function dukpdNe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function MtoMSun(k)

    implicit none

    real(8) :: MtoMSun
    real(8) :: k
    real(8), parameter :: MP = 1.0d0
    real(8), parameter :: Mpc = 3.8082082148306456d56/MP
    real(8), parameter :: gamma = 0.2d0
    real(8), parameter :: gstar = 107.5d0

    MtoMSun = 3.68d0*(gamma/0.2d0)*((gstar/10.75d0)**(-1.0d0/6.0d0))/ &
    ((k*Mpc)/1.0d6)**2

end function MtoMSun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
