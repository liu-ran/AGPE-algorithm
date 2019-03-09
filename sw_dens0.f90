function sw_dens0(S,T)

implicit none

! SW_DENS0   Denisty of sea water at atmospheric pressure
!=========================================================================
! SW_DENS0  $Id: sw_dens0.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
!           Copyright (C) CSIRO, Phil Morgan 1992
!
! USAGE:  dens0 = sw_dens0(S,T)
!
! DESCRIPTION:
!    Density of Sea Water at atmospheric pressure using
!   UNESCO 1983 (EOS 1980) polynomial.
!
! INPUT:  (all must have same dimensions)
!   S = salinity    [psu      (PSS-78)]
!   T = temperature [degree C (ITS-90)]
!
! OUTPUT:
!   dens0 = density  [kg/m^3] of salt water with properties S,T,
!           P=0 (0 db gauge pressure)

    real(kind=8),intent(in) :: S,T


    real(kind=8) :: T68,b0,b1,b2,b3,b4,c0,c1,c2,d0,sw_dens0
    real(kind=8) :: sw_smow


    T68 = T * 1.00024d0

    !     UNESCO 1983 eqn(13) p17.

    b0 =  8.24493d-1
    b1 = -4.0899d-3
    b2 =  7.6438d-5
    b3 = -8.2467d-7
    b4 =  5.3875d-9

    c0 = -5.72466d-3
    c1 = +1.0227d-4
    c2 = -1.6546d-6

    d0 = 4.8314d-4
    sw_dens0 = sw_smow(T) + (b0 + (b1 + (b2 + (b3 + b4*T68)*T68)*T68)*T68)*S + (c0 + (c1 + c2*T68)*T68)*S*sqrt(S) + d0*S**2.0
    
    end function sw_dens0       
    
    