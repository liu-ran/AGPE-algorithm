function sw_smow(T)

! SW_SMOW    Denisty of standard mean ocean water (pure water)
!=========================================================================
! SW_SMOW  $Id: sw_smow.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
!          Copyright (C) CSIRO, Phil Morgan 1992.
!
! USAGE:  dens = sw_smow(T)
!
! DESCRIPTION:
!    Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.
!
! INPUT:
!   T = temperature [degree C (ITS-90)]
!
! OUTPUT:
!   dens = density  [kg/m^3]

    implicit none

    real(kind=8),intent(in) :: T


    real(kind=8) :: T68,a0,a1,a2,a3,a4,a5,sw_smow



    a0 = 999.842594d0
    a1 =   6.793952d-2
    a2 =  -9.095290d-3
    a3 =   1.001685d-4
    a4 =  -1.120083d-6
    a5 =   6.536332d-9

    T68 = T * 1.00024d0
    sw_smow = a0 + (a1 + (a2 + (a3 + (a4 + a5*T68)*T68)*T68)*T68)*T68

    end function sw_smow