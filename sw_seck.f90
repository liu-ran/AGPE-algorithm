function sw_seck(S,T,P)

! SW_SECK    Secant bulk modulus (K) of sea water
!=========================================================================
! SW_SECK  $Id: sw_seck.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
!          Copyright (C) CSIRO, Phil Morgan 1992.
!
! USAGE:  dens = sw_seck(S,T,P)
!
! DESCRIPTION:
!    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980.
!    UNESCO polynomial implementation.
!
! INPUT:  (all must have same dimensions)
!   S = salinity    [psu      (PSS-78) ]
!   T = temperature [degree C (ITS-90)]
!   P = pressure    [db]
!       (alternatively, may have dimensions 1*1 or 1*n where n is columns in S)
!
! OUTPUT:
!   K = Secant Bulk Modulus  [bars]

    implicit none

    real(kind=8),intent(in) :: S,T,P


    real(kind=8) :: P10,T68,h3,h2,h1,h0,AW,k2,k1,k0,BW,e4,e3,e2,e1,e0,KW,j0,i2,i1,i0,SR,A,m2,m1,m0,B,f3,f2,f1,f0,g2,g1,g0,KK0
    real(kind=8),parameter  :: pi=3.1415926
    real(kind=8) :: sw_seck


    !--------------------------------------------------------------------
! COMPUTE COMPRESSION TERMS
!--------------------------------------------------------------------
    P10 = P/10.0d0  !convert from db to atmospheric pressure units
    T68 = T * 1.00024d0

! Pure water terms of the secant bulk modulus at atmos pressure.
! UNESCO eqn 19 p 18

    h3 = -5.77905d-7
    h2 = +1.16092d-4
    h1 = +1.43713d-3
    h0 = +3.239908d0   ![-0.1194975];

    AW  = h0 + (h1 + (h2 + h3*T68)*T68)*T68

    k2 =  5.2787d-8
    k1 = -6.12293d-6
    k0 =  +8.50935d-5  ![+3.47718E-5];

    BW  = k0 + (k1 + k2*T68)*T68

    e4 = -5.155288d-5
    e3 = +1.360477d-2
    e2 = -2.327105d0
    e1 = +148.4206d0
    e0 = 19652.21d0    ![-1930.06];

    KW  = e0 + (e1 + (e2 + (e3 + e4*T68)*T68)*T68)*T68   ! eqn 19

!--------------------------------------------------------------------
! SEA WATER TERMS OF SECANT BULK MODULUS AT ATMOS PRESSURE.
!--------------------------------------------------------------------
    j0 = 1.91075d-4

    i2 = -1.6078d-6
    i1 = -1.0981d-5
    i0 =  2.2838d-3

    SR = sqrt(S)

    A  = AW + (i0 + (i1 + i2*T68)*T68 + j0*SR)*S


    m2 =  9.1697d-10
    m1 = +2.0816d-8
    m0 = -9.9348d-7

    B = BW + (m0 + (m1 + m2*T68)*T68)*S   ! eqn 18

    f3 =  -6.1670d-5
    f2 =  +1.09987d-2
    f1 =  -0.603459d0
    f0 = +54.6746d0

    g2 = -5.3009d-4
    g1 = +1.6483d-2
    g0 = +7.944d-2

    KK0 = KW + (  f0 + (f1 + (f2 + f3*T68)*T68)*T68 + (g0 + (g1 + g2*T68)*T68)*SR )*S    ! eqn 16
              

    sw_seck = KK0 + (A + B*P10)*P10;  !eqn 15

    end function sw_seck