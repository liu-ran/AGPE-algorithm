function sw_adtg(S,T,P)

! SW_ADTG    Adiabatic temperature gradient
!===========================================================================
! SW_ADTG   $Id: sw_adtg.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
!           Copyright (C) CSIRO, Phil Morgan  1992.
!
! adtg = sw_adtg(S,T,P)
!
! DESCRIPTION:
!    Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
!
! INPUT:  (all must have same dimensions)
!   S = salinity    [psu      (PSS-78) ]
!   T = temperature [degree C (ITS-90)]
!   P = pressure    [db]
!       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
!
! OUTPUT:
!   ADTG = adiabatic temperature gradient [degree_C/db]

    
    implicit none

    real(kind=8),intent(in) :: S,T,P



    real(kind=8) :: T68,a0,a1,a2,a3,b0,b1,c0,c1,c2,c3,d0,d1,e0,e1,e2,sw_adtg




!-------------
! BEGIN
!-------------

    T68 = 1.00024d0 * T;

    a0 =  3.5803d-5;
    a1 = +8.5258d-6;
    a2 = -6.836d-8;
    a3 =  6.6228d-10;

    b0 = +1.8932d-6;
    b1 = -4.2393d-8;

    c0 = +1.8741d-8;
    c1 = -6.7795d-10;
    c2 = +8.733d-12;
    c3 = -5.4481d-14;

    d0 = -1.1351d-10;
    d1 =  2.7759d-12;

    e0 = -4.6206d-13;
    e1 = +1.8676d-14;
    e2 = -2.1687d-16;

    sw_adtg = a0+(a1+(a2+a3*T68)*T68)*T68+ (b0+b1*T68)*(S-35.0)+ ((c0+(c1+(c2+c3*T68)*T68)*T68)+ (d0+d1*T68)*(S-35.0))*P+ (e0+(e1+e2*T68)*T68)*P*P

    end function sw_adtg
         
     
     
         