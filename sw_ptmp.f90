function sw_ptmp(S,T,P,PR)

! SW_PTMP    Potential temperature
!===========================================================================
! SW_PTMP  $Id: sw_ptmp.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
!          Copyright (C) CSIRO, Phil Morgan 1992.
!
! USAGE:  ptmp = sw_ptmp(S,T,P,PR)
!
! DESCRIPTION:
!    Calculates potential temperature as per UNESCO 1983 report.
!
! INPUT:  (all must have same dimensions)
!   S  = salinity    [psu      (PSS-78) ]
!   T  = temperature [degree C (ITS-90)]
!   P  = pressure    [db]
!   PR = Reference pressure  [db]
!        (P & PR may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
!
! OUTPUT:
!   ptmp = Potential temperature relative to PR [degree C (ITS-90)]

    implicit none

    real(kind=8),intent(in) :: S,T,P,PR



    real(kind=8) :: del_P,del_th,th,q,sw_adtg,sw_ptmp
    
!------
! BEGIN
!------

! theta1
    del_P  = PR - P
    del_th = del_P*sw_adtg(S,T,P)
    th     = T * 1.00024d0 + 0.5d0*del_th
    q      = del_th

! theta2
    del_th = del_P*sw_adtg(S,th/1.00024d0,P+0.5d0*del_P)
    th     = th + (1.0d0- 1.0d0/sqrt(2.0d0))*(del_th - q)
    q      = (2.0d0-sqrt(2.0d0))*del_th + (-2.0d0+3.0d0/sqrt(2.0d0))*q

! theta3
    del_th = del_P*sw_adtg(S,th/1.00024d0,P+0.5d0*del_P)
    th     = th + (1.0d0 + 1.0d0/sqrt(2.0d0))*(del_th - q)
    q      = (2.0d0 + sqrt(2.0d0))*del_th + (-2.0d0-3.0d0/sqrt(2.0d0))*q

! theta4
    del_th = del_P*sw_adtg(S,th/1.00024d0,P+del_P)
    sw_ptmp     = (th + (del_th - 2.0d0*q)/6.0d0)/1.00024d0

    end function sw_ptmp