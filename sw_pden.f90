function sw_pden(S,T,P,PR)

! SW_PDEN    Potential density
!===========================================================================
! SW_PDEN  $Id: sw_pden.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
!          Copyright (C) CSIRO, Phil Morgan  1992.
!
! USAGE:  pden = sw_pden(S,T,P,PR)
!
! DESCRIPTION:
!    Calculates potential density of water mass relative to the specified
!    reference pressure by pden = sw_dens(S,ptmp,PR).
!
! INPUT:  (all must have same dimensions)
!   S  = salinity    [psu      (PSS-78) ]
!  T  = temperature [degree C (ITS-90)]
!   P  = pressure    [db]
!   PR = Reference pressure  [db]
!       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
!
! OUTPUT:
!   pden = Potential denisty relative to the ref. pressure [kg/m^3]
    
    implicit none

    real(kind=8),intent(in) :: S,T,P,PR



    real(kind=8) :: ptmp,sw_ptmp,sw_dens,sw_pden



    ptmp = sw_ptmp(S,T,P,PR);
    sw_pden = sw_dens(S,ptmp,PR);

    end function sw_pden