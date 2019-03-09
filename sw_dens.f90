function sw_dens(S,T,P)

implicit none

! SW_DENS    Density of sea water

! DESCRIPTION:
!    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.
!
! INPUT:  (all must have same dimensions)
!   S = salinity    [psu      (PSS-78)]
!   T = temperature [degree C (ITS-90)]
!   P = pressure    [db]
!       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
!
! OUTPUT:
!   dens = density  [kg/m^3]
    real(kind=8),intent(in) :: S,T,P


    real(kind=8) :: densP0,K,P10,sw_dens
    real(kind=8),parameter  :: pi=3.1415926
    real(kind=8) :: sw_dens0,sw_seck


    densP0 = sw_dens0(S,T)
    K      = sw_seck(S,T,P)
    P10      = P/10;  ! convert from db to atm pressure units
    sw_dens   = densP0/(1-P10/K)

    end function sw_dens