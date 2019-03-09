function sw_temp(S,T,P,PR)

! SW_TEMP    Temperature from potential temperature
! DESCRIPTION:
!    Calculates temperature from potential temperature at the reference
!    pressure PR and in-situ pressure P.
!
! INPUT:  (all must have same dimensions)
!   S     = salinity              [psu      (PSS-78) ]
!   PTMP  = potential temperature [degree C (ITS-90)]
!   P     = pressure              [db]
!   PR    = Reference pressure    [db]
!           (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
!
! OUTPUT:
!   temp = temperature [degree C (ITS-90)]

    implicit none

    real(kind=8),intent(in) :: S,T,P,PR
    real(kind=8) :: sw_temp,sw_ptmp

    sw_temp = sw_ptmp(S,T,PR,P)

    end function sw_temp
