function sw_pres(depth,latitude)

implicit none

! SW_PRES    Pressure from depth

! INPUT:  (all must have same dimensions)
!   depth = depth [metres]
!   lat   = Latitude in decimal degress north [-90..+90]
!           (LAT may have dimensions 1x1 or 1xn where depth(mxn) )
!
! OUTPUT:
!  pres   = Pressure    [db]

    real(kind=8),intent(in) :: latitude,depth


    real(kind=8) :: DEG2RAD,X,C1,sw_pres
    real(kind=8),parameter  :: pi=3.1415926



    DEG2RAD = pi/180.
    X       = sin(abs(latitude)*DEG2RAD)  ! convert to radians
    C1      = 5.92d-3+X**2*5.25d-3
    sw_pres    = ((1.0d0-C1)-sqrt(((1.0d0-C1)**2)-(8.84d-6*depth)))/4.42d-6

    end function sw_pres