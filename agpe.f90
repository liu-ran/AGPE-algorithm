!  agpe.f90 
!
!  FUNCTIONS:
!  agpe - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: agpe
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program agpe

   implicit none

    ! Variables
    integer :: imt,jmt,kmt
    real(kind=8) :: dpa,agpe_ave,ave_dmass
    real(kind=8),dimension(35) :: ave_dmass_all
    real(kind=8),dimension(380,329,50) :: t , s
    real(kind=8),dimension(380) :: lon
    real(kind=8),dimension(329) :: lat
    real(kind=8),dimension(50) :: lev
    character(len=23) :: savename,savename1
    character(len=23) :: t_nameall(48),s_nameall(48)
    integer :: i,nn,mk


    ! CMIP5 Body of AGPE_Lr
    open (unit=111,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_3dim\t_filename.dat',status='old',action='read')  !!3,4为rcp
    open (unit=112,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_3dim\s_filename.dat',status='old',action='read')
    read(111,*) t_nameall
    read(112,*) s_nameall
    !print *, t_nameall(17)
    !print *, s_nameall(17)
    close(111)
    close(112)

    open (unit=3,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_3dim\soda_lat.dat',status='old',action='read')
    open (unit=4,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_3dim\soda_lon.dat',status='old',action='read')
    open (unit=5,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_3dim\soda_lev.dat',status='old',action='read')
    read(3,*) lat
    read(4,*) lon
    read(5,*) lev
    close(3)
    close(4)
    close(5)

!    open (unit=6,file='D:\AGPE\CMIP5\cmip5_3dim\dmass.dat',status='old',action='read')
!    read(6,*) ave_dmass_all
!    close(6)

    do nn=19,30
        open (unit=1,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_t_dat\'//t_nameall(nn),status='old',action='read')
        open (unit=2,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_s_dat\'//s_nameall(nn),status='old',action='read')
        do i=1,125020
            read(1,*) t(ceiling(real(i)/329.),mod(i-1,329)+1,:)
            read(2,*) s(ceiling(real(i)/329.),mod(i-1,329)+1,:)
        end do
        jmt=size(t,1)
        imt=size(t,2)
        kmt=size(t,3)
        close(1) 
        close(2)
        savename1=t_nameall(nn)
        savename='soda_'//savename1(14:)
        print *, kmt,imt,jmt
        print *, savename

    !open (unit=1,file='t.dat',status='old',action='read')
    !open (unit=2,file='s.dat',status='old',action='read')
    !open (unit=3,file='lat.dat',status='old',action='read')
    !open (unit=4,file='lon.dat',status='old',action='read')
    !open (unit=5,file='lev.dat',status='old',action='read')
    !do i=1,64800
       !read(1,*) t(ceiling(i/180.),mod(i-1,180)+1,:)
       !read(2,*) s(ceiling(i/180.),mod(i-1,180)+1,:)
    !end do
    !jmt=size(t,1)
    !imt=size(t,2)
    !kmt=size(t,3)
    !read(3,*) lat
    !read(4,*) lon
    !read(5,*) lev
    !close(1) 
    !close(2)
    !close(3)
    !close(4)
    !close(5)

    !!open (unit=99,file='sst.dat',status='new',action='write')                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!! 测试存读
    !!do i=1,imt
    !!    write(99,'(360f12.4)') t(:,i,1)
    !!end do
    !!close (99)  

    !write (*,*) 'please enter input file name'
    !read (*,*) savename
    !write (*,1000) savename
    !1000 format (' ','the input file name is : ',A) 

    !write (*,*) 'please enter output dpa (my suggest is 10 ) '
    !read (*,*) dpa



    !write (*,*) 'please enter the average mass loss of bossinessq model  '
    !read (*,*) ave_dmass

    dpa=10.
    ave_dmass=0.0
    !ave_dmass=ave_dmass_all(nn)
    print *, dpa
    print *, ave_dmass
    call AGPE_CALCULATE(ave_dmass,savename,dpa,imt,jmt,kmt,lon,lat,lev,t,s,agpe_ave)

    end do

    end program agpe




   subroutine AGPE_CALCULATE(ave_dmass,savename,dpa,imt,jmt,kmt,lon,lat,lev,t,s,ave_agpe)

    implicit none
    !!!!!!!!!!!!!计算可利用重力势能
    !!!!!!!!!!!!! t(j,i,k)
    real(kind=8), intent(in) :: ave_dmass
    character(len=*),intent(in) :: savename
    integer,intent(in) :: imt,jmt,kmt
    real(kind=8),intent(in) :: dpa
    real(kind=8),intent(in),dimension(jmt) :: lon
    real(kind=8),intent(in),dimension(imt) :: lat
    real(kind=8),intent(in),dimension(kmt) :: lev
    real(kind=8),intent(in),dimension(jmt,imt,kmt) :: t
    real(kind=8),intent(in),dimension(jmt,imt,kmt) :: s
    real(kind=8),intent(out) :: ave_agpe
    real(kind=8),allocatable,dimension(:) :: dpb
    real(kind=8),allocatable,dimension(:) :: gpe1,pre_rho,mass1,s1,t1,pden1,pden1_o,pden1_r,pden1_c,h,volume_ref1
    real(kind=16),allocatable,dimension(:) :: dh,dp
    integer,allocatable,dimension(:) :: index_c,index1,index1_o
    integer,allocatable,dimension(:,:) :: index_o
    real(kind=8),dimension(kmt) :: z,dz,zt,zb,zzc,zbk,area
    logical, dimension(kmt):: mask
    real(kind=8),dimension(imt,jmt) :: dxyj,dxyj1
    real(kind=8),dimension(kmt,imt,jmt) :: volume,volume_ref,den,gpe,mass,agpe,gper,zr,denref
    integer :: water_num,ijkm,i,j,k,nmt,n,kk,kc,num1,jd,it(1),x,t_potential,jcheck
    real(kind=8) :: pb0,pres,sum_agpe,sum_volume,  test,dist, zzzzz,t_p
    real(kind=8),parameter :: g=9.81
    real(kind=8),parameter  :: r=6.371d6
    real(kind=8),parameter  :: pi=3.1415926
    real(kind=8) :: sw_pres,sw_dens,sw_pden,sw_temp                 !!!!!!!!!!!!!!!!!!!此行 定义自编函数

    t_potential=1
    !print *, ave_dmass
    !write (*,*) 'T is potential temp ? (1 is yes or 0 is no) '
    !read (*,*) t_potential
    !write (*,*) t_potential

    water_num=0
    do i=1,imt
        do j=1,jmt
            do k=1,kmt
                if (t(j,i,k)<9999.0d0 ) then
                    water_num=water_num+1
                end if
            end do
        end do
    end do
    ijkm=water_num

    z=lev
    dz(1)=0.5d0*(z(2)-z(1))+z(1)
    zt(1)=0.0d0;
    zb(1)=z(1)+0.5d0*(z(2)-z(1))
    do k=2,kmt
        if (k==kmt) then
            dz(k)=dz(k-1)          !单元水体的厚度
            zt(k)=0.5d0*(z(k)+z(k-1))!单元水体的顶点的深度
            zb(k)=zt(k)+dz(k);     !单元水体的底部的深度
        else
            dz(k)=0.5d0*(z(k+1)-z(k-1))
            zt(k)=0.5d0*(z(k)+z(k-1))
            zb(k)=zt(k)+dz(k)
        end if
    end do
    do k=1,kmt
        zzc(k)=zb(kmt)-0.5d0*(zt(k)+zb(k))  !每个水体单元重心到海底的高度
    end do
    do k=kmt,1,-1
        zbk(k)=zb(kmt)-zb(k)  !每个水体单元底部到海底的高度
    end do

    area=0.                             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!用来累加的变量，定义双精度后，初始值不是0

    do j=1,jmt
        do i=1,imt
            if (j>1 .and. i>1 .and. j<jmt .and. i<imt) then
                dxyj(i,j)=((lon(j+1)-lon(j-1))*pi*r*cosd(lat(i))/360.)*((lat(i+1)-lat(i-1))*pi*r/360.) !计算单元格点面积
            elseif (j==1 .and. i>1 .and. i<imt) then
                dxyj(i,j)=((lon(j+1)-lon(j))*pi*r*cosd(lat(i))/180.)*((lat(i+1)-lat(i-1))*pi*r/360.)
            elseif (j==1 .and. i==imt) then
                dxyj(i,j)=((lon(j+1)-lon(j))*pi*r*cosd(lat(i))/180.)*((lat(i)-lat(i-1))*pi*r/180.d0)
            elseif (j==1 .and. i==1) then
                dxyj(i,j)=((lon(j+1)-lon(j))*pi*r*cosd(lat(i))/180.)*((lat(i+1)-lat(i))*pi*r/180.)
            elseif (j>1 .and. j<jmt .and. i==1) then
                dxyj(i,j)=((lon(j+1)-lon(j-1))*pi*r*cosd(lat(i))/360.)*((lat(i+1)-lat(i))*pi*r/180.)
            elseif (j==jmt .and. i==1) then
                dxyj(i,j)=((lon(j)-lon(j-1))*pi*r*cosd(lat(i))/180.)*((lat(i+1)-lat(i))*pi*r/180.)
            elseif (j>1 .and. j<jmt .and. i==imt) then
                dxyj(i,j)=((lon(j+1)-lon(j-1))*pi*r*cosd(lat(i))/360.)*((lat(i)-lat(i-1))*pi*r/180.)
            elseif (j==jmt .and. i==imt) then
                dxyj(i,j)=((lon(j)-lon(j-1))*pi*r*cosd(lat(i))/180.)*((lat(i)-lat(i-1))*pi*r/180.)
            elseif (j==jmt .and. i>1 .and. i<imt) then
                dxyj(i,j)=((lon(j)-lon(j-1))*pi*r*cosd(lat(i))/180.)*((lat(i+1)-lat(i-1))*pi*r/360.)
            end if
        end do
    end do


    do k=1,kmt
        do j=1,jmt
            do i=1,imt
                if (t(j,i,k)<9999.0) then
                    area(k)=area(k)+dxyj(i,j) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 定义变量初始值问题
                end if
             end do
         end do
    end do

    do k=kmt,1,-1
        if (area(k)>0.0d0)then
            kk=k
            kc=k
        exit
        end if
    end do

    pb0=ceiling(sw_pres(zb(kc),30.0d0))
    nmt=ceiling(pb0/dpa)+1
    allocate(dpb(nmt-1))
    dpb(1:nmt-2)=dpa*10000.0d0
    dpb(nmt-1)=(pb0-dpa*(nmt-2))*10000.0d0
    allocate(pre_rho(nmt))
    do n=1,nmt
        if (n==1) then
            pre_rho(n)=pb0
        else
            pre_rho(n)=pre_rho(n-1)-dpb(n-1)/10000.0d0
        end if
    end do




    !print *, pb0
    !print *, dpb
    !print *,pre_rho
    !stop

    !!!!计算现场势能
    if (t_potential==0) then
        do i=1,imt
            do j=1,jmt
                do k=1,kmt
                    if (t(j,i,k)<9999.0d0) then
                    pres=sw_pres(lev(k),lat(i))
                    den(k,i,j)=sw_dens(s(j,i,k),t(j,i,k),pres)
                    gpe(k,i,j)=den(k,i,j)*dxyj(i,j)*dz(k)*g*zzc(k)+ave_dmass*g*zzc(k)
                    mass(k,i,j)=den(k,i,j)*dxyj(i,j)*dz(k)
                    volume(k,i,j)=dxyj(i,j)*dz(k)
                    else
                    gpe(k,i,j)=9999.9999d0
                    mass(k,i,j)=9999.9999d0
                    volume(k,i,j)=9999.9999d0
                    den(k,i,j)=9999.9999d0
                    end if
                end do
            end do
        end do
    elseif (t_potential==1) then
        do i=1,imt
            do j=1,jmt
                do k=1,kmt
                    if (t(j,i,k)<9999.0d0) then
                    pres=sw_pres(lev(k),lat(i))
                    t_p=sw_temp(s(j,i,k),t(j,i,k),pres,0.0d0)
                    den(k,i,j)=sw_dens(s(j,i,k),t_p,pres)
                    gpe(k,i,j)=den(k,i,j)*dxyj(i,j)*dz(k)*g*zzc(k)+ave_dmass*g*zzc(k)
                    mass(k,i,j)=den(k,i,j)*dxyj(i,j)*dz(k)
                    volume(k,i,j)=dxyj(i,j)*dz(k)
                    else
                    gpe(k,i,j)=9999.9999d0
                    mass(k,i,j)=9999.9999d0
                    volume(k,i,j)=9999.9999d0
                    den(k,i,j)=9999.9999d0
                    end if
                end do
            end do
        end do
    end if

    !!!!!!!!!!!!计算参考面势能
    num1=0
    jd=1
    allocate(gpe1(ijkm),mass1(ijkm),s1(ijkm),t1(ijkm),index_o(ijkm,3),pden1(ijkm),pden1_r(ijkm),pden1_o(ijkm),pden1_c(ijkm),dh(ijkm),dp(ijkm),index1(ijkm),index1_o(ijkm),index_c(ijkm),volume_ref1(ijkm))
    do i=1,imt
        do j=1,jmt
            do k=1,kmt
                if (t(j,i,k)<9999.0d0) then
                    num1=num1+1
                    gpe1(num1)=gpe(k,i,j)
                    mass1(num1)=mass(k,i,j)
                    s1(num1)=s(j,i,k)
                    t1(num1)=t(j,i,k)
                    index_o(num1,1)=i
                    index_o(num1,2)=j
                    index_o(num1,3)=k
                end if
            end do
        end do
    end do

    dist=dz(kk)
    jcheck=1
    if (t_potential==0) then

    do n=1,nmt
        print *, n
            do i=1,num1
                if (t1(i)<9999.0d0) then
                    pres=sw_pres(lev(index_o(i,3)),lat(index_o(i,1)))
                    pden1(i)=sw_pden(s1(i),t1(i),pres,pre_rho(n))
                else
                    pden1(i)=9999.9999d0
                end if
            end do

        call indexx(num1,pden1,index1_o)
        pden1_o=pden1
        call sort(num1,pden1_o)
        do i=1,num1                                                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 转换为从大到小
            pden1_r(i)=pden1_o(num1+1-i)
            index1(i)=index1_o(num1+1-i)
        end do

        do j=jd,num1
            dh(j)=mass1(index1(j))/pden1_r(j)/area(kk)
            dp(j)=pden1_r(j)*g*dh(j)
            
            !if (pden1_r(j)>9999.  .or.  pden1_r(j)<1000.) then ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         检错
            !print *, 'error on pden1_r is '
            !print *, pden1_r(j)
            !stop
            !end if

            !if (sum(dh(1:j))>dist) then
             !   if (kk>1) then
             !       kk=kk-1
             !   else
             !       kk=1
             !   end if
              !  dist=dist+dz(kk)
            !end if

            if (sum(dh(jcheck:j))>=dist) then
                if (kk>1) then
                    kk=kk-1
                else
                    kk=1
                end if
                dist=dz(kk)
                jcheck=j+1
            end if


            !test=sum(dh(1:j))
            !print *, n, test, kk, dist, j, num1 , mass1(index1(j)), pden1_r(j)      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   检错

            !do k=1,kmt
             !   dist(k)=zb(kc)-sum(dh(1:j))-zt(k)
              !  if (dist(k)<0) then
               !     mask(k)=.false.
                !else
                 !   mask(k)=.true.
                !end if
            !end do
            !it=minloc(dist,mask)
            !if (it(1)<kk) then
             !   kk=it(1)
           ! end if

            if (n<nmt) then
                if (sum(dp(jd:j))>=dpb(n)) then
                    x=j
                    exit
                end if
            else
                x=num1
                exit
            end if
        end do
        pden1_c(jd:x)=pden1_r(jd:x)
        index_c(jd:x)=index1(jd:x)

        s1(index1(jd:x))=9999.9999d0
        t1(index1(jd:x))=9999.9999d0
        jd=x+1
    end do

    else if (t_potential==1) then

        do n=1,nmt
        print *, n
            do i=1,num1
                if (t1(i)<9999.0d0) then
                    pres=sw_pres(lev(index_o(i,3)),lat(index_o(i,1)))
                    t_p=sw_temp(s1(i),t1(i),pres,0.d0)
                    pden1(i)=sw_pden(s1(i),t_p,pres,pre_rho(n))

                else
                    pden1(i)=9999.9999d0
                end if
            end do
            
        call indexx(num1,pden1,index1_o)
        pden1_o=pden1
        call sort(num1,pden1_o)
        do i=1,num1                                                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 转换为从大到小
            pden1_r(i)=pden1_o(num1+1-i)
            index1(i)=index1_o(num1+1-i)
        end do

        do j=jd,num1
            dh(j)=mass1(index1(j))/pden1_r(j)/area(kk)
            dp(j)=pden1_r(j)*g*dh(j)
            volume_ref1(j)=dh(j)*area(kk)
            
            !if (pden1_r(j)>9999.  .or.  pden1_r(j)<1000.) then ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         检错
            !print *, 'error on pden1_r is '
            !print *, pden1_r(j)
            !stop
            !end if

           ! if (sum(dh(1:j))>dist) then
             !   if (kk>1) then
                !    kk=kk-1
             !   else
             !       kk=1
             !   end if
             !   dist=dist+dz(kk)
            !end if



            if (sum(dh(jcheck:j))>=dist) then
                if (kk>1) then
                    kk=kk-1
                else
                    kk=1
                end if
                dist=dz(kk)
                jcheck=j+1
            end if


            !test=sum(dh(1:j))
            !print *, n, test, kk, dist, j, num1 , mass1(index1(j)), pden1_r(j)      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   检错

            !do k=1,kmt
             !   dist(k)=zb(kc)-sum(dh(1:j))-zt(k)
              !  if (dist(k)<0) then
               !     mask(k)=.false.
                !else
                 !   mask(k)=.true.
                !end if
            !end do
            !it=minloc(dist,mask)
            !if (it(1)<kk) then
             !   kk=it(1)
           ! end if

            if (n<nmt) then
                if (sum(dp(jd:j))>=dpb(n)) then
                    x=j
                    exit
                end if
            else
                x=num1
                exit
            end if
        end do
        pden1_c(jd:x)=pden1_r(jd:x)
        index_c(jd:x)=index1(jd:x)

        s1(index1(jd:x))=9999.9999d0
        t1(index1(jd:x))=9999.9999d0
        jd=x+1
    end do
    end if

    allocate(h(num1))
    h(1)=0.5d0*dh(1)
    do j=2,num1
        h(j)=h(j-1)+0.5d0*(dh(j)+dh(j-1))
    end do

        do i=1,num1
        do j=1,num1
            if (index_c(j)==i) then
                gper(index_o(i,3),index_o(i,1),index_o(i,2))=mass1(index_c(j))*g*h(j)+ave_dmass*g*h(j)
                zr(index_o(i,3),index_o(i,1),index_o(i,2))=h(j)
                denref(index_o(i,3),index_o(i,1),index_o(i,2))=pden1_c(j)
                volume_ref(index_o(i,3),index_o(i,1),index_o(i,2))=volume_ref1(j)
            end if
        end do
    end do
    do i=1,imt
        do j=1,jmt
            do k=1,kmt
                if (t(j,i,k)>9999.0d0) then
                    gper(k,i,j)=9999.9999d0
                    zr(k,i,j)=9999.9999d0
                    denref(k,i,j)=9999.9999d0
                    volume_ref(k,i,j)=9999.9999d0
                end if
            end do
        end do
    end do
    
    do i=1,imt
        do j=1,jmt
            do k=1,kmt
                if (t(j,i,k)>9999.0d0) then
                        agpe(k,i,j)=9999.9999d0
                    else
                        agpe(k,i,j)=gpe(k,i,j)-gper(k,i,j)
                end if
            end do
        end do
    end do

    sum_agpe=0.d0
    sum_volume=0.d0
    do i=1,imt
        do j=1,jmt
            do k=1,kmt
                if (t(j,i,k)<9999.d0) then
                    sum_agpe=sum_agpe+agpe(k,i,j)
                    sum_volume=sum_volume+volume(k,i,j)
                    else
                    sum_agpe=sum_agpe
                    sum_volume=sum_volume
                end if
            end do
        end do
    end do

    ave_agpe=sum_agpe/sum_volume

    write (*,*) ' agpe_ave= '
    write (*,'(ES16.7)'), ave_agpe
    write (*,*) ' hmax= '
    write (*,'(ES16.7)'), h(num1)

    open (unit=25,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\agpe_'//savename,status='new',action='write')
    do k=1,kmt
        do i=1,imt
            write(25,'(360ES16.7)') agpe(k,i,:)
        end do
    end do
    close (25)   

    open (unit=25,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\denref_'//savename,status='new',action='write')
    do k=1,kmt
        do i=1,imt
            write(25,'(360ES16.7)') denref(k,i,:)
        end do
    end do
    close (25)  
    
    open (unit=25,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\gper_'//savename,status='new',action='write')
    do k=1,kmt
        do i=1,imt
            write(25,'(360ES16.7)') gper(k,i,:)
        end do
    end do
    close (25)  
    
    open (unit=25,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\zr_'//savename,status='new',action='write')
    do k=1,kmt
        do i=1,imt
            write(25,'(360ES16.7)') zr(k,i,:)
        end do
    end do
    close (25)  
    
    open (unit=26,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\volume_'//savename,status='new',action='write')
    do k=1,kmt
        do i=1,imt
            write(26,'(360ES16.7)') volume(k,i,:)
        end do
    end do
    close (26)   

    open (unit=26,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\volume_ref_'//savename,status='new',action='write')
    do k=1,kmt
        do i=1,imt
            write(26,'(360ES16.7)') volume_ref(k,i,:)
        end do
    end do
    close (26)  
    
    open (unit=27,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\gpe_'//savename,status='new',action='write')
    do k=1,kmt
        do i=1,imt
            write(27,'(360ES16.7)') gpe(k,i,:)
        end do
    end do
    close (27)   

    open (unit=28,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\dxyj_'//savename,status='new',action='write')
    do i=1,imt
        write(28,'(360ES16.7)') dxyj(i,:)
    end do
    close (28)   

    open (unit=29,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\dz_'//savename,status='new',action='write')
    do k=1,kmt
        write(29,'(f11.4)') dz(k)
    end do
    close (29)   

    open (unit=29,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\zzc_'//savename,status='new',action='write')
    do k=1,kmt
        write(29,'(f11.4)') zzc(k)
    end do
    close (29)   
    
    open (unit=44,file='G:\LIURAN_SCIENCE\SODA_pacific\soda_agpe\den_'//savename,status='new',action='write')
    do k=1,kmt
        do i=1,imt
            write(44,'(360ES16.7)') den(k,i,:)
        end do
    end do
    close (44)  

    end subroutine AGPE_CALCULATE
