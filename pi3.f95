program picalc
    implicit none
    integer count,i
    real*16::pi,x,x2,Newton,ll,sl,regular_polygon
    real*16::Accelerated,GrL,Sharp,EM,GaL
    integer::ti,tf,tr

    do
        write (*,*) 'Loop'
        read(*,*) count
        if(count >= 0) then
            exit
        end if
    end do

    do i=1,3,1
        write (*,*) ''
    end do

    pi=acos(-1q0)
    write (*,'(a30,f33.31,a,f7.3,a)') 'pi = ',pi
    write (*,*) ''

    call system_clock(ti)

    x=3q0
    do i=0,count-1,1
        x=Newton(x)
    end do

    call system_clock(tf,tr)

    write (*,'(a30,f33.31,a,f7.3,a)') 'Newton method : ',x,' (',(tf-ti)/dble(tr),' sec)'
    write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
    write (*,*) ''

    call system_clock(ti)

    ll=4q0
    do i=1,count,1
        ll=regular_polygon(ll,i)
    end do
    sl=ll/sqrt(1.0q0+(ll/2.0q0**(count+2))**2)

    call system_clock(tf,tr)

    write (*,'(a30,f33.31,a,f33.31,a,f7.3,a)') 'Polygon : ',sl,' < pi < ',ll,' (',(tf-ti)/dble(tr),' sec)'
    write (*,'(a25,es10.3,a6)') '( ',2d0**(count+2),'-gon )'
    write (*,'(a30,es15.8,a6)') '( error (Circumscribe) : ',abs(ll-pi),' )'
    write (*,'(a30,es15.8,a6)') '( error (Inscribe) : ',abs(sl-pi),' )'
    write (*,*) ''
          
    if (count <= 8000) then
       call system_clock(ti)

        x=Accelerated(count)

        call system_clock(tf,tr)

        write (*,'(a30,f33.31,a,f7.3,a)') 'Acceleration method : ',x,' (',(tf-ti)/dble(tr),' sec)'
        write (*,'(A50)') '(Takahiro Takebe (1664 - 1739))'
        write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
        write (*,*) ''

    end if

    call system_clock(ti)

    x=Sharp(count)

    call system_clock(tf,tr)

    write (*,'(a30,f33.31,a,f7.3,a)') 'Sharp method : ',x,' (',(tf-ti)/dble(tr),' sec)'
    write (*,'(A30)') '(1699)'
    write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
    write (*,*) ''

    call system_clock(ti)

    x=EM(count)

    call system_clock(tf,tr)

    write (*,'(a30,f33.31,a,f7.3,a)') 'Euler-Machin method : ',x,' (',(tf-ti)/dble(tr),' sec)'
    write (*,'(A30)') '(1706)'
    write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
    write (*,*) ''

    call system_clock(ti)

    x=GaL(count)

    call system_clock(tf,tr)

    write (*,'(a30,f33.31,a,f7.3,a)') 'Gauss-Legendre algorithm : ',x,' (',(tf-ti)/dble(tr),' sec)'
    write (*,'(A30)') '(1976)'
    write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
    write (*,*) ''

    call system_clock(ti)

    x=GrL(count)

    call system_clock(tf,tr)

    write (*,'(a30,f33.31,a,f7.3,a)') 'Gregory-Leibniz method : ',x,' (',(tf-ti)/dble(tr),' sec)'
    write (*,'(A30)') '(1671)'
    write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
    write (*,*) ''


    write (*,*) ''
    write (*,'(a30,i10)') 'Loop : ',count

    do i=1,3,1
        write (*,*) ''
    end do
    stop

end program

function Newton(x)
    implicit none
    real*16::Newton,x
    Newton=x-sin(x)/cos(x)
end function

function regular_polygon(ll,n)
    implicit none
    integer n
    real*16::regular_polygon,ll,sl
    sl=ll/sqrt(1.0q0+(ll/2.0q0**(n+1))**2)
    regular_polygon=2.0q0/(1.0q0/sl+1.0q0/ll)
end function

function Accelerated(count)
    implicit none
    integer count,i,j
    real*16::Accelerated,fact,x
    x=0q0
    do i=0,count-1,1
        fact=1q0
        do j=i,1,-1
            fact=fact*(j*j/(2q0+2q0*j)/(1q0+2q0*j))
        end do
        x=x+fact
    end do
    Accelerated=sqrt(x*9q0)
end function

function Sharp(count)
    implicit none
    integer count,i
    real*16::Sharp,x
    x=0q0
    do i=0,count-1,1
        x=x+(-1q0/3q0)**i/(2q0*i+1q0)
    end do
    Sharp=x*sqrt(3q0)*2q0
end function

function GrL(count)
    implicit none
    integer count,i
    real*16::GrL,x
    x=0q0
    do i=0,count-1,1
        x=x+(-1q0)**i/(2q0*i+1q0)
    end do
    GrL=x*4q0
end function

function EM(count)
    implicit none
    integer count,i
    real*16::EM,x,y
    x=0q0
    y=0q0
    do i=0,count-1,1    
        x=x+((1q0/5q0)**(1q0+2q0*i))/(2q0*i+1q0)*(-1q0)**i
        y=y+((1q0/239q0)**(1q0+2q0*i))/(2q0*i+1q0)*(-1q0)**i
    end do
    EM=(x*4q0-y)*4q0
end function


function GaL(count)
    implicit none
    integer count,i
    real*16::GaL,a,b,t,p,a_,b_,t_,p_
    a=1q0
    b=1q0/sqrt(2q0)
    t=1q0/4q0
    p=1q0
    a_=a
    b_=b
    t_=t
    p_=p
    do i=0,count-1,1
        a=(a_+b_)/2.0
        b=sqrt(a_*b_)
        t=t_-p_*(a_-a)**2q0
        p=2q0*p_
        a_=a
        b_=b
        t_=t
        p_=p
    end do
    GaL=(a_+b_)**2/(4q0*t_)
end function
