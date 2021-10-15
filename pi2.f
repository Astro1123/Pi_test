      program picalc
          implicit none
          integer count,i
          double precision::pi,x,x2,Newton,ll,sl,regular_polygon
          double precision::Accelerated,GrL,Sharp,EM,GaL
          integer::ti,tf,tr
          
  100     write (*,*) 'Loop'
          read(*,*) count
          if(count<0) then
            goto 100
          end if

          do i=1,3,1
          write (*,*) ''
          end do

          pi=acos(-1d0)
          write (*,200) 'pi = ',pi
          write (*,*) ''

          call system_clock(ti)

          x=3d0
          do i=0,count-1,1
            x=Newton(x)
          end do

          call system_clock(tf,tr)

          write (*,200) 'Newton method : ',x,' ('
     &       ,(tf-ti)/dble(tr),' sec)'
          write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
          write (*,*) ''

          call system_clock(ti)

          ll=4d0
          do i=1,count,1
            ll=regular_polygon(ll,i)
          end do
          sl=ll/sqrt(1d0+(ll/2d0**(count+2))**2)

          call system_clock(tf,tr)

          write (*,400) 'Polygon : ',sl,ll,' ('
     &       ,(tf-ti)/dble(tr),' sec)'
          write (*,'(a25,es10.3,a6)') '( ',2d0**(count+2),'-gon )'
          write (*,'(a30,es15.8,a6)') '( error (Circumscribe) : '
     &       ,abs(ll-pi),' )'
          write (*,'(a30,es15.8,a6)') '( error (Inscribe) : '
     &       ,abs(sl-pi),' )'
          write (*,*) ''

          if (count>8000) then
            goto 150
          end if
          call system_clock(ti)

          x=Accelerated(count)

          call system_clock(tf,tr)

          write (*,200) 'Acceleration method : ',x,' ('
     &       ,(tf-ti)/dble(tr),' sec)'
          write (*,'(A50)') '(Takahiro Takebe (1664 - 1739))'
          write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
          write (*,*) ''

  150     call system_clock(ti)

          x=Sharp(count)

          call system_clock(tf,tr)

          write (*,200) 'Sharp method : ',x,' ('
     &       ,(tf-ti)/dble(tr),' sec)'
          write (*,'(A30)') '(1699)'
          write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
          write (*,*) ''

          call system_clock(ti)

          x=EM(count)

          call system_clock(tf,tr)

          write (*,200) 'Euler-Machin method : ',x,' ('
     &       ,(tf-ti)/dble(tr),' sec)'
          write (*,'(A30)') '(1706)'
          write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
          write (*,*) ''

          call system_clock(ti)

          x=GaL(count)

          call system_clock(tf,tr)

          write (*,200) 'Gauss-Legendre algorithm : ',x,' ('
     &       ,(tf-ti)/dble(tr),' sec)'
          write (*,'(A30)') '(1976)'
          write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
          write (*,*) ''

          call system_clock(ti)

          x=GrL(count)

          call system_clock(tf,tr)

          write (*,200) 'Gregory-Leibniz method : ',x,' ('
     &       ,(tf-ti)/dble(tr),' sec)'
          write (*,'(A30)') '(1671)'
          write (*,'(a30,es15.8,a6)') '( error : ',abs(x-pi),' )'
          write (*,*) ''


          write (*,*) ''
          write (*,300) 'Loop : ',count

          do i=1,3,1
          write (*,*) ''
          end do
          stop

  200 format (a30,f33.31,a,f7.3,a)
  300 format (a30,i10)
  400 format (a30,f33.31,' < pi < ',f33.31,a,f7.3,a)
      end program

      function Newton(x)
          implicit none
          double precision::Newton,x
          Newton=x-sin(x)/cos(x)
      end function

      function regular_polygon(ll,n)
          implicit none
          integer n
          double precision::regular_polygon,ll,sl
	  sl=ll/sqrt(1d0+(ll/2d0**(n+1))**2)
	  regular_polygon=2d0/(1d0/sl+1d0/ll)
      end function

      function Accelerated(count)
          implicit none
          integer count,i,j
          double precision::Accelerated,fact,x
          x=0d0
          do i=0,count-1,1
	      fact=1
              do j=i,1,-1
		  fact=fact*(j*j/(2d0+2d0*j)/(1d0+2d0*j))
              end do
              x=x+fact
	  end do
	  Accelerated=sqrt(x*9d0)
      end function

      function Sharp(count)
          implicit none
          integer count,i
          double precision::Sharp,x
          x=0d0
          do i=0,count-1,1
              x=x+(-1d0/3d0)**i/(2d0*i+1d0)
	  end do
	  Sharp=x*sqrt(3d0)*2
      end function

      function GrL(count)
          implicit none
          integer count,i
          double precision::GrL,x
          x=0d0
          do i=0,count-1,1
              x=x+(-1d0)**i/(2d0*i+1d0)
	  end do
	  GrL=x*4d0
      end function

      function EM(count)
          implicit none
          integer count,i
          double precision::EM,x,y
          x=0d0
          y=0d0
          do i=0,count-1,1    
	      x=x+((1d0/5d0)**(1d0+2d0*i))/(2d0*i+1d0)*(-1d0)**i
	      y=y+((1d0/239d0)**(1d0+2d0*i))/(2d0*i+1d0)*(-1d0)**i
	  end do
	  EM=(x*4d0-y)*4d0
      end function


      function GaL(count)
          implicit none
          integer count,i
          double precision::GaL,a,b,t,p,a_,b_,t_,p_
          a=1d0
          b=1d0/sqrt(2d0)
          t=1d0/4d0
          p=1d0
          a_=a
          b_=b
          t_=t
          p_=p
          do i=0,count-1,1
              a=(a_+b_)/2d0
              b=sqrt(a_*b_)
              t=t_-p_*(a_-a)**2d0
              p=2d0*p_
              a_=a
              b_=b
              t_=t
              p_=p
	  end do
	  GaL=(a_+b_)**2/(4d0*t_)
      end function
