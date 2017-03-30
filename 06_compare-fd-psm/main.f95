program main

	real, parameter		:: pi = 3.1415926
	real, parameter		:: pi2 = (pi*2.0)
	complex, parameter	:: jj = (0.0,1.0)

	integer, parameter 	:: N = 64
	real, dimension(N)	:: x,yana,yfd,ypsm
	integer, dimension(N)	:: pl,pr

	real :: dt=0.1, dx=pi2/real(N), c=0.1
	integer :: i

	real 	:: onecycle, numcycle

	!print *, "Courant Number =", (dt*c)/dx

	onecycle	= pi2/(dt*c)
	numcycle	= 10

	! init
	do i	= 1, N
		x(i)	= (0.5+i-1)*dx
		pl(i)	= i-1
		pr(i)	= i+1
	end do
	pl(1)	= N
	pr(N)	= 1

	call sety(0,x,yana)
	call sety(0,x,yfd)
	call sety(0,x,ypsm)

	! loop
	do i	= 1, int(onecycle*numcycle)
		!Analytic
		call sety(i,x,yana)

		!FD
		call heun2_fd(yfd)

		!pesudo-SM
		call heun2_psm(ypsm)

		!output
	end do

	do i = 1, N
		print "(4f12.8)", x(i), yana(i),yfd(i),ypsm(i)
	end do

contains

	subroutine sety(step,x,y)
		implicit none
		integer, intent(in) 		:: step
		real, dimension(N), intent(in)	:: x
		real, dimension(N), intent(out) :: y
		integer :: i
		real	:: xmid
		xmid	= pi+c*dt*step
		if (xmid.gt.pi2) then
			xmid	= xmid-(pi2*int(xmid/pi2))
		end if
		do i	= 1, N
			y(i)	= gauss(x(i),0.5,xmid)
		end do
	end subroutine

	subroutine psm(y,tend)
		implicit none
		real, dimension(N), intent(in)	:: y
		real, dimension(N), intent(out)	:: tend
		complex, dimension(N)		:: y_cplx
		complex, dimension(-N/2:N/2-1)	:: ck
		integer :: k
		y_cplx 	= cmplx(y)
		call dft_1d_cplx(y_cplx,N,ck)
		do k = -N/2, N/2-1
			!ck(k) = -1.0*c*ck(k)*(jj*pi2*k)/N
			ck(k) = -1.0*c*ck(k)*(jj*k)
		end do
		call idft_1d_cplx(ck,N,y_cplx)
		tend	= real(y_cplx)
	end subroutine

	subroutine heun2_psm(y)
		implicit none
		real, dimension(N), intent(inout) :: y
		real, dimension(N) :: yest, tend1, tend2
		call psm(y,tend1)
		yest	= y + dt*tend1
		call psm(yest,tend2)
		y	= y + 0.5*dt*(tend1+tend2)
	end subroutine

	subroutine fd(y,tend)
		implicit none
		real, dimension(N), intent(in)  :: y
		real, dimension(N), intent(out) :: tend
		integer :: i
		do i = 1, N
			tend(i) = -1.0*c*(y(pr(i))-y(pl(i)))/(2.0*dx)
		end do
	end subroutine

	subroutine heun2_fd(y)
		implicit none
		real, dimension(N), intent(inout) :: y
		real, dimension(N) :: yest, tend1, tend2
		call fd(y,tend1)
		yest 	= y + dt*tend1
		call fd(yest,tend2)
		y	= y + 0.5*dt*(tend1+tend2)
	end subroutine

	function gauss(x,sigma,mean)
		implicit none
		real, intent(in) :: x, sigma, mean
		real :: gauss
		gauss = 1.0/sqrt(2.0*pi*sigma**2) * exp((x-mean)**2/(-2.0*sigma**2))
	end function

	subroutine trunc(N,coef,kmax,y)
		implicit none
		integer, intent(in) :: N
		complex, dimension(-N/2:N/2-1), intent(in) :: coef
		integer, intent(in) :: kmax
		real, dimension(N), intent(out) :: y
		complex, dimension(-N/2:N/2-1) :: coeftmp, ytmp

		coeftmp					= coef
		coeftmp(kmax+1:N/2-1)	= (0,0)
		coeftmp(-N/2:-1-kmax)	= (0,0)
		call idft_1d_cplx(coeftmp,N,ytmp)
		y	= real(ytmp)
	end subroutine

	subroutine dft_1d_cplx(x,N,coef)
		implicit none
		integer, intent(in)				:: N
		complex, dimension(1:N), intent(in)		:: x
		complex, dimension(-N/2:N/2-1), intent(out) 	:: coef

		complex, parameter	:: jj = (0,1)
		real, parameter		:: Pi = 3.1415926
		integer :: i, k

		do k = -N/2, N/2-1
			coef(k) = (0,0)
			do i = 1, N
				coef(k) = coef(k) + x(i)*exp(-2.0*jj*Pi*k*i/N)
			end do
			coef(k) = coef(k)/N
		end do
	end subroutine

	subroutine idft_1d_cplx(coef,N,x)
		implicit none
		integer, intent(in)				:: N
		complex, dimension(-N/2:N/2-1), intent(in) 	:: coef
		complex, dimension(1:N), intent(out)		:: x

		complex, parameter	:: jj = (0,1)
		real, parameter		:: Pi = 3.1415926
		integer :: i, k

		do i = 1, N 
			x(i) = (0,0)
			do k = -N/2, N/2-1
				x(i) = x(i) + coef(k)*exp(2.0*jj*Pi*k*i/N)
			end do
		end do
	end subroutine

	subroutine dft_1d_real(x,N,a,b)
		!-- interface --
		integer, intent(in) :: N
		real, dimension(N), intent(in) :: x
		real, dimension(0:N/2), intent(out) :: a, b
		!-- local --
		real, parameter		:: Pi = 3.1415926
		integer :: i, k

		do k = 0, N/2
			a(k) = 0.0
			b(k) = 0.0
			do i = 1, N		! Sigma
				a(k) = a(k) + x(i)*cos(2.0*Pi*k*i/N)
				b(k) = b(k) + x(i)*sin(2.0*Pi*k*i/N)
			end do
			a(k) = 2.0*a(k)/N
			b(k) = 2.0*b(k)/N
		end do
		a(0) 	= a(0)/2.0
		a(N/2) 	= a(N/2)/2.0
	end subroutine

	subroutine idft_1d_real(a,b,N,x)
		!-- interface --
		integer, intent(in) :: N
		real, dimension(0:N/2), intent(in) :: a, b
		real, dimension(1:N), intent(out) :: x
		!-- local --
		real, parameter		:: Pi = 3.1415926
		integer :: i, k

		do i = 1, N
			x(i) = 0.0
			do k = 0, N/2
				x(i) = x(i) + a(k)*cos(2.0*Pi*k*i/N) + b(k)*sin(2.0*Pi*k*i/N)
			end do
		end do 
	end subroutine

end program
