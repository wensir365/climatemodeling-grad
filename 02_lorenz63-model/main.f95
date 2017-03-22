! Xinyu Wen, Peking Univ, Mar/20/2017

program main

	implicit none

	integer, parameter      :: numvar=3		! number of state variables
	real, dimension(numvar) :: v			! state vector
	integer :: i							! loop var

	! init									! need high-order time scheme if far from attractors
	v(1) = 5.0
	v(2) = 5.0
	v(3) = 5.0
	print *, v

	! main loop
	do i = 1, 10000
		call heun2(v)
		print *, v
	end do

contains

	subroutine heun2(x)
		implicit none
		real, dimension(numvar), intent(inout) :: x
		real, dimension(numvar) :: xest, tend1, tend2
		real, parameter :: dt=0.01

		call computetend(x,tend1)
		xest	= x + dt*tend1

		call computetend(xest,tend2)
		x		= x + 0.5*dt*(tend1+tend2)
	end subroutine

	subroutine computetend(x,tend)
		implicit none
		real, dimension(numvar), intent(in)  :: x
		real, dimension(numvar), intent(out) :: tend
		! Lorenz 63 Model
		tend(1) = 10.0*(x(2)-x(1))			! dx/dt = 10(y-x)
		tend(2) = x(1)*(28.0-x(3))-x(2)		! dy/dt = x(28-z)-y
		tend(3) = x(1)*x(2)-8.0/3.0*x(3)	! dz/dt = xy-8z/3
	end subroutine

end program
