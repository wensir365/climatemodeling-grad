! Xinyu Wen, Peking Univ, Mar/20/2017

program main

	implicit none

	integer, parameter :: numvar=2		! number of state variables in dynamical system
	integer, parameter :: numlines=1000	! number of testing trajactories

	real, dimension(numvar) :: v		! state vector
	integer :: i, lines					! for loop


	do lines = 1, numlines
		call random_number(v)			! random starting point at [0,1]
		v	= (v-0.5)*10				! adjust starting point from [0,1] to [-5,5]
		print *, v

		do i = 1, 100
			call heun2(v)				! old v in, new v out
			print *, v
		end do
		print *, " "					! print a blank line to separate trajactories
	end do

contains

	subroutine heun2(x)
		implicit none
		real, dimension(numvar), intent(inout) :: x
		real, dimension(numvar) :: xest, tend1, tend2
		real, parameter :: dt = 0.05

		call computetend(x,tend1)
		xest 	= x + dt*tend1

		call computetend(xest,tend2)
		x	= x + 0.5*dt*(tend1+tend2)
	end subroutine

	subroutine computetend(x,tend)
		implicit none
		real, dimension(numvar), intent(in)  :: x
		real, dimension(numvar), intent(out) :: tend

		tend(1) = 1.0-x(1)*x(1)			! dx/dt = 1-x^2
		tend(2) = 1.0-x(1)*x(2)			! dy/dt = 1-xy
	end subroutine

end program
