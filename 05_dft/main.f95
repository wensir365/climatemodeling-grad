program main

	integer, parameter 	:: N = 144
	real, dimension(N)	:: z = (/ &
	   5745.00,  5730.00,  5712.00,  5692.00,  5671.00,  5647.00,  5620.00,  5586.00,&
	   5545.00,  5497.00,  5448.00,  5402.00,  5366.00,  5344.00,  5338.00,  5348.00,&
	   5367.00,  5391.00,  5413.00,  5430.00,  5439.00,  5443.00,  5442.00,  5438.00,&
	   5431.00,  5422.00,  5412.00,  5405.00,  5405.00,  5418.00,  5443.00,  5478.00,&
	   5516.00,  5550.00,  5574.00,  5585.00,  5586.00,  5579.00,  5568.00,  5553.00,&
	   5534.00,  5509.00,  5480.00,  5449.00,  5420.00,  5398.00,  5384.00,  5375.00,&
	   5368.00,  5360.00,  5348.00,  5335.00,  5323.00,  5315.00,  5309.00,  5305.00,&
	   5301.00,  5294.00,  5287.00,  5280.00,  5275.00,  5273.00,  5274.00,  5276.00,&
	   5279.00,  5284.00,  5290.00,  5299.00,  5309.00,  5318.00,  5323.00,  5326.00,&
	   5326.00,  5327.00,  5331.00,  5340.00,  5353.00,  5366.00,  5377.00,  5382.00,&
	   5379.00,  5369.00,  5352.00,  5330.00,  5306.00,  5282.00,  5262.00,  5249.00,&
	   5247.00,  5255.00,  5273.00,  5295.00,  5317.00,  5339.00,  5358.00,  5378.00,&
	   5400.00,  5424.00,  5448.00,  5469.00,  5482.00,  5485.00,  5477.00,  5461.00,&
	   5444.00,  5426.00,  5410.00,  5392.00,  5367.00,  5331.00,  5285.00,  5234.00,&
	   5187.00,  5153.00,  5137.00,  5139.00,  5154.00,  5175.00,  5194.00,  5205.00,&
	   5208.00,  5202.00,  5190.00,  5178.00,  5171.00,  5177.00,  5201.00,  5243.00,&
	   5303.00,  5375.00,  5453.00,  5530.00,  5601.00,  5661.00,  5708.00,  5741.00,&
	   5761.00,  5771.00,  5772.00,  5770.00,  5767.00,  5764.00,  5761.00,  5755.00 /), z0,z1,z2,z3,z5,z10,z20,z30

	complex, dimension(N) :: x,y
	complex, dimension(-N/2:N/2-1)	:: c, c_tr

	real :: dx
	integer :: i
	
	x 		= cmplx(z)
	call dft_1d_cplx(x,N,c)

	!--- Truncation
	call trunc(N,c,0,z0)
	call trunc(N,c,1,z1)
	call trunc(N,c,2,z2)
	call trunc(N,c,3,z3)
	call trunc(N,c,5,z5)
	call trunc(N,c,10,z10)
	call trunc(N,c,20,z20)
	call trunc(N,c,30,z30)

	!--- output
	dx	= 360.0/real(N)
	do i = 1, N
		print "(10f10.2)", (real(i)*dx-dx/2.0), z(i), z0(i),z1(i),z2(i),z3(i),z5(i),z10(i),z20(i),z30(i)
	end do

contains

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
		integer, intent(in)							:: N
		complex, dimension(1:N), intent(in)			:: x
		complex, dimension(-N/2:N/2-1), intent(out) :: coef

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
		integer, intent(in)							:: N
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

end program
