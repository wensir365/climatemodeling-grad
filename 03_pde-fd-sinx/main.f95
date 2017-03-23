! Xinyu Wen, Peking Univ, Mar/23/2017
program main

	implicit none

	real, parameter 	:: PI = 3.1415926
	integer, parameter 	:: N = 16			! number of grids

	real, dimension(N) 	:: x,sinx,cosx			! x, N(x), and ANA solution
	integer, dimension(N)	:: pr, pl			! pointer to left & right

	integer :: i						! for loop
	real	:: dx						! grid space

	real	:: fw1, bw1, cen2, cen4				! NUM solution
	real	:: diff_fw1, diff_bw1, diff_cen2, diff_cen4	! difference of NUM & ANA

	! init
	dx	= 2.0*PI/real(N)
	do i	= 1, N
		x(i)	= (real(i)-0.5)*dx			! x
		sinx(i)	= sin(x(i))				! N(x)
		cosx(i) = cos(x(i))				! ANA solution = pN/px
		pl(i)	= i-1					! to-left grid pointer
		pr(i)	= i+1					! to-right grid pointer
	end do
	pl(1) = N
	pr(N) = 1

	! compute
	do i 	= 1, N
		fw1	= (sinx(pr(i))-sinx(i))/dx
		bw1	= (sinx(i)-sinx(pl(i)))/dx
		cen2	= (sinx(pr(i))-sinx(pl(i)))/(2.0*dx)
		cen4 	= (sinx(pl(pl(i)))-8.0*sinx(pl(i))+ 		&
		           8.0*sinx(pr(i))-sinx(pr(pr(i))))/(12.0*dx)

		diff_fw1	= fw1-cosx(i)
		diff_bw1	= bw1-cosx(i)
		diff_cen2	= cen2-cosx(i)
		diff_cen4	= cen4-cosx(i)

		! to print x, N, pN/px, and NUM solutions
		print "(7f10.4)", x(i), sinx(i), cosx(i), fw1, bw1, cen2, cen4

		! to print x, N, pN/px, and differences which=(NUM-ANA)
		!print "(7f10.4)", x(i), sinx(i), cosx(i), diff_fw1, diff_bw1, diff_cen2, diff_cen4

	end do

end program
