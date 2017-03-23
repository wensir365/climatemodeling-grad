program bigplussmall

	implicit none

	real :: big, small, plus, minus
	integer :: i

	big = 1.23456789
	small = 0.00000001
	plus = big + small
	minus = big - small

	print "(3f20.10)", big, small, plus
	print "(3f20.10)", big, small, minus
	print "(1f20.10)", 1.2345679
	print *, "Another Test:"
	print "(1f20.10)", plus+small*100000

	do i = 1, 100000
		plus = plus + small
	end do

	print "(1f20.10)", plus

end program
