program smalldiff

	implicit none

	real :: a, b, diff

	a = 1.2345678
	b = 1.2345673
	diff = a-b

	print "(3f20.10)", a,b,diff

end program
