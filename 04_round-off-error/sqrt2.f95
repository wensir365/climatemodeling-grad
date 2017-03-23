program sqrt2

	implicit none

	real :: a, b, c, x

	a = sqrt(2.0)
	b = sqrt(2.0)
	c = a*b
	x = 2.0

	print "(f20.10, f20.10, f20.10)", a, b, c
	print "(1f20.10)", x

end program
