set term pdf
set output "plot-diff.pdf"

set xlabel "x"
set ylabel "sin(x)"

set yrange [-0.4:0.4]
set xrange [0:6.28]

set title "Difference between Numerical and Analytical solutions of d(sinx)/dx"

plot "out.txt" u 1:4 t "Forward" w l lw 4, "out.txt" u 1:5 t "Backward" w l lw 4, "out.txt" u 1:6 t "Central-2" w l lw 4, "out.txt" u 1:7 t "Central-4" w l lw 4

#pause -1

