set term pdf
set output "plot-absvalue.pdf"

set xlabel "x"
set ylabel "sin(x)"

set yrange [-1.1:1.1]
set xrange [0:6.28]

set title "Analytical and Numercial solutions of d(sinx)/dx"

plot "out.txt" u 1:2 t "sin(x)" w l lw 8, "out.txt" u 1:3 t "cos(x)" w l lw 8, "out.txt" u 1:4 t "Forward" w l lw 2, "out.txt" u 1:5 t "Backward" w l lw 2, "out.txt" u 1:6 t "Central-2" w l lw 2, "out.txt" u 1:7 t "Central-4" w l lw 2

#pause -1
