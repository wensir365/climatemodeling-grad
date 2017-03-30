#set term pdf
#set output "plot-absvalue.pdf"

set xlabel "x"
set ylabel "y"

set yrange [-0.2:1.0]
set xrange [0:2*pi]

set title "1D Advection Eq Solutions"

plot "out" u 1:2 t "Analytic" w l lw 10, "out" u 1:3 t "FD" w l lw 2, "out" u 1:4 t "Pesudo-SM" w l lw 4

pause 5
