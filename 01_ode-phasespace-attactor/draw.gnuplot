set term pdf
set output "plot.pdf"

set xlabel "X"
set ylabel "Y"

set xrange [-5:5]
set yrange [-5:5]

set title "2D Phase Space"

set size square
#plot "out.txt" t "dx/dt=1-x^2, dy/dt=1-xy" w l linecolor rgb "blue" lw 0.5
plot "out.txt" t "dx/dt=1-x^2, dy/dt=1-xy" w dots linecolor rgb "blue"
#pause -1

