set term pdf
set output "plot.pdf"

set xlabel "X"
set ylabel "Y"
set zlabel "Z"

set title "3D Phase Space"
splot "out.txt" w l lw 1 linecolor rgb "red" t "Lorenz 63 Model"
#pause -1
