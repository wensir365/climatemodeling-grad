set term pdf
set output "plot.pdf"

set xlabel "Longitude (Deg_East)"
set ylabel "Z500 (gpm)"

set yrange [5000:6000]
set xrange [0:360]

set title "Truncation on Geopotential Height at 500mb (Mar/23/2017)"

#plot "out" u 1:2 w l lw 12 t "Z500", "out" u 1:3 w l t "W# 0", "out" u 1:4 w l t "W# 1", "out" u 1:5 w l t "W# 2", "out" u 1:6 w l t "W# 3", "out" u 1:7 w l t "W# 5", "out" u 1:8 w l lw 2 t "W# 10","out" u 1:9 w l lw 2 t "W# 20", "out" u 1:10 w l lw 2 t "W# 30"

plot "out" u 1:2 w l lw 20 linecolor rgb "cyan" t "Z500", "out" u 1:3 w l linecolor rgb "black" t "W# 0", "out" u 1:4 w l linecolor rgb "green" t "W# 1", "out" u 1:6 w l linecolor rgb "purple" t "W# 3", "out" u 1:7 w l linecolor rgb "brown" t "W# 5", "out" u 1:8 w l lw 2 linecolor rgb "red" t "W# 10","out" u 1:9 w l lw 2 linecolor rgb "blue" t "W# 20"

#pause 5
