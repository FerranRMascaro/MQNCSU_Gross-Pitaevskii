set term png

set title ""

set xrange [0.1:6]
set yrange [0:3]
# set xlabel ""
# set ylabel ""

# set key top/bottom left/right

set output "fort9.png"

plot "fort.9" index 0 u 1:2 t"" w lp pt 7

set output "fort10.png"

set xrange [0.1:6]
set yrange []

plot "fort.10" index 0 u 1:2 t"" w lp pt 7