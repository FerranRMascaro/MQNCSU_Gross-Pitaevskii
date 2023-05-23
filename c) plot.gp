set term png

set title ""

set xrange [0:9]
set yrange [0:0.7]
# set xlabel ""
# set ylabel ""

set key top right

set output "c.png"

plot "1000.9" u 1:2 t"N=1000" w lp pt 7, \
     "100000.9" u 1:2 t"N=100000" w lp pt 7