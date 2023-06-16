set term png

set title ""

set xrange [0:8]
set yrange [0:0.7]
set xlabel "r"
set ylabel "ρ"

set key top right

set output "c1000.png"

plot "1000_den.txt" u 1:2 t"N=1000" w l lw 2.5, \
     "b_results1000_den.txt" u 1:2 t"Thomas-Fermi" w l lw 2.5


set output "c100000.png"

set yrange [0:0.04]

plot "100000_den.txt" u 1:2 t"N=100000" w l lw 2.5, \
     "b_results100000_den.txt" u 1:2 t"Thomas-Fermi" w l lw 2.5