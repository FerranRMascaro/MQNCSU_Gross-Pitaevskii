set term png font ",16"
set datafile separator ","

set title ""

input = "b.csv"
set output "b2.png"

set xrange [100:500]
set yrange [0:0.02]
unset logscale x
unset logscale y

plot input u 1:($3/$1) t"E_T (N^{2/5})" w lp pt 7