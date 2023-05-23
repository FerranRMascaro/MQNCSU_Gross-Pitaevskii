set term png
set datafile separator ";"

set title ""

set xrange [100:1000000]
set yrange [0:45]
# set xlabel ""
# set ylabel ""

set key top left

input = "a)_outputs.csv"
set output "a.png"

plot input u 1:2 t"Energia" w lp pt 7, \
    input u 1:3 t"mu" w lp pt 7, \
    input u 1:4 t"Cinètica" w lp pt 7, \
    input u 1:5 t"Potencial" w lp pt 7, \
    input u 1:6 t"potho" w lp pt 7, \
    input u 1:7 t"potint" w lp pt 7, \
    input u 1:8 t"Radi" w lp pt 7, \
    input u 1:9 t"Radi2" w lp pt 7