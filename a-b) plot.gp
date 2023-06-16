set term png font ",16"
set datafile separator ","

set title ""

set xrange [100:1000000]
# set yrange [0:1]
set yrange [0.04:50]
set format x '10^{%T}'
set format y '10^{%T}'
set logscale x
set logscale y
set xlabel "N"
set ylabel "E_i , μ_i"

set key top left

input = "results.csv"
set output "a.png"

# plot input u 1:($4/42.11922607982219) t"μ/μ_{10^6}" w lp pt 7, \
#     input u 1:($3/30.120102109535193) t"E/E_{10^6}" w lp pt 7, \
#     input u 1:($5/0.6564914216057008) t"E^{kin}/E^{kin}_{10^2}" w lp pt 7, \
#     input u 1:($6/30.058862952430413) t"V/V_{10^6}" w lp pt 7, \
#     input u 1:($7/18.059738982143504) t"V^{osc}/V^{osc}_{10^6}" w lp pt 7, \
#     input u 1:($8/11.99912397028691) t"V^{int}/V^{int}_{10^6}" w lp pt 7, \
#     input u 1:($9/6.00994824971788) t"r/r_{10^6}" w lp pt 7, \
#     input u 1:($10/36.11947796428701) t"r2/r2_{10^6}" w lp pt 7

plot input u 1:4 t"μ" w lp pt 7, \
    input u 1:6 t"μ_T" w lp pt 7, \
    input u 1:3 t"E_T" w lp pt 7, \
    input u 1:5 t"E_{cin.}" w lp pt 7, \
    input u 1:7 t"E_{osc.}" w lp pt 7, \
    input u 1:8 t"μ_{int.}" w lp pt 7
    # input u 1:9 t"r" w lp pt 7, \
    # input u 1:10 t"r2" w lp pt 7


input = "b.csv"
set output "b.png"

set yrange [0.1:50]

f(x) = a*x**c
fit f(x) input using 1:3 via a,c

plot input u 1:4 t"μ" w lp pt 7, \
    input u 1:6 t"μ_T" w lp pt 7, \
    input u 1:3 t"E_T" w lp pt 7, \
    input u 1:7 t"E_{osc.}" w lp pt 7, \
    input u 1:8 t"μ_{int.}" w lp pt 7
    # input u 1:5 t"E_{cin.}" w lp pt 7, \
    # input u 1:9 t"r" w lp pt 7, \
    # input u 1:10 t"r2" w lp pt 7