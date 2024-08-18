#Gráfico da Percolação
set terminal pngcairo dashed
set title "N = 10^5, L = 10^3"
set output "percolation.png"
set xlabel "time [iterations]"
set ylabel "P_c / N"
set xrange [200:400]
set yrange [0:1]
set datafile missing NaN
unset key
plot "percolation.out" title 'Pc / N' with lines lt -1