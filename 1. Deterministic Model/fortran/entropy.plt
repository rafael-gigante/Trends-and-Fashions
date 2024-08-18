#Gráfico da Entropia
set terminal pngcairo dashed 
set title "N = 10^5, L = 10^3"
set output "entropy.png"
set xlabel "time [iterations]"
set ylabel "S"
set xrange [0:400]
set datafile missing NaN
plot "entropy.out" title 'entropy' with lines lt -1, \
     "log.out" title 'ln(L)'with points ps 0.2 pt 7 lc rgb "#FF1100"