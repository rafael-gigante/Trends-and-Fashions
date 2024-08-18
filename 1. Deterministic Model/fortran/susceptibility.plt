#Gráfico da suscetibilidade
set terminal pngcairo dashed
set title "N = 10^5, L = 10^3"
set output "susceptibility.png"
set xlabel "time [iterations]"
set ylabel "S_c / N²"
set xrange [200:400]
set yrange [0:0.1]
set datafile missing NaN
unset key
plot "susceptibility.out" title 'Pc / N' with lines lt -1