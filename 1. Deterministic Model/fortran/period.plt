#Gráfico do decaimento das labels
set terminal pngcairo dashed 
set output "period.png"
set xlabel "índice de período"
set ylabel "período de oscilação"
set xrange [0:35]
set yrange[15:25]
unset key
plot "period.out" with points ps 1 pt 7 lc rgb "#000000"