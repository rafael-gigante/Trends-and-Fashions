#Gráfico do decaimento das labels
set terminal pngcairo dashed 
set output "trend_decay.png"
set xlabel "tempo"
set ylabel "taxa de ocupação da tendência dominante"
set xrange [0:100]
set yrange [0:1]
unset key
plot "label_decay.out" with linespoints ps 1 pt 7 lc rgb "#000000"