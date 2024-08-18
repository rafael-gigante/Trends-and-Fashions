#Gerador do gif
set terminal gif animate delay 50
set output "modelo.gif"
stats "data.out" name "A" nooutput
set xrange [1:1000]
set yrange[0:2000]
set xlabel "índice da tendência"
set ylabel "número de pessoas"
do for [i=0:(A_blocks)]{
   plot 'data.out' index i w impulses lw 2.5 lc 'black' notitle
}