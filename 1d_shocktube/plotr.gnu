set term png giant truecolor \
    font arial size 1200,800

set output "density.png"

#set palette gray
#MAX=GPVAL_Y_MAX
#MIN=GPVAL_Y_MIN

set xrange[-5:5]
#set yrange [MIN-(MAX-MIN)*0.05:MAX+(MAX-MIN)*0.05]
#set yrange[-1:4]

set xlabel "x"
set ylabel "density"

set pointsize 0.5
set grid

plot 'exact.dat' using 1:3 w l lw 2.0 lc rgb 'red' title "exact", 'output.dat' using 1:3 pt 6 ps 0.5 lc rgb 'blue' title "numerical"
