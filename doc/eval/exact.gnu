set terminal postscript eps enhanced colour solid rounded
set output 'exact.eps'

# change the angle the plot is shown
set view ,, 0.7, 1

set isosample 50, 50

set yrange [-2:0.5]
set xrange [-0.5:2]

set grid
set key off

set xlabel "X"
set ylabel "Y"

splot cos(x + y) * sin(x - y)
reset
