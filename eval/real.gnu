set terminal postscript eps enhanced colour solid rounded
set output 'real.eps'

# change the angle the plot is shown
set view 60, 30, 1, 1

# set ranges
set isosample 50, 50

set yrange [-2:0.5]
set xrange [-0.5:2]

set grid
set key off

set xlabel "x"
set ylabel "y"

set title "Exact Solution of the Poisson Equation"

splot cos(x + y) * sin(x - y)
reset