set terminal postscript eps enhanced colour solid rounded
set output 'real.eps'

# output of the surface
#set pm3d# at bs
#set hidden3d offset 1 trianglepattern 3 undefined 1 altdiagonal #bentover
set style data lines

# remove the box with the meaning of the color
unset colorbox

# set ranges
set isosample 50, 50

set yrange [-2:0.5]
set xrange [-0.5:2]

set grid
set key off

# remove tics on the z-axis
unset ztics


set xlabel "x"
set ylabel "y"

set xtics axis
set ytics axis

set title "Exact Solution of the Poisson Equation"

splot cos(x + y) * sin(x - y)
reset