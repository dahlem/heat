set terminal postscript eps enhanced colour solid rounded
set output 'poiss.eps'

# input file format is comma-separated values
set datafile separator ","

# choose quality
#set samples 50, 50
#set isosamples 100

# change the angle the plot is shown
#set view 55, 45, 1, 2

# output of the surface
#set pm3d# at bs
#set hidden3d offset 1 trianglepattern 3 undefined 1 altdiagonal #bentover
set style data lines

# remove the box with the meaning of the color
unset colorbox

# set ranges
#set zrange [-1.1:5.8]
#set yrange [-1:1]
#set xrange [-2:2]

set grid
set key off

# remove tics on the z-axis
unset ztics


set xlabel "x"
set ylabel "y"

set xtics axis
set ytics axis

set title "Poisson Equation"

splot 'poiss.dat' w l ls 1#, 'poiss.dat' w pm3d
reset