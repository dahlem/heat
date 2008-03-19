set terminal postscript eps enhanced colour solid rounded
set output 'poiss.eps'

# input file format is comma-separated values
set datafile separator ","

# change the angle the plot is shown
set view 60, 30, 1, 1

# output of the surface
set pm3d at s explicit
set grid
set key off

set xlabel "x"
set ylabel "y"

set title "Poisson Equation"

splot 'poiss.dat' u 1:2:3:4 w pm3d
reset
