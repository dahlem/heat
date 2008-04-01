set terminal postscript eps enhanced colour solid rounded
set output 'poiss.eps'

# input file format is comma-separated values
set datafile separator ","

# change the angle the plot is shown
set view 75, 30, 1, 0.5

# output of the surface
set pm3d at s explicit
set grid
set key off
unset hidden3d

set palette rgbformulae 22,13,-31

set xlabel "X"
set ylabel "Y"

set size 1.3,1

set multiplot layout 1,2
set origin 0,0
set colorbox user size 0.02,0.2 origin 0.5,0.3
set title "SERIAL"
splot 'serial.dat' u 1:2:3:4 w pm3d

set cblabel "ERR"
set origin 0.65,0
set colorbox user size 0.02,0.2 origin 1.15,0.3
set title "PARALLEL"
splot 'parallel.dat' u 1:2:3:4 w pm3d

unset multiplot
reset
