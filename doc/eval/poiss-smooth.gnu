load "common.gnu"

# change the angle the plot is shown
set view ,, 0.7, 1

set output 'poiss-smooth.eps'
splot 'serial.dat' u 1:2:3:4 w pm3d
reset
