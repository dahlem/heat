load "common.gnu"

# change the angle the plot is shown
set view ,, 0.7, 1

set output 'poiss_1e-12.eps'
set cblabel "ERR"
splot 'serial_1e-12.dat' every 10:10 u 1:2:3:4 w pm3d
reset
