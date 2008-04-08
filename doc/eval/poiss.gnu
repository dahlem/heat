load "common.gnu"

set output 'poiss.eps'

set multiplot layout 1,2 scale 1,0.7
set title "SERIAL"
splot 'serial.dat' every 10:10 u 1:2:3:4 w pm3d

set cblabel "ERR"
set title "PARALLEL"
splot 'parallel.dat' every 10:10 u 1:2:3:4 w pm3d

unset multiplot
reset
