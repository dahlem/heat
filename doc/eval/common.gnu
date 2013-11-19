set terminal postscript eps enhanced colour solid rounded

# input file format is comma-separated values
set datafile separator ","

# output of the surface
set pm3d at s explicit
set grid
set key off
unset hidden3d

set palette rgbformulae 22,13,-31

set xlabel "X"
set ylabel "Y"
