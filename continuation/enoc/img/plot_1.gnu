reset

#set terminal epslatex color '"default" 11'
#set output '/home/dave/TeX/ENOC/mag/bif2.eps'

set nokey

set mxtics 2
set mytics 2

set size ratio 9/16

set xtics nomirror
set ytics nomirror

# set xrange [0:0.035]

set xlabel '$\lambda$'
set ylabel '$D$'

plot \
'../data/b.mag.4' u 5:12 w l 1,\
'../data/b.mag.4-' u 5:12 w l 1,\
'../data/b.mag.4+' u 5:12 w l 1
