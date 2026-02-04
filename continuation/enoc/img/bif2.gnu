reset

set terminal epslatex monochrome "default" 11
set output '/home/dave/TeX/ENOC/mag/bif2_1.eps'

set nokey

set mxtics 2
set mytics 2

set size ratio 9/16

set xtics nomirror
set ytics nomirror

set xrange [-0.05 : 0.15]

set xlabel '$\lambda$'
set ylabel '$D$'

plot \
'../data/b.mag.4' u 5:12 w l lt 1 lw 2 lc 1,\
'../data/b.mag.4-' u 5:12 w l lt 1 lw 2 lc 1,\
'../data/b.mag.4+' u 5:12 w l lt 1 lw 2 lc 1
