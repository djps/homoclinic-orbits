# plot \
#'../data/b.mag.5_new' u 5:12 w l lt 1,\
#'../data/b.mag.5_new-' u 5:12 w l lt 1,\
#'../data/b.mag.5_new+' u 5:12 w l lt 1,\
#'../data/b.mag.5x' u 5:12 w l lt 2,\
#'../data/b.mag.5x-' u 5:12 w l lt 2 

set terminal epslatex monochrome "default" 11
set output '/home/dave/TeX/ENOC/mag/stationary_1.eps'

set nokey

set mxtics 2
set mytics 2

set size ratio 9/16

set xtics nomirror
set ytics nomirror

#set xrange [-0.04 : 0.04]

set xlabel '$\lambda$'
set ylabel '$D$'

plot '../data/b.stationary.dat' u 5:12 w l lt 1 lw 2

#plot '../data/codim22-.dat' u 5:12 w l lt 1, '../data/codim22+.dat' u 5:12 w l lt 1, '../data/codim22.dat' u 5:12 w l lt 1
