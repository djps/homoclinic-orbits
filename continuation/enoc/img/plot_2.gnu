# plot \
#'../data/b.mag.5_new' u 5:12 w l lt 1,\
#'../datab.mag.5_new-' u 5:12 w l lt 1,\
#'../datab.mag.5_new+' u 5:12 w l lt 1,\
#'../datab.mag.5x' u 5:12 w l lt 2,\
#'../datab.mag.5x-' u 5:12 w l lt 2 

set terminal epslatex color "default" 11
set output '/home/dave/TeX/ENOC/mag/codim22.eps'

set nokey

set mxtics 2
set mytics 2

set size ratio 9/16

set xtics nomirror
set ytics nomirror

# set xrange [0:0.035]

set xlabel '$\lambda$'
set ylabel '$D$'

#plot '../data/b.stationary.dat' u 5:12 w l lt 1 

plot '../data/codim22-.dat' u 5:12 w l lt 1, '../data/codim22+.dat' u 5:12 w l lt 1, '../data/codim22.dat' u 5:12 w l lt 1
