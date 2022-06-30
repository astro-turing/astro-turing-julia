# ~\~ language=Gnuplot filename=scripts/advection.gnuplot
# ~\~ begin <<lit/upwind-scheme.md|scripts/advection.gnuplot>>[init]
set term svg
set yrange [-0.05:1.05]
plot 'data/advection.out' i 0 t"t=0.0" w l lc 1, \
     '' i 1 t"t=0.5" w l lc 2, '' i 2 t"t=1.0" w l lc 3
# ~\~ end
