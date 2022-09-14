# ~\~ language=Gnuplot filename=scripts/adv-diff.gnuplot
# ~\~ begin <<lit/fiadeiro-veronis.md|scripts/adv-diff.gnuplot>>[init]
set term svg
set multiplot
set lmargin at screen 0.1
set size 1, 0.35
set yrange [-4:0]
set ylabel "log_{10}(c)"
set xlabel "x"
set key box opaque bottom
plot 'data/adv-diff.out' i 0 u 1:(log10($2)/2) w l lc 6 t'diffusion', \
     '' i 0 u 1:(log10($3)) w l lc 7 t'advection'

set origin 0, 0.32
set size 1, 0.7
unset logscale y
set yrange [-0.1:1.05]
set ylabel "y"
unset xlabel
set label "← advection dominated" at -0.95,-0.05 left
set label "diffusion dominated →" at 0.95,-0.05 right
unset xtics
set key box opaque top
plot 'data/adv-diff.out' i 1 t"t=0.0" w l lc 1, \
     '' i 2 t"t=0.5" w l lc 2, '' i 3 t"t=1.0" w l lc 3
unset multiplot
# ~\~ end
