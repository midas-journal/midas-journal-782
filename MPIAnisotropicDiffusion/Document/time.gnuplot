set terminal postscript eps color solid rounded 
set output "time.eps"
set xlabel "Nodes"
set ylabel "Time (s)"
plot \
  "time_50slices.txt"  using 1:2 title "50 Slices"  w l, \
  "time_100slices.txt" using 1:2 title "100 Slices" w l, \
  "time_200slices.txt" using 1:2 title "200 Slices" w l
