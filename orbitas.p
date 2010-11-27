      # Gnuplot script file for plotting data in file "force.dat"
      # This file is called   force.p
      set   autoscale                        # scale axes automatically
      set title "Orbitas"
      set xlabel "X"
      set ylabel "Y"
     # plot     "out" using 5:6 with lines title "moon"
      plot    "outFey" using 3:4 with lines title "Earth", "outFey" using 5:6 with lines title "moon", "outFey" using 7:8 with lines title "venus"
