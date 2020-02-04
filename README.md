# Ising-Model
Monte Carlo simulation of the Ising Model of Ferromagnetism

Project to make a parallel MC simulation of the Ising Model of Ferromagnetism. 

Dependencies: pcgbasic.h + pcgbasic.o and associated makefile, gnuplot and associated plotting scripts included here. If running on windows I recommend cygwin64 to compile from the command line and Xlaunch to add gnuplot functionality. Depending on the system you may need to set the display before plotting via e.g. "export DISPLAY=:0"

Progress notes: barebones simulation currently present. Progression from random spin state to statistical equilibrium can be simulated and charted for variable input parameters. Detection of equilibrium method work well for B>0, some issues for B=0 meaning the sim invariably proceeds to max iterations. Next step to fix this, parallelize and begin investiagting accuracy of statistical properties. 

