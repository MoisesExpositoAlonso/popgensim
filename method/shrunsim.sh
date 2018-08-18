#!/bin/bash

#cd ..

# Figure 2
Rscript method/runsimulation.R --p 0.5 --b 1 --a 0.1 --ss 0 --svar 0.1 --mu 1 --m 100 --epi 1 --hidden 0.2 --training 0.8 --iterations 1e5 &

# Figure S3  -> no hidden
Rscript method/runsimulation.R --p 0 --b 1 --a 0.1 --ss 0 --svar 0.1 --mu 1 --m 100 --epi 1 --hidden 0. --training 0.8 --iterations 1e5 &

# Figure S3  -> no mortaility
Rscript method/runsimulation.R --p 0 --b 1 --a 0.1 --ss 0 --svar 0.1 --mu 1 --m 100 --epi 1 --hidden 0.2 --training 0.8 --iterations 1e5 &


# Figure S4 -> stronger SNP effects
Rscript method/runsimulation.R --p 0.5 --b 1 --a 0.1 --ss 0.8 --svar 0.2 --mu 1 --m 100 --epi 1 --hidden 0.2 --training 0.8 --iterations 1e5 &


# Figure S5 -> stronger SNP effects and nomortality
Rscript method/runsimulation.R --p 0 --b 1 --a 0.1 --ss 0.8 --svar 0.2 --mu 1 --m 100 --epi 1 --hidden 0.2 --training 0.8 --iterations 1e5 &

# Figure S6 -> With epistatic effect
Rscript method/runsimulation.R --p 0.5 --b 1 --a 0.1 --ss 0.8 --svar 0.2 --mu 1 --m 100 --epi 1.5 --hidden 0.2 --training 0.8 --iterations 1e5 &

# Figure S7 -> With many loci
Rscript method/runsimulation.R --p 0.5 --b 1 --a 0.1 --ss 0.8 --svar 0.1 --mu 1 --m 1000 --epi 1 --hidden 0.2 --training 0.8 --iterations 1e5 &

# Figure S7 -> With many many loci
Rscript method/runsimulation.R --p 0.5 --b 1 --a 0.1 --ss 0.8 --svar 0.1 --mu 1 --m 5000 --epi 1 --hidden 0.2 --training 0.8 --iterations 1e5 &

