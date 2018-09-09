#!/bin/bash

cd ..

Rscript method/sim_gws.R --genome top --sparsity 0.9 --phenotype additive --heritability 1 & 
Rscript method/sim_gws.R --genome top --sparsity 0.9 --phenotype epistatic --heritability 1 & 
Rscript method/sim_gws.R --genome top --sparsity 0.5 --phenotype additive --heritability 1 & 
Rscript method/sim_gws.R --genome top --sparsity 0.5 --phenotype epistatic --heritability 1 & 

