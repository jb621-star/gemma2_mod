#!/bin/bash

#module load R
export R_LIBS=~/Rlibs:${R_LIBS}
#Rscript -e "install.packages('readxl', '~/Rlibs', 'https://cran.r-project.org')"
#Rscript -e "install.packages('factoextra', '~/Rlibs', 'https://cran.r-project.org')"
#Rscript -e "install.packages('gemma2', '~/Rlibs', 'https://cran.r-project.org')"
#Rscript -e "install.packages('car', '~/Rlibs', 'https://cran.r-project.org')"


#sbatch -p scavenger --mem=100G --out=calcQi_mod.out --wrap="module load R/4.1.1-rhel8 && Rscript --verbose --vanilla calcQi_mod.R"

#sbatch -p scavenger --mem=150G --out=calcXHiY.out --wrap="module load R/4.1.1-rhel8 && Rscript --verbose --vanilla calcXHiY.R"

#sbatch -p scavenger --mem=150G --out=MphCalcLogL.out --wrap="module load R/4.1.1-rhel8 && Rscript --verbose --vanilla MphCalcLogL.R"

#sbatch -p scavenger --mem=150G --out=update_U_E_quants.out --wrap="module load R/4.1.1-rhel8 && Rscript --verbose --vanilla update_U_E_quants.R"

#sbatch -p scavenger --mem=150G --out=calc_epsilon.out --wrap="module load R/4.1.1-rhel8 && Rscript --verbose --vanilla calc_epsilon.R"

sbatch -p scavenger --mem=150G --out=update_v.out --wrap="module load R && Rscript --verbose --vanilla update_v.R"

#sbatch -p scavenger --mem=150G --out=MphCalcP.out --wrap="module load R && Rscript --verbose --vanilla MphCalcP.R"

#sbatch -p scavenger --mem=200G --out=process_geno.out --wrap="module load R/4.1.1-rhel8 && Rscript --verbose --vanilla process_genotypes_n_counts.R"

#sbatch -p scavenger --mem=200G --out=prune_geno.out --wrap="module load R/4.1.1-rhel8 && Rscript --verbose --vanilla prune_genotypes.R"

#sbatch -p scavenger --mem=200G --out=strain_freq-inf.out --wrap="module load R/4.1.1-rhel8 && Rscript --verbose --vanilla strain_frequency_inference.R"
