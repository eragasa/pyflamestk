mass 1 24.305
mass 2 15.999

group magnesium type 1
group oxygen type 2

set group magnesium charge 1.6925590324649487
set group oxygen charge -1.6925590324649487

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1040.798071088434 0.2873457264072148 0.0 ${R_cut}
pair_coeff 2 2 5822.672917103254 0.24989971725184199 73.18991654066915 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
