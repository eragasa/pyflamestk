mass 1 24.305
mass 2 15.999

group magnesium type 1
group oxygen type 2

set group magnesium charge 1.571726425841593
set group oxygen charge -1.571726425841593

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1251.1791881737313 0.3070349710035269 0.0 ${R_cut}
pair_coeff 2 2 23368.223007503744 0.2392353000140214 54.26020541205216 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
