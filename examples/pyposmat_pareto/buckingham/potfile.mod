mass 0 24.305
mass 1 15.999

group magnesium type 0
group oxygen type 1

set group magnesium charge 2.0
set group oxygen charge -2.0

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 0.29896507169690045 14631.148352824932 42.22881193884092 ${R_cut}
pair_coeff 2 2 0.1554572309938736 3598.2776740679305 13.014032120857236 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
