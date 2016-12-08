mass 1 24.305
mass 2 15.999

group magnesium type 1
group oxygen type 2

set group magnesium charge 2.1380080337988705
set group oxygen charge -2.1380080337988705

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1193.8938078034512 0.31984544298092676 7.113576676484863 ${R_cut}
pair_coeff 2 2 24649.704429204125 0.3426215989680028 70.95271797556802 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
