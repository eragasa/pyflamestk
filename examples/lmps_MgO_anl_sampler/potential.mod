mass 1 24.305
mass 2 15.999

group magnesium type 1
group oxygen type 2

set group magnesium charge 2.0611295138069776
set group oxygen charge -2.0611295138069776

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 1199.8568533318069 0.2917778937980834 0.0 ${R_cut}
pair_coeff 2 2 4745.2067211230715 0.34039368784625 38.18657310809654 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes