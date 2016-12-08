mass 0 24.305
mass 1 15.999

group magnesium type 0
group oxygen type 1

set group magnesium charge 2.0
set group oxygen charge -2.0

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 0.1069023957054892 17943.296064682618 65.44108967241823 ${R_cut}
pair_coeff 2 2 0.34019719310916274 17898.401701213337 50.17309589384038 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
