mass 1 24.305
mass 2 15.999

group magnesium type 1
group oxygen type 2

set group magnesium charge 2.1371937311693414
set group oxygen charge -2.1371937311693414

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 831.2251655175153 0.30799390650175645 0.0 ${R_cut}
pair_coeff 2 2 23138.612152302983 0.26365421131443306 64.82394130917247 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
