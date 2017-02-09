mass 1 24.305
mass 2 15.999

group magnesium type 1
group oxygen type 2

set group magnesium charge 2.1335249445945137
set group oxygen charge -2.1335249445945137

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 849.3091599660628 0.3025359510724032 0.0 ${R_cut}
pair_coeff 2 2 23489.37888205462 0.3558498499632836 76.7448985072756 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
# setup minimization style
min_style cg
min_modify dmax ${dmax} line quadratic

# setup output
thermo 1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
