mass 1 24.305
mass 2 15.999

group magnesium type 1
group oxygen type 2

set group magnesium charge 1.7
set group oxygen charge -1.7

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 26007.0 0.21938 299.981 ${R_cut}
pair_coeff 2 2 22370.0 0.28323 17.6997 ${R_cut}

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
