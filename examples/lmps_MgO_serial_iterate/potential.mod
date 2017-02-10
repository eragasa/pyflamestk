mass 1 24.305
mass 2 15.999

group magnesium type 1
group oxygen type 2

set group magnesium charge 1.7842000914946492
set group oxygen charge -1.7842000914946492

variable R_cut equal 10.0

pair_style buck/coul/long ${R_cut}
pair_coeff 1 1 0.0 0.5 0.0 ${R_cut}
pair_coeff 1 2 668.1794481015199 0.30803597064512583 0.0 ${R_cut}
pair_coeff 2 2 50062.959616012406 0.1571533490914843 46.64299548222826 ${R_cut}

kspace_style pppm 1.0e-5

neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
