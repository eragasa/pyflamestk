import pyflamestk.gulp as gulp

sim_name = 'test'
sim_dir = 'test_dir'
structure_name = 'structure_name'
fname_structure_vasp = 'MgO_NaCl_111.vasp'

# read lines
f = open('regress_lammps.out')
lines = f.readlines()
f.close()

# process each linels

names = [v.strip() for v in lines[0].split(',')]
types = [v.strip() for v in lines[1].split(',')]
values = [lines[i] for i in range(2,len(lines))]
for i,l in enumerate(values):
    values[i] = [float(v.strip()) for v in l.split(',')]
# <-------
gulp_pot = gulp.BuckinghamPotential(["Mg","O"])
for i,t in enumerate(types):
    if t == 'param':
        name = names[i]
        value = values[0][i]
        gulp_pot.param_dict[name] = value

#print(gulp_pot.parameter_names)
#print(gulp_pot.get_potential_string())

gulp_sim = gulp.Simulation(sim_name,
                           sim_dir,
                           structure_name)
gulp_sim.potential = gulp_pot
gulp_sim.fname_structure_vasp = fname_structure_vasp
gulp_sim.write_input_file()
#gulp_sim.create()
