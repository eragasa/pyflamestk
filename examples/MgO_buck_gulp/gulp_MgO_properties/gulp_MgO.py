import pyflamestk.gulp as gulp

fname_potential_file = 'pot_param.dat'
gulp_bin = '/Users/eugeneragasa/bin/gulp'
sim_name = 'test'
sim_dir = 'test_dir'
structure_name = 'structure_name'
fname_structure_vasp = 'MgO_NaCl_111.vasp'

# read in potentials
f = open(fname_potential_file)
lines = f.readlines()
f.close()

# process each linels
# for the input file,
# line 1: contains the names of the values
# line 2: contains the type of the values
# line 3: contains the value of the value

names = [v.strip() for v in lines[0].split(',')] # read in first line
types = [v.strip() for v in lines[1].split(',')] # read in second line
values = [lines[i] for i in range(2,len(lines))] # read in third line

# turn values into datatypes we can actually use
for i,v in enumerate(values):
    print(values[i])
    values[i] = v.strip().split(',')
    print(values[i])
    for j,val in enumerate(values[i]):
       if types[j] in ['param','qoi','err']:
          values[i][j] = float(val)
       elif types[j] == 'elements':
          values[i][j] = val.split(' ')
    print(values[i])

#
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
