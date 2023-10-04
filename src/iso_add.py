from lammps import LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR

# fxn for adding an isolobal atom to each possible isolobal site

def iso_add(x, y, z, new_id, lmp):
   
# create dummy atom at x,y,z

   lmp.create_atoms(1, [new_id], [3], [x,y,z])

# recalc global properties after adding new atom with a run of 0 timesteps

   lmp.commands_string('''
   compute TotalPE all pe
   run 0                 
   ''')

# extract total potential compute

   nrg = lmp.numpy.extract_compute('TotalPE', LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)

# delete compute before calling again

   lmp.command('uncompute TotalPE')

# return new total pe after adding isolobal atom

   return nrg