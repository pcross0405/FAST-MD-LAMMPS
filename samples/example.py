from FAST_MD import 
from lammps import lammps

# create lammps object

lmp = lammps(cmdargs = ['-screen', 'none']) 

# load and run input file

lmp.file('IrSi.in')

# find number of isolobal bonds needed for while loop
# initialize iso_limit variable

iso_limit = 0

# check for atom type = 1 since transition metal is type 1

for i in lmp.numpy.extract_atom('type'):
   if i == 1:
      iso_limit += 1

# define counter for while loop

n = 0

# define list to log which pairs have an isolobal atom in between

pairs = []

# define list to log id of all isolobal atoms added

added_atoms = []

# find all pairs

pair_ids = find_pairs(pairs, lmp)

# move cursor to correct spot for print progress bar

print('\n\n')

# initial progress bar at 0 completion

progress_bar(0, iso_limit, 0, 0, 0, prefix = 'Progress', suffix = '')

# minimize box before starting

box_relax(lmp)

# loop to add isolobal atoms until the 18 - n limit is reached

while n < iso_limit/2:
   
# create new id number for isolobal atom 

   new_id = lmp.get_natoms() + 1

# find coordinates for isolobal atoms

   positions = find_centers(pair_ids, lmp)

# find coordinate of lowest energy isolobal atom

   iso_coord_and_nrg = iso_loop(new_id, positions, lmp, n, iso_limit)

# iso_loop function returns [[x,y,z],min_nrg]
# iso_coord_and_nrg[0] = [x,y,z] aka coordinates of added isolobal atom

   iso_coord = iso_coord_and_nrg[0]

# iso_coord_and_nrg[1] = min_nrg aka current potential energy w/ isolobal atom added in
# needed for progress bar

   min_nrg = iso_coord_and_nrg[1]

# add lowest energy atom to box

   lmp.create_atoms(1, [new_id], [3], [iso_coord[0], iso_coord[1], iso_coord[2]])   

# add pair ids to list of pairs with isolobal atom in between

   pairs.append(pair_ids[positions.index(iso_coord)])

# remove transition metal atom pair id 

   pair_ids.remove(pair_ids[positions.index(iso_coord)])

# update added_atoms list with newly added isolobal atom id

   added_atoms.append(new_id)

# update lammps group to incorporate all added isolobal atoms

   lmp.command('group AddedAtoms type 3') 

# remove forces from isolobal atoms

   if n == 0:
      lmp.command('fix NoForce AddedAtoms setforce 0.0 0.0 0.0')

# minimize current box     

   lmp.commands_string('''
   min_style quickmin                    
   minimize 0 1e-10 1000 10000
   ''')

# get current step for progress bar

   step = lmp.get_thermo('step')

# make sure after final minimization that all isolobal atoms are centered between transition metal pairs

   check_center(pairs, added_atoms, lmp, n, min_nrg, iso_limit, new_id)     

# increment counter

   n += 1

# update progress as loop iterations finish

   progress_bar(n, iso_limit/2, int(step), min_nrg, new_id, prefix = 'Progress', suffix = '')

# close out lammps object 
                  
lmp.close()
