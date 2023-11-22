from lammps import lammps
from fast_md import find_centers, find_pairs, iso_loop, progress_bar, check_center

# create lammps object

lmp = lammps(cmdargs = ['-screen', 'none']) 

# load and run input file

lmp.file('NiSi2.in')

# initialize variable for number of isolobal bonds

iso_limit = 0

# count number of type 1 atoms to find total number of Ni atoms
# number Ni atoms = number of desired isolobal bonds

for i in lmp.numpy.extract_atom('type'):
   if i == 1:
      iso_limit += 1

# initialize variable for counting loop iterations
# initialize list for adding Ni atom pairs to after an isolobal atom has been placed inbetween Ni atoms
# initialize list for tracking isolobal atoms added

n = 0
pairs_w_iso_atom = []
added_atoms = []

# find all Ni pairs within specified distance

min_dist = 3.1
max_dist = 3.9
pair_ids = find_pairs.find_pairs(pairs_w_iso_atom, min_dist, max_dist, lmp)

# move cursor to correct spot for print progress bar

print('\n\n')

# initial progress bar at 0 completion

progress_bar.progress_bar(0, iso_limit, 0, 0, 0, prefix = 'Progress', suffix = '')

# loop until all desired isolobal atoms are added

while n < iso_limit:
   
   # create new id number for isolobal atom 

   new_id = lmp.get_natoms() + 1

   # find coordinates for isolobal atoms

   positions = find_centers.find_centers(pair_ids, lmp)

   # find coordinate of lowest energy isolobal atom

   iso_coord_and_nrg = iso_loop.iso_loop(new_id, positions, lmp, n, iso_limit)
   iso_coord = iso_coord_and_nrg[0]
   min_nrg = iso_coord_and_nrg[1]

   # add lowest energy atom to box

   lmp.create_atoms(1, [new_id], [3], [iso_coord[0], iso_coord[1], iso_coord[2]])   

   # add Ni pair with isolobal atom to list

   pairs_w_iso_atom.append(pair_ids[positions.index(iso_coord)])

   # remove Ni metal atom pair id 

   pair_ids.remove(pair_ids[positions.index(iso_coord)])

   # update added atoms list with newly added isolobal atom

   added_atoms.append(new_id)

   # remove forces from isolobal atoms

   lmp.command('group AddedAtoms type 3') 
   if n == 0:
      lmp.command('fix NoForce AddedAtoms setforce 0.0 0.0 0.0')

   # minimize current box     

   lmp.commands_string('''
   min_style quickmin                    
   minimize 0 1e-10 1000 10000
   ''')
   step = lmp.get_thermo('step')

   # check if isolobal atoms remained centered
   
   check_center.check_center(pairs_w_iso_atom, added_atoms, lmp, n, min_nrg, iso_limit, new_id)

   # increment counter

   n += 1

   # update progress as loop iterations finish

   progress_bar.progress_bar(n, iso_limit, int(step), min_nrg, new_id, prefix = 'Progress', suffix = '')

