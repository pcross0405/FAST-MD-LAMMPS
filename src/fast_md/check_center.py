import numpy as np
from lammps import LMP_TYPE_VECTOR, LMP_STYLE_ATOM
from . import find_centers, progress_bar

# function for checking if isolobal atoms remain in center of transition metals

def check_center(tm_pairs, added_atoms, lmp, n, min_nrg, iso_limit, new_id):
   
# find current centers of pairs of transition metals with isolobal atom in between

   current_centers = find_centers.find_centers(tm_pairs, lmp)

# loop through pair centers

   for i, center in enumerate(current_centers):
      
# compute x,y,z coordinates of isolobal atom in between transition metal pair

        lmp.commands_string(f'''
        group IsolobalAtom id {added_atoms[i]}
        compute IsolobalAtomX IsolobalAtom property/atom x
        compute IsolobalAtomY IsolobalAtom property/atom y
        compute IsolobalAtomZ IsolobalAtom property/atom z
        run 0
        ''')

# extract coordinates of isolobal atom

        iso_atom_x = [x for x in lmp.numpy.extract_compute('IsolobalAtomX', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if x != 0][0]
        iso_atom_y = [y for y in lmp.numpy.extract_compute('IsolobalAtomY', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if y != 0][0]
        iso_atom_z = [z for z in lmp.numpy.extract_compute('IsolobalAtomZ', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if z != 0][0]

# list comprehension above removes zeros, add 0 back if coordinate happens to be zero

        if iso_atom_x == []:
            iso_atom_x = [0.0]
        if iso_atom_y == []:
            iso_atom_y = [0.0]
        if iso_atom_z == []:
            iso_atom_z = [0.0]

# delete computes and group before next loop iteration

        lmp.commands_string('''
        uncompute IsolobalAtomX
        uncompute IsolobalAtomY
        uncompute IsolobalAtomZ
        group IsolobalAtom delete                                                                               
        ''')

# find difference between true pair center and isolobal atom coordinates

        dx = center[0] - iso_atom_x
        dy = center[1] - iso_atom_y
        dz = center[2] - iso_atom_z

# compute distance between atom and true center

        check_distance = ((dx)**2 + (dy)**2 + (dz)**2)**(1/2)

# check if distance is greater than cutoff

        if check_distance > 0.01:

# move isolobal atom back to true center if distance is greater than cutoff
# minimize after moving atom 

            lmp.commands_string(f'''
            group MoveAtom id {added_atoms[i]}
            displace_atoms MoveAtom move {dx} {dy} {dz}
            min_style sd
            minimize 0 1e-10 1000 10000
            group MoveAtom delete
            ''')

# get current step

            step = lmp.get_thermo('step')

# update progress bar

            progress_bar.progress_bar(n, iso_limit, int(step), min_nrg, new_id, prefix = 'Progress', suffix = '')

# rerun function to see if atom stayed centered after minimization

            check_center([tm_pairs[i]], [added_atoms[i]], lmp, n, min_nrg, iso_limit, new_id)
