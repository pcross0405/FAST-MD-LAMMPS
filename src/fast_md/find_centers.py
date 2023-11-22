import numpy as np
from lammps import LMP_TYPE_VECTOR, LMP_STYLE_ATOM

# function for finding all positions at the center of transition metal contacts

def find_centers(pair_ids, lmp):

# grab lengths of simulation box

   xmax = lmp.extract_box()[1][0]
   ymax = lmp.extract_box()[1][1]
   zmax = lmp.extract_box()[1][2]

# initialize list

   iso_positions = []

# loop through all pairs of transition metal contacts

   for atom in pair_ids:

# compute distance components between atoms in pair and lower id coordinates

        lmp.commands_string(f'''
        group Atom1 id {int(atom[0])}
        group Atom2 id {int(atom[1])}
        compute Atom1X Atom1 property/atom x
        compute Atom1Y Atom1 property/atom y
        compute Atom1Z Atom1 property/atom z
        compute Atom2X Atom2 property/atom x
        compute Atom2Y Atom2 property/atom y
        compute Atom2Z Atom2 property/atom z      
        run 0
        ''')
# extract distance components and coordinates

        coord1X = [x for x in lmp.numpy.extract_compute('Atom1X', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if x != 0]
        coord1Y = [y for y in lmp.numpy.extract_compute('Atom1Y', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if y != 0]
        coord1Z = [z for z in lmp.numpy.extract_compute('Atom1Z', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if z != 0]
        coord2X = [x for x in lmp.numpy.extract_compute('Atom2X', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if x != 0]
        coord2Y = [y for y in lmp.numpy.extract_compute('Atom2Y', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if y != 0]
        coord2Z = [z for z in lmp.numpy.extract_compute('Atom2Z', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if z != 0]

# list comprehension above removes zeros, add zeros back if coordinate happens to be zero

        if coord1X == []:
            coord1X = [0.0]
        if coord1Y == []:
            coord1Y = [0.0]
        if coord1Z == []:
            coord1Z = [0.0]
        if coord2X == []:
            coord2X = [0.0]
        if coord2Y == []:
            coord2Y = [0.0]
        if coord2Z == []:
            coord2Z = [0.0]

# delete computes and groups before next iteration

        lmp.commands_string('''
        uncompute Atom1X
        uncompute Atom1Y
        uncompute Atom1Z
        uncompute Atom2X
        uncompute Atom2Y
        uncompute Atom2Z
        group Atom1 delete
        group Atom2 delete
        ''')

# compute coordinates as an angle on a circle

        angle1X = coord1X[0]/xmax*2*np.pi
        angle1Y = coord1Y[0]/ymax*2*np.pi
        angle1Z = coord1Z[0]/zmax*2*np.pi
        angle2X = coord2X[0]/xmax*2*np.pi
        angle2Y = coord2Y[0]/ymax*2*np.pi
        angle2Z = coord2Z[0]/zmax*2*np.pi

# find average position on circle

        x_p1_avg = np.cos(angle1X) + np.cos(angle2X)
        x_p2_avg = np.sin(angle1X) + np.sin(angle2X)
        y_p1_avg = np.cos(angle1Y) + np.cos(angle2Y)
        y_p2_avg = np.sin(angle1Y) + np.sin(angle2Y)
        z_p1_avg = np.cos(angle1Z) + np.cos(angle2Z)
        z_p2_avg = np.sin(angle1Z) + np.sin(angle2Z)

# find average angle from average positions

        x_ang_avg = np.arctan2(-x_p2_avg, -x_p1_avg) + np.pi
        y_ang_avg = np.arctan2(-y_p2_avg, -y_p1_avg) + np.pi
        z_ang_avg = np.arctan2(-z_p2_avg, -z_p1_avg) + np.pi

# find center of atoms

        x_center = xmax*x_ang_avg/(2*np.pi)
        y_center = ymax*y_ang_avg/(2*np.pi)
        z_center = zmax*z_ang_avg/(2*np.pi)

# append center coords to list

        iso_positions.append([x_center, y_center, z_center])
        
   return iso_positions
