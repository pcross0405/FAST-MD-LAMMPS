#----------------------------------------------------------------------------------------------------------------#
#------------------ FAST-MD module containing all functions used in the FAST-MD implementation ------------------#
#----------------------------------------------------------------------------------------------------------------#

import numpy as np
from lammps import LMP_TYPE_VECTOR, LMP_STYLE_ATOM, LMP_TYPE_ARRAY, LMP_STYLE_LOCAL, LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR, LMP_VAR_EQUAL

#----------------------------------------------------------------------------------------------------------------#

# function for optimizing simulation box

def box_relax(lmp, x_press, y_press, z_press, xy_press = None, xz_press = None, yz_press = None):

    if xy_press != None and xz_press != None and yz_press != None:

        # define fix to optimize box during minimization run

        lmp.commands_string(f'''
        fix BoxRelax all box/relax x {x_press} y {y_press} z {z_press} xy {xy_press} xz {xz_press} yz {yz_press} nreset 1
        min_style sd
        minimize 0 1e-10 1000 10000
        unfix BoxRelax
        min_style sd
        minimize 0 1e-10 1000 10000
        ''')

    else:

        # define fix to optimize box during minimization run

        lmp.commands_string(f'''
        fix BoxRelax all box/relax x {x_press} y {y_press} z {z_press} nreset 1
        min_style sd
        minimize 0 1e-10 1000 10000
        unfix BoxRelax
        min_style sd
        minimize 0 1e-10 1000 10000
        ''')

    # extract box lengths after first optimization

    xmax1 = lmp.extract_box()[1][0]
    ymax1 = lmp.extract_box()[1][1]
    zmax1 = lmp.extract_box()[1][2]

    xmin1 = lmp.extract_box()[0][0]
    ymin1 = lmp.extract_box()[0][1]
    zmin1 = lmp.extract_box()[0][2]

    if xy_press != None and xz_press != None and yz_press != None:

        # define fix to optimize box during minimization run

        lmp.commands_string(f'''
        fix BoxRelax all box/relax x {x_press} y {y_press} z {z_press} xy {xy_press} xz {xz_press} yz {yz_press} nreset 1
        min_style sd
        minimize 0 1e-10 1000 10000
        unfix BoxRelax
        min_style sd
        minimize 0 1e-10 1000 10000
        ''')

    else:

        # define fix to optimize box during minimization run

        lmp.commands_string(f'''
        fix BoxRelax all box/relax x {x_press} y {y_press} z {z_press} nreset 1
        min_style sd
        minimize 0 1e-10 1000 10000
        unfix BoxRelax
        min_style sd
        minimize 0 1e-10 1000 10000
        ''')

    # extract box lengths after second optimization

    xmax2 = lmp.extract_box()[1][0]
    ymax2 = lmp.extract_box()[1][1]
    zmax2 = lmp.extract_box()[1][2]

    xmin2 = lmp.extract_box()[0][0]
    ymin2 = lmp.extract_box()[0][1]
    zmin2 = lmp.extract_box()[0][2]

    # if box lengths are not within 0.01%, rerun

    if xmin1/xmin2 < 0.9999 or xmin1/xmin2 > 1.0001:
        box_relax(lmp, x_press, y_press, z_press, xy_press, xz_press, yz_press)
    if ymin1/ymin2 < 0.9999 or ymin1/ymin2 > 1.0001:
        box_relax(lmp, x_press, y_press, z_press, xy_press, xz_press, yz_press)
    if zmin1/zmin2 < 0.9999 or zmin1/zmin2 > 1.0001:
        box_relax(lmp, x_press, y_press, z_press, xy_press, xz_press, yz_press)
    if xmax1/xmax2 < 0.9999 or xmax1/xmax2 > 1.0001:
        box_relax(lmp, x_press, y_press, z_press, xy_press, xz_press, yz_press)
    if ymax1/ymax2 < 0.9999 or ymax1/ymax2 > 1.0001:
        box_relax(lmp, x_press, y_press, z_press, xy_press, xz_press, yz_press)
    if zmax1/zmax2 < 0.9999 or zmax1/zmax2 > 1.0001:
        box_relax(lmp, x_press, y_press, z_press, xy_press, xz_press, yz_press)

    # if box lengths are within 0.01%, return

    else:
        return

#----------------------------------------------------------------------------------------------------------------#

# function for checking if isolobal atoms remain in center of transition metals

def check_center(tm_pairs, added_atoms, lmp, n, min_nrg, iso_limit, new_id, min_style = 'sd'):
   
    # find current centers of pairs of transition metals with isolobal atom in between

    current_centers = find_centers(tm_pairs, lmp)

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

        iso_atom_x = [x for x in lmp.numpy.extract_compute('IsolobalAtomX', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if x != 0]
        iso_atom_y = [y for y in lmp.numpy.extract_compute('IsolobalAtomY', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if y != 0]
        iso_atom_z = [z for z in lmp.numpy.extract_compute('IsolobalAtomZ', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if z != 0]

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

        dx = center[0] - iso_atom_x[0]
        dy = center[1] - iso_atom_y[0]
        dz = center[2] - iso_atom_z[0]

        # move isolobal atom back to true center
        # minimize after moving atom 

        lmp.commands_string(f'''
        group MoveAtom id {added_atoms[i]}
        displace_atoms MoveAtom move {dx} {dy} {dz}
        min_style {min_style}
        minimize 0 1e-10 1000 10000
        group MoveAtom delete
        ''')

        # recompute distance from center after minimization
        # compute x,y,z coordinates of isolobal atom in between transition metal pair

        lmp.commands_string(f'''
        group IsolobalAtom id {added_atoms[i]}
        compute IsolobalAtomX IsolobalAtom property/atom x
        compute IsolobalAtomY IsolobalAtom property/atom y
        compute IsolobalAtomZ IsolobalAtom property/atom z
        run 0
        ''')

        # extract coordinates of isolobal atom

        iso_atom_x = [x for x in lmp.numpy.extract_compute('IsolobalAtomX', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if x != 0]
        iso_atom_y = [y for y in lmp.numpy.extract_compute('IsolobalAtomY', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if y != 0]
        iso_atom_z = [z for z in lmp.numpy.extract_compute('IsolobalAtomZ', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if z != 0]

        # iso_atom_x looks like [0, 0, 0, ...., X, ....., 0, 0, 0] where 0's are put in for all atoms not part of the lammps compute
        # list comprehension above removes 0's, add 0 back if coordinate happens to be zero

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

        dx = center[0] - iso_atom_x[0]
        dy = center[1] - iso_atom_y[0]
        dz = center[2] - iso_atom_z[0]

        # compute distance between atom and true center

        check_distance = ((dx)**2 + (dy)**2 + (dz)**2)**(1/2)

        # check if distance is greater than cutoff

        if check_distance > 0.01:

            # get current step

            step = lmp.get_thermo('step')

            # update progress bar

            progress_bar(n, iso_limit, int(step), min_nrg, new_id, prefix = 'Progress', suffix = '')

            # rerun function to see if atom stayed centered after minimization

            check_center(tm_pairs, added_atoms, lmp, n, min_nrg, iso_limit, new_id)
    
    return

#----------------------------------------------------------------------------------------------------------------#

# function for finding all positions at the center of transition metal contacts

def find_centers(pair_ids, lmp):

    # grab lengths of simulation box

    xmax = lmp.extract_box()[1][0] - lmp.extract_box()[0][0]
    ymax = lmp.extract_box()[1][1] - lmp.extract_box()[0][1]
    zmax = lmp.extract_box()[1][2] - lmp.extract_box()[0][2]

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

        x_p1_avg = (np.cos(angle1X) + np.cos(angle2X))/2
        x_p2_avg = (np.sin(angle1X) + np.sin(angle2X))/2
        y_p1_avg = (np.cos(angle1Y) + np.cos(angle2Y))/2
        y_p2_avg = (np.sin(angle1Y) + np.sin(angle2Y))/2
        z_p1_avg = (np.cos(angle1Z) + np.cos(angle2Z))/2
        z_p2_avg = (np.sin(angle1Z) + np.sin(angle2Z))/2

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

#----------------------------------------------------------------------------------------------------------------#

# function for finding all transition metal pairs within certain radius

def find_pairs(pairs, min_dist, max_dist, lmp):
   
    # find all transition metal pairs

    lmp.commands_string('''
    group               TransitionMetals type 1 
    compute             TransitionMetalDists TransitionMetals pair/local dist
    compute             atomIDs TransitionMetals property/local patom1 patom2
    run                 0
    ''')

    # extract distances to dists and atom ids to ids

    dists = lmp.numpy.extract_compute('TransitionMetalDists', LMP_STYLE_LOCAL, LMP_TYPE_VECTOR).astype(np.float64)
    ids = lmp.numpy.extract_compute('atomIDs', LMP_STYLE_LOCAL, LMP_TYPE_ARRAY).astype(np.float64)

    # delete computes and group before calling function again

    lmp.commands_string('''
    uncompute TransitionMetalDists
    uncompute atomIDs
    group TransitionMetals delete
    ''')

    # add all transition metal atom pairs within first coordination sphere to new list

    pair_ids = []
    for i, distance in enumerate(dists):
        if distance > min_dist and distance < max_dist:
            pair_ids.append(list(ids[i]))

    # if a pair already had an isolobal atom placed in between then remove that pair from list
    # the pairs list is a list of ids that have an isolobal atom 

    pair_ids = [ids for ids in pair_ids if ids not in pairs]
    return pair_ids

#----------------------------------------------------------------------------------------------------------------#

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

#----------------------------------------------------------------------------------------------------------------#

# fxn for looping through isolobal postions

def iso_loop(new_id, positions, lmp, n, iso_limit):
   
    # initialize minimum energy variable

    min_nrg = 10**10

    # extract simulation box lengths

    x_len = lmp.extract_box()[1][0] - lmp.extract_box()[0][0]
    y_len = lmp.extract_box()[1][1] - lmp.extract_box()[0][1]
    z_len = lmp.extract_box()[1][2] - lmp.extract_box()[0][2]
    box_lengths = [x_len, y_len, z_len]

    # loop through all possible isolobal sites

    for atom_pos in positions:
      
        # check if an isolobal atom has already been added to a site
        # define lammps regions to check if an atom overlaps with region
        # region1 will have coordinates equal to the to-be-added atom's position
        # region2 will be the periodic image of region if region1 happens to be on boundary

        region2 = list(atom_pos)

        for i, val in enumerate(region2):
            region2[i] = val/box_lengths[i]*2*np.pi

        for i, val in enumerate(region2):
            if val > (box_lengths[i] - 0.5)/box_lengths[i]*2*np.pi:
                region2[i] = val - 2*np.pi

            elif val < (0.5)/box_lengths[i]*2*np.pi:
                region2[i] = val + 2*np.pi
                
        for i, val in enumerate(region2):
            region2[i] = box_lengths[i]*(np.arctan2(-np.sin(region2[i]),-np.cos(region2[i])) + np.pi)/(2*np.pi)

        # make lammps regions and groups to check if atom overlaps in region

        lmp.commands_string(f'''
        region AtomSphere1 sphere {atom_pos[0]} {atom_pos[1]} {atom_pos[2]} 0.5
        region AtomSphere2 sphere {region2[0]} {region2[1]} {region2[2]} 0.5
        group AtomCheck1 region AtomSphere1
        group AtomCheck2 region AtomSphere2
        variable AtomCount1 equal count(AtomCheck1,AtomSphere1)
        variable AtomCount2 equal count(AtomCheck2,AtomSphere2)
        run 0
        ''')

        # extract number of atoms in region

        atom_count1 = lmp.extract_variable('AtomCount1', None, LMP_VAR_EQUAL)
        atom_count2 = lmp.extract_variable('AtomCount2', None, LMP_VAR_EQUAL)

        # delete groups and regions before next iteration

        lmp.commands_string('''
        group AtomCheck1 delete
        group AtomCheck2 delete
        region AtomSphere1 delete
        region AtomSphere2 delete
        variable AtomCount1 delete
        variable AtomCount2 delete                                                    
        ''')

        # if any atoms exist in region then move to next iteration

        if atom_count1 + atom_count2 > 0:
            continue
      
        # otherwise test isolobal atom 

        else:
            
            # find total pe at each site

            nrg = iso_add(atom_pos[0], atom_pos[1], atom_pos[2], new_id, lmp)

            # get current step for progress bar

            step = lmp.get_thermo('step')

            # update progress bar

            progress_bar(n, iso_limit, int(step), nrg, new_id, prefix = 'Progress', suffix = '')
            
            # find site with lowest potential energy after addition of isolobal atom

            if nrg < min_nrg:
                min_nrg = nrg

                # save coordinates of lowest energy isolobal atom

                coord = atom_pos

            # delete atom created by iso_add before testing next isolobal atom

            lmp.commands_string(f'''
            group isolobal_atom id {new_id}
            delete_atoms group isolobal_atom compress no
            ''')

    # return coordinates of lowest energy isolobal atom
    # energy is for updating progress bar later

    return [coord, min_nrg]

#----------------------------------------------------------------------------------------------------------------#

# progress bar code credit: https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters

def progress_bar(iteration, total, step, nrg, atoms, prefix = '', suffix = '', length = 100, fill = '█', up = '\x1B[3A', clr = '\x1B[0K'):
    percent = ('{0:.' + '1' + 'f}').format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'{up}{prefix} |{bar}| {percent}% {suffix}{clr}\nStep: {step}, Potential Energy: {nrg} eV, Atoms: {atoms}{clr}\n')

    # print new line when complete

    if iteration == total: 
        print('\nCOMPLETE!')
