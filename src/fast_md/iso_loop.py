import numpy as np
from lammps import LMP_VAR_EQUAL
from . import iso_add, progress_bar

# fxn for looping through isolobal postions

def iso_loop(new_id, positions, lmp, n, iso_limit):
   
# initialize minimum energy variable

   min_nrg = 10**10

# extract simulation box lengths

   box_lengths = [lmp.extract_box()[1][0], lmp.extract_box()[1][1], lmp.extract_box()[1][2]]

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

         progress_bar(n, iso_limit/2, int(step), nrg, new_id, prefix = 'Progress', suffix = '')
         
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
