from lammps import LMP_TYPE_ARRAY, LMP_STYLE_LOCAL, LMP_TYPE_VECTOR
import numpy as np

# function for finding all transition metal pairs within certain radius

def find_pairs(pairs, lmp):
   
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
      if distance > 3.1 and distance < 3.3:
         pair_ids.append(list(ids[i]))

# if a pair already had an isolobal atom placed in between then remove that pair from list
# the pairs list is a list of ids that have an isolobal atom 

   pair_ids = [ids for ids in pair_ids if ids not in pairs]
   return pair_ids