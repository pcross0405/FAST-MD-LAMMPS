import subprocess as sp
import py4vasp as p4v
from lammps import LMP_STYLE_ATOM, LMP_TYPE_VECTOR
import numpy as np
import time
import sys
from . import create_vasp
from . import fastmdtools as fmd

class ManageVasp:

    def __init__(self, working_dir, job_name, remote):
        self.working_dir = working_dir
        self.job_name = job_name
        self.remote = remote
        self.user = sp.run(['echo $USER'], shell = True, capture_output = True).stdout.decode().strip()

    #-----------------------------------------------------------------------------------------#
    #---------------------------------- FUNCTION DEFINTIONS ----------------------------------#
    #-----------------------------------------------------------------------------------------#

    # define functions to run while loop to check status of VASP job periodically

    def normal_check(self):

        oszicar_check = self.working_dir + '/' + 'OSZICAR'

        while True:
            
            status = sp.run([f'ssh {self.user}@{self.remote} "qstat | grep {self.job_name}"'], shell = True, capture_output = True).stdout.decode().split(' ')
            if status != ['']:
                break

        while True:

            # grep to find status of job
            # remove spaces from status list

            status = sp.run([f'ssh {self.user}@{self.remote} "qstat | grep {self.job_name}"'], shell = True, capture_output = True).stdout.decode().split(' ')
            status = [x for x in status if x != '']
            job_number = [x for x in status[0] if x == 0 or x == 1 or x == 2 or x == 3 or x == 4 or x == 5 or x == 6 or x == 7 or x == 8 or x == 9]
            job_number = ''.join(job_number)

            # check if max number of electronic iterations has been reached

            if sp.run([f'grep "DAV:  60" {oszicar_check}'], shell = True, capture_output = True).stdout.decode() != '':
                print('Maximum number of eletronic iterations has been reached, cancelling job, exiting script')
                sp.run([f'ssh {self.user}@{self.remote} "qdel {job_number}"'], shell = True)
                sys.exit()
            
            # check if job has completed
            # break if it has, or wait 30 sec and check again

            if status[-3] == 'C':

                break

            else:

                time.sleep(30)

        return

    #------------------------------------------------------------------------------------------#

    def exact_check(self):

        oszicar_check = self.working_dir + '/' + 'OSZICAR'

        while True:
            
            status = sp.run([f'ssh {self.user}@{self.remote} "qstat | grep {self.job_name}"'], shell = True, capture_output = True).stdout.decode().split(' ')
            if status != ['']:
                break

        while True:

            # grep to find status of job
            # remove spaces from status list

            status = sp.run([f'ssh {self.user}@{self.remote} "qstat | grep {self.job_name}"'], shell = True, capture_output = True).stdout.decode().split(' ')
            status = [x for x in status if x != '']
            job_number = [x for x in status[0] if x == 0 or x == 1 or x == 2 or x == 3 or x == 4 or x == 5 or x == 6 or x == 7 or x == 8 or x == 9]
            job_number = ''.join(job_number)

            # check if max number of electronic iterations has been reached

            if sp.run([f'grep "DIA:  60" {oszicar_check}'], shell = True, capture_output = True).stdout.decode() != '':
                print('Maximum number of eletronic iterations has been reached, cancelling job, exiting script')
                sp.run([f'ssh {self.user}@{self.remote} "qdel {job_number}"'], shell = True)
                sys.exit()
            
            # check if job has completed
            # return if it has, or wait 30 sec and check again

            if status[-3] == 'C':

                break

            else:

                time.sleep(30)

        return

    #------------------------------------------------------------------------------------------#

    def fast_check(self):

        oszicar_check = self.working_dir + '/' + 'OSZICAR'

        while True:

            status = sp.run([f'ssh {self.user}@{self.remote} "qstat | grep {self.job_name}"'], shell = True, capture_output = True).stdout.decode().split(' ')
            if status != ['']:
                break

        while True:

            # grep to find status of job
            # remove spaces from status list

            status = sp.run([f'ssh {self.user}@{self.remote} "qstat | grep {self.job_name}"'], shell = True, capture_output = True).stdout.decode().split(' ')
            status = [x for x in status if x != '']
            job_number = [x for x in status[0] if x == 0 or x == 1 or x == 2 or x == 3 or x == 4 or x == 5 or x == 6 or x == 7 or x == 8 or x == 9]
            job_number = ''.join(job_number)

            # check if max number of electronic iterations has been reached

            if sp.run([f'grep "RMM:  60" {oszicar_check}'], shell = True, capture_output = True).stdout.decode() != '':
                print('Maximum number of eletronic iterations has been reached, cancelling job, exiting script')
                sp.run([f'ssh {self.user}@{self.remote} "qdel {job_number}"'], shell = True)
                sys.exit()

            # check if job has completed
            # return if it has, or wait 30 sec and check again

            if status[-3] == 'C':

                break

            else:

                time.sleep(30)

        return

    #------------------------------------------------------------------------------------------#

    def veryfast_check(self):

        oszicar_check = self.working_dir + '/' + 'OSZICAR'

        while True:

            status = sp.run([f'ssh {self.user}@{self.remote} "qstat | grep {self.job_name}"'], shell = True, capture_output = True).stdout.decode().split(' ')
            if status != ['']:
                break

        while True:

            # grep to find status of job
            # remove spaces from status list

            status = sp.run([f'ssh {self.user}@{self.remote} "qstat | grep {self.job_name}"'], shell = True, capture_output = True).stdout.decode().split(' ')
            status = [x for x in status if x != '']
            job_number = [x for x in status[0] if x == 0 or x == 1 or x == 2 or x == 3 or x == 4 or x == 5 or x == 6 or x == 7 or x == 8 or x == 9]
            job_number = ''.join(job_number)

            # check if max number of electronic iterations has been reached

            if sp.run([f'grep "RMM:  60" {oszicar_check}'], shell = True, capture_output = True).stdout.decode() != '':
                print('Maximum number of eletronic iterations has been reached, cancelling job, exiting script')
                sp.run([f'ssh {self.user}@{self.remote} "qdel {job_number}"'], shell = True)
                sys.exit()

            # check if job has completed
            # return if it has, or wait 30 sec and check again

            if status[-3] == 'C':

                break

            else:

                time.sleep(30)

        return

    #------------------------------------------------------------------------------------------#

    def check_vasp(self):

        incar_check = self.working_dir + '/' + 'INCAR'

        algo = sp.run([f'grep ALGO {incar_check}'], shell = True, capture_output = True).stdout.decode().split(' ')
        algo = algo[-1]

        # check if algorithm is Normal or Exact
        # algo determines how to check for convergence

        if algo == 'Normal\n':

            self.normal_check()
                
            # give some time to ensure output files are completely printed out

            time.sleep(30)

        elif algo == 'Exact\n':

            self.exact_check()

            # give some time to ensure output files are completely printed out

            time.sleep(30)

        elif algo == 'Fast\n':

            self.fast_check()

            # give some time to ensure output files are completely printed out

            time.sleep(30)

        elif algo == 'VeryFast\n':

            self.veryfast_check()

            # give some time to ensure output files are completely printed out

            time.sleep(30)

        else:

            print('''ALGO in INCAR must be set to one of Normal, Fast, VeryFast, or Exact
        ---------------------------------------------------
        DO NOT INCLUDE ANY SPACES AFTER ALGO SETTING
        ''')
            sys.exit()
        
        return
    
    #------------------------------------------------------------------------------------------#
    
    def run_vasp(self, lmp, job_count, job_subcount):

        # initialize job_dir variable to cd into and submit VASP job in

        current_dir = sp.run(['pwd'], shell = True, capture_output = True).stdout.decode().strip()
        job_dir = current_dir + '/' + self.working_dir

        sp.run([f'ssh {self.user}@{self.remote} "cd {job_dir} ; qvasp6 -np 12 fast_md{job_count}_{job_subcount}" > /dev/null'], shell = True)

        # after VASP job has been submitted, check its status

        new_job_check = ManageVasp(self.working_dir, f'fast_md{job_count}_{job_subcount}', self.remote)
        new_job_check.check_vasp()

        # once VASP has finished, process results with py4vasp
        # fetch forces and stress on cell to use for updating LAMMPS

        calc = p4v.Calculation.from_path(f'{self.working_dir}')
        forces = calc.force.to_dict()['forces']

        # find total number of LAMMPS atom types from lammps.data file

        num_types = sp.run(['grep "atom types" *data'], shell = True, capture_output = True).stdout.decode().strip().split(' ')
        num_types = int(num_types[0]) - 1

        # create LAMMPS group and compute commands for each atom type

        group_commands = []
        compute_commands = []
        group_del_commands = []
        uncompute_commands = []
        for i in range(num_types):
            group_commands.append(f'group Type{i + 1}Atoms type {i + 1}')
            compute_commands.append(f'compute Type{i + 1}IDs Type{i + 1}Atoms property/atom id')
            group_del_commands.append(f'group Type{i + 1}Atoms delete')
            uncompute_commands.append(f'uncompute Type{i + 1}IDs')

        # run LAMMPS group and compute commands

        for i in range(num_types):
            lmp.commands_string('''
            {}
            {}                                
            '''.format(group_commands[i], compute_commands[i]))

        # create dictionary for extracting
        # IDs from compute commands to
        # remove any 0 enteries         

        type_ids = {}
        all_ids = []
        for i in range(num_types):
            key = f'type_{i + 1}_ids'
            type_ids[key] = lmp.numpy.extract_compute(f'Type{i + 1}IDs', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64)
            type_ids[key] = [x for x in type_ids[key] if x != 0]
            for j in type_ids[key]:
                all_ids.append(j)

        # initialize dictionary with IDs each matching to an empty list
        # list will be updated with cooresponding atom's coordinates

        id_and_coord = {id_val : [] for id_val in all_ids}

        # run LAMMPS group delete and uncompute commands

        for i in range(num_types):
            lmp.commands_string('''
            {}
            {}                                
            '''.format(uncompute_commands[i], group_del_commands[i]))

        # get coordinate of each atom from LAMMPS using compute
        # add VASP forces with LAMMPS addforce to small spherical
        # zone centered about each atom   

        count = 0
        total_sum = 0
        for type_val in type_ids:
            if count > 0:
                total_sum = total_sum + len(type_ids[f'type_{count}_ids'])
            for i, atom_id in enumerate(type_ids[type_val]):
                lmp.commands_string('''
                group AtomToForce id {}
                compute ForceAtomX AtomToForce property/atom x
                compute ForceAtomY AtomToForce property/atom y
                compute ForceAtomZ AtomToForce property/atom z
                '''.format(int(atom_id)))
                coord1X = [x for x in lmp.numpy.extract_compute('ForceAtomX', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if x != 0]
                coord1Y = [y for y in lmp.numpy.extract_compute('ForceAtomY', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if y != 0]
                coord1Z = [z for z in lmp.numpy.extract_compute('ForceAtomZ', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if z != 0]
                if coord1X == []:
                    coord1X = [0.0]
                if coord1Y == []:
                    coord1Y = [0.0]
                if coord1Z == []:
                    coord1Z = [0.0]
                id_and_coord[atom_id] = [coord1X[0], coord1Y[0], coord1Z[0]]
                lmp.commands_string('''
                uncompute ForceAtomX
                uncompute ForceAtomY
                uncompute ForceAtomZ
                group AtomToForce delete
                region ForceAtom{} sphere {} {} {} 0.5
                fix Atom{} all addforce {} {} {} region ForceAtom{}
                '''.format(int(i + total_sum), coord1X[0], coord1Y[0], coord1Z[0], int(i + total_sum), forces[i + total_sum][0], forces[i + total_sum][1], forces[i + total_sum][2], int(i + total_sum)))
            count += 1

        return  type_ids, id_and_coord

    #------------------------------------------------------------------------------------------#

    def run_check(self, lmp, pairs_w_iso_atom, added_atoms, iso_limit, min_nrg, new_id, job_count, job_subcount, min_style = '', run_box_relax = ''):

        # generate POSCAR from current geometry

        vasp_poscar = create_vasp.CreateVasp('NiSi2', f'{job_count}_{job_subcount}_geo')
        vasp_poscar.make_poscar()
        sp.run([f'mv POSCAR {self.working_dir}'], shell = True)
        sp.run([f'mv {job_count}_{job_subcount}_geo geo_dir/'], shell = True)

        # run VASP calc

        type_ids, id_and_coord_1 = self.run_vasp(lmp, job_count, job_subcount)

        # get stress tensor from VASP calc

        calc = p4v.Calculation.from_path(f'{self.working_dir}')
        stress = calc.stress.to_dict()['stress']
        stress_vals = [stress[0][0], stress[1][1], stress[2][2], stress[1][0], stress[2][0], stress[2][1]]

        # minimize current ionic positions

        lmp.commands_string('''
        min_style {}
        min_modify dmax 0.025
        minimize 0 1e-10 100 10000
        '''.format(min_style))

        # check if isolobal atoms remained centered

        fmd.check_center(pairs_w_iso_atom, added_atoms, lmp, job_count, min_nrg, iso_limit, new_id)

        # find box lengths

        x_len1 = lmp.extract_box()[1][0] - lmp.extract_box()[0][0]
        y_len1 = lmp.extract_box()[1][1] - lmp.extract_box()[0][1]
        z_len1 = lmp.extract_box()[1][2] - lmp.extract_box()[0][2]

        # check to see if box_relax is requested

        if run_box_relax == '':
            pass

        # check for triclinic box

        elif sp.run(['grep xz *data'], shell = True, capture_output = True).stdout.decode() != '': 

            # relax box including tilt factors

            fmd.box_relax(lmp, stress_vals[0]*1000, stress_vals[1]*1000, stress_vals[2]*1000, stress_vals[3]*1000, stress_vals[4]*1000, stress_vals[5]*1000)

        else:

            # relax box including only x, y, z components of stress tensor

            fmd.box_relax(lmp, stress_vals[0]*1000, stress_vals[1]*1000, stress_vals[2]*1000)

        # find new box lengths

        x_len2 = lmp.extract_box()[1][0] - lmp.extract_box()[0][0]
        y_len2 = lmp.extract_box()[1][1] - lmp.extract_box()[0][1]
        z_len2 = lmp.extract_box()[1][2] - lmp.extract_box()[0][2]

        # create new id and coordinate dictionary to compare against

        id_and_coord_2 = dict(id_and_coord_1)

        # find new coordinates after minimization

        for atom_id in id_and_coord_1:
            lmp.commands_string('''
            group AtomToForce id {}
            compute ForceAtomX AtomToForce property/atom x
            compute ForceAtomY AtomToForce property/atom y
            compute ForceAtomZ AtomToForce property/atom z
            '''.format(int(atom_id)))
            coord1X = [x for x in lmp.numpy.extract_compute('ForceAtomX', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if x != 0]
            coord1Y = [y for y in lmp.numpy.extract_compute('ForceAtomY', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if y != 0]
            coord1Z = [z for z in lmp.numpy.extract_compute('ForceAtomZ', LMP_STYLE_ATOM, LMP_TYPE_VECTOR).astype(np.float64) if z != 0]
            if coord1X == []:
                coord1X = [0.0]
            if coord1Y == []:
                coord1Y = [0.0]
            if coord1Z == []:
                coord1Z = [0.0]
            id_and_coord_2[atom_id] = [coord1X[0], coord1Y[0], coord1Z[0]]
            lmp.commands_string('''
            uncompute ForceAtomX
            uncompute ForceAtomY
            uncompute ForceAtomZ
            group AtomToForce delete
            ''')

        # compare how much atoms have moved between minimizations

        for atom_id in id_and_coord_1:

            # calculate relative positions of each atom before an after minimization

            rel_pos1 = [id_and_coord_1[atom_id][0]/x_len1, id_and_coord_1[atom_id][1]/y_len1, id_and_coord_1[atom_id][2]/z_len1]
            rel_pos2 = [id_and_coord_2[atom_id][0]/x_len2, id_and_coord_2[atom_id][1]/y_len2, id_and_coord_2[atom_id][2]/z_len2]

            # find difference between relative positions in angular units

            x_ang_diff = rel_pos1[0]*2*np.pi - rel_pos2[0]*2*np.pi
            y_ang_diff = rel_pos1[1]*2*np.pi - rel_pos2[1]*2*np.pi
            z_ang_diff = rel_pos1[2]*2*np.pi - rel_pos2[2]*2*np.pi

            # find total difference in linear units by converting back to linear from angular

            x_diff = np.arctan2(np.sin(x_ang_diff), np.cos(x_ang_diff))/(2*np.pi)
            y_diff = np.arctan2(np.sin(y_ang_diff), np.cos(y_ang_diff))/(2*np.pi)
            z_diff = np.arctan2(np.sin(z_ang_diff), np.cos(z_ang_diff))/(2*np.pi)

            # if relative difference in any direction is greater than 0.1%, rerun

            if x_diff > 0.001 or y_diff > 0.001 or z_diff > 0.001:

                # update sub-job count (this is number of reruns)

                job_subcount += 1

                # initialize variables for clearing current LAMMPS group and computes

                count = 0
                total_sum = 0

                # loop over id dictionary

                for type_val in type_ids:
                    if count > 0:
                        total_sum = total_sum + len(type_ids[f'type_{count}_ids'])
                    for i, atom_id in enumerate(type_ids[type_val]):
                        lmp.commands_string('''
                        unfix Atom{}
                        region ForceAtom{} delete
                        '''.format(int(i + total_sum), int(i + total_sum)))
                    count += 1

                # clear VASP directory before rerunning

                rmvasp = 'rm CHG CHGCAR CONTCAR *log DOSCAR EIGENVAL IBZKPT OSZICAR OUTCAR PCDAT WAVECAR XDATCAR vasprun.xml vaspout.h5 REPORT'                
                sp.run([f'cd {self.working_dir} ; {rmvasp} ; cd ../'], shell = True)

                # generate lammps data file with current geometry

                lmp.command(f'write_data {job_count}_{job_subcount}_geo')
                time.sleep(3)

                # edit the lammps write data file to usable format
                
                vasp_poscar.lammps_data_edit(f'{job_count}_{job_subcount}_geo')

                # rerun

                self.run_check(lmp, pairs_w_iso_atom, added_atoms, iso_limit, min_nrg, new_id, job_count, job_subcount, min_style)
                break

        return type_ids
    
    #-----------------------------------------------------------------------------------------#
    #-------------------------------- FUNCTION DEFINTIONS END --------------------------------#
    #-----------------------------------------------------------------------------------------#
