import subprocess as sp
import numpy as np

class CreateVasp:

#---------------------------- INIT ---------------------------#

    def __init__(self, compound_name, lammps_data):

        # initialize class attributes that are provided
        # when class is created 

        self.compound = compound_name

        with open(lammps_data, 'r') as lammps_file:
            self.lines = lammps_file.readlines()

        for i, line in enumerate(self.lines):
            if i == 0:
                elements = line.strip().split(' ')
                del elements[0]
                break
        
        self.elements = elements

#---------------------------- ATTRIBUTES ---------------------------#

    # list of all elements to check if lammps_data file has valid elements in header
    
    periodic_table = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
                    'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
                    'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
                    'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
                    'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
                    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
                    'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

#---------------------------- METHODS ---------------------------#

    def make_incar(self, algo):

        if algo != 'Normal' and algo != 'Exact' and algo != 'Fast' and algo != 'VeryFast':
            raise SystemExit('ALGO must be set to one of the following: Normal, Fast, VeryFast, or Exact')

        incar = open('INCAR', 'x')

        print('''SYSTEM = {}

    ISMEAR = 1         ! Methfessel-Paxton smearing
    SIGMA  = 0.2       ! smearing in eV
    ALGO   = {}
    PREC   = Accurate  ! precision
    ISYM   = -1        ! no symmetry imposed
    IBRION = -1        ! no relaxtion 
    LREAL  = A         ! automatic optimization of real space projectors
    EDIFF  = 1E-6      ! convergence criteria, difference in energy
    EDIFFG = -1E-5     ! convergence criteria, differnce in force
    LWAVE  = .FALSE.   ! no WAVECAR printed
    LCHARG = .FALSE.   ! no CHGCAR printed
    LAECHG = .FALSE.   ! all-electron density not printed
    '''.format(self.compound, algo), file = incar)
        
        incar.close()
        
        return

#--------------------------------------------------------------#

    def make_kpoints(self):

        kpoints = open('KPOINTS', 'x')

        print('''Gamma point only
    0
    Gamma 
    1 1 1
    0 0 0 
    ''', file = kpoints)
        
        kpoints.close()

        return

#--------------------------------------------------------------#

    def make_poscar(self):

        # create dictionary with key for each element and list for each value
        
        atom_coord = {element: [] for element in self.elements}

        for i, line in enumerate(self.lines):

            # find the box lengths

            if i == 5:
                x_max = line.strip().split(' ')[1]
                x_min = float(line.strip().split(' ')[0])

            if i == 6:
                y_max = line.strip().split(' ')[1]
                y_min = float(line.strip().split(' ')[0])

            if i == 7:
                z_max = line.strip().split(' ')[1]
                z_min = float(line.strip().split(' ')[0])

            # append coordinates to dictionary 
            # lammps_data file stores atom coordinates as:
            # "atom index" "atom type" "atom_x" "atom_y" "atom_z"
            # get element name from atom type by listing dictionary keys
            # if atom type is greater than total number of elements then
            # it is an isolobal atom and will not be in the poscar
            # remove "atom index" and "atom type"
            # append x, y, z to value list in dictionary

            if i > 10:
                coord = line.strip().split(' ')
                element_type = int(coord[1]) - 1
                if int(coord[1]) > len(self.elements):
                    continue
                del coord[0:2]
                element = list(atom_coord)[element_type]
                coord = np.array(list(map(float, coord)))
                org_shift = np.array([x_min, y_min, z_min])
                coord = list(coord - org_shift)
                atom_coord[element].append(coord)

        # check if lammps_data file has correct header
            
        for atom in self.elements:
            if atom not in CreateVasp.periodic_table:
                raise SystemExit('First line of lammps.data file must only contain a "#" and a list of elements\ni.e. For NiSi2 the first line would read: "# Ni Si" (without the quotations)')

        poscar = open('POSCAR', 'x')

        print('''{}                  
    1.00000000000000     
        {} 0.0 0.0
        0.0 {} 0.0
        0.0 0.0 {}      
    '''.format(self.compound, x_max, y_max, z_max), end = '', file = poscar)
        
        for atom in self.elements:
            print(atom, end = ' ', file = poscar)

        print('\n', end = '', file = poscar)

        for atom in self.elements:
            print(len(atom_coord[atom]), end = ' ', file = poscar)

        print('\nCartesian', file = poscar)

        for atom in atom_coord:
            for coord in atom_coord[atom]:
                print('  {} {} {}'.format(coord[0], coord[1], coord[2]), file = poscar)    

        poscar.close()   

        return
    
#--------------------------------------------------------------#

    def make_potcar(self, path_to_potentials):

        ptp = path_to_potentials.strip()

        print('Select desired potential (PBE/GGA)')
        
        while True:
            fxnal = input()
            if fxnal != 'PBE' and fxnal != 'pbe' and fxnal != 'GGA' and fxnal != 'gga':
                print('Invalid option, please select either PBE or GGA')
            else:
                break

        fxnal = str(fxnal).upper()

        potentials = []

        for atom in self.elements:

            check_semicore = sp.run([f'ls {ptp}/potpaw_{fxnal}/{atom}_d'], shell = True, capture_output = True).stdout.decode()
            if check_semicore != '':
                print(f'Would you like to use d semicore states for {atom}? (Y/N)')
                check = input()
                if check == 'Y' or check == 'y' or check == 'Yes' or check == 'YES' or check == 'yes':
                    potentials.append(str(atom) + '_d')
                    continue

            check_semicore = sp.run([f'ls {ptp}/potpaw_{fxnal}/{atom}_pv'], shell = True, capture_output = True).stdout.decode()
            if check_semicore != '':
                print(f'Would you like to use p semicore states for {atom}? (Y/N)')
                check = input()
                if check == 'Y' or check == 'y' or check == 'Yes' or check == 'YES' or check == 'yes':
                    potentials.append(str(atom) + '_pv')
                    continue

            check_semicore = sp.run([f'ls {ptp}/potpaw_{fxnal}/{atom}_sv'], shell = True, capture_output = True).stdout.decode()
            if check_semicore != '':
                print(f'Would you like to use s semicore states for {atom}? (Y/N)')
                check = input()
                if check == 'Y' or check == 'y' or check == 'Yes' or check == 'YES' or check == 'yes':
                    potentials.append(str(atom) + '_sv')
                    continue

            else:
                potentials.append(atom)

        zcat_str = ''

        for potential in potentials:
            zcat_str = zcat_str + f"{ptp}/potpaw_{fxnal}/{potential}/POTCAR.Z "

        sp.run([f'zcat {zcat_str}> POTCAR'], shell = True)

        return
    
#--------------------------------------------------------------#

    def lammps_data_edit(self, data_to_edit):

        data_file = open(data_to_edit, 'r+')

        data_lines = data_file.readlines()

        data_file.seek(0)
        data_file.truncate()

        data_lines[0] = '# '

        for atom in self.elements:
            data_lines[0] = data_lines[0] + f'{atom} '

        data_lines[0] = data_lines[0] + '\n'

        for i, line in enumerate(data_lines):
            if 'Masses' == line.strip():
                a = i
            if 'Atoms' == line.strip().split(' ')[0]:
                b = i
            if 'Velocities' == line.strip():
                c = i
                break

        d = len(data_lines)

        for ind in range(b + 2, c - 1):
            vals = data_lines[ind].strip().split(' ')
            del vals[5:8]
            data_lines[ind] = ' '.join(vals) + '\n'

        for ind in range(a,b):
            del data_lines[a]           
        
        for ind in range(c, d):
            del data_lines[c - (b - a)]

        for i, line in enumerate(data_lines):
            if i > 10 and line == '\n':
                continue
            else:
                data_file.write(line)

        data_file.close()

        return
    
#--------------------------------------------------------------#
