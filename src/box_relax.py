# function for optimizing simulation box

def box_relax(lmp):
   
# define fix to optimize box during minimization run
   lmp.commands_string('''
   fix BoxRelax all box/relax aniso 0.0 vmax 0.001
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

# define second fix to optimize box

   lmp.commands_string('''
   fix BoxRelax all box/relax aniso 0.0 vmax 0.001
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

# if box lengths are within 0.01% move on

   if 0.9999 < xmax1/xmax2 < 1.0001 and 0.9999 < ymax1/ymax2 < 1.0001 and 0.9999 < zmax1/zmax2 < 1.0001:
      return

# if box lengths are not within 0.01%, rerun until within 0.01%   

   else:
      box_relax(lmp)