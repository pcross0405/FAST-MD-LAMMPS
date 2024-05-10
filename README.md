-------------------------------------------------------------------------------------------------------  
FAST-MD (22 Nov. 2023)
-------------------------------------------------------------------------------------------------------  

A package that enables molecular dynamic simulations to incorporate the FAST principle via LAMMPS.  

-------------------------------------------------------------------------------------------------------  
--------------------------------------- INSTALLATION INSTRUCTIONS -------------------------------------  
-------------------------------------------------------------------------------------------------------  

1) Make a directory to build the package in

2) Inside that directory type on the command line  
   "git clone https://github.com/pcross0405/FAST-MD-LAMMPS.git"

3) Make sure you have python's build tool up to date with  
   "python3 -m pip install --upgrade build"

4) Once up to date type  
   "python3 -m build"

5) This should create a "dist" directory with a .whl file inside

6) On the command line type  
   "pip install dist/'*.whl'" 

-------------------------------------------------------------------------------------------------------  
---------------------------------------- INSTRUCTIONS FOR USE -----------------------------------------  
-------------------------------------------------------------------------------------------------------  

To use any of the modules in this package in a python script type  
"name-of-module.name-of-module(args)"

Of course, replace name-of-module with the name of the module you wish to use  

The arguments "(args)" that a module requires can be found by typing  
"vi /path/to/python-directory/lib/python-version-number/site-packages/fast_md/name-of-module"  

Replace /path/to/python-directory/ with the path to your python interpreter  

Replace the "version-number" in /python-version-number/ with the version number of your python installation  
	(python version number can be found by typing "python3 --version", generally this will print something  
	 like "Python x.y.z"; you need to change /python-version-number/ -> /python-x.y/)  

Replace /name-of-module with the name of the module you wish to use (i.e. /name-of-module -> /check_center.py)  

---------------------------------------------------------------------------------------------------------  
------------------------------------------- REPORTING ISSUES --------------------------------------------  
---------------------------------------------------------------------------------------------------------  

Please report any issues to "https://github.com/pcross0405/FAST-MD-LAMMPS/issues"  

-------------------------------------------------------------------------------------------------------------------------  
-- SEE SAMPLES DIRECTORY FOR AN EXAMPLE OF HOW TO IMPLEMENT THE PACKAGE IN A PYTHON SCRIPT THAT INTERFACES WITH LAMMPS --  
-------------------------------------------------------------------------------------------------------------------------  
