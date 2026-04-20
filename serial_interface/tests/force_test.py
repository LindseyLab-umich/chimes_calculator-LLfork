#! /usr/bin/python3
#
#  Compares forces in the specified tests to finite differences of energy by running the C++ serial interface.
#  Run this program on a dedicated node only, it will use all CPU's.
#
import os
from os.path import exists
import glob
import sys
from multiprocessing import Pool

# Step size in finite difference calculations.
dx = 0.00025
#tests=[['published_params.liqCO.2+3b.cubic.txt','CO.2.5gcc_6500K.OUTCAR_force_#000.xyz','CO.2.5gcc_6500K.OUTCAR_force_#001.xyz','CO.2.5gcc_6500K.OUTCAR_force_#002.xyz','1'],
#       ["validated_params.TiO2.2+3b.Tersoff.txt", "TiO2.unitcell_arbrot_force_#000.xyz", "TiO2.unitcell_arbrot_force_#001.xyz", "TiO2.unitcell_arbrot_force_#002.xyz","0"]]
#tests=[["validated_params.TiO2.2+3b.Tersoff.txt", "TiO2.unitcell_arbrot_force_#000.xyz","0"],
#       ['published_params.liqCO.2+3b.cubic.txt','CO.2.5gcc_6500K.OUTCAR_force_#000.xyz',"1"],
#       ["published_params.liqCO.2+3b.cubic.txt","diam.2_#000.xyz","1"]
#       ]

# Tests to run.  Format is the same as run_tests.sh
tests=[
["published_params.liqC.2b.cubic.txt","liqC.2.5gcc_6000K.OUTCAR_#000.xyz","0"],
["published_params.liqC.2+3b.cubic.txt","liqC.2.5gcc_6000K.OUTCAR_#000.xyz","0"],
["published_params.liqCO.2+3b.cubic.txt","CO.2.5gcc_6500K.OUTCAR_#000.xyz","1"],
["validated_params.CO2400K.2+3+4b.Tersoff.special.offsets.relabel.txt","CO.2.5gcc_6500K.OUTCAR_#000.relabel.xyz","0"],
["published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt","CO.2.5gcc_6500K.OUTCAR_#000.scramble.xyz","0"],
["published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt","CO.2.5gcc_6500K.OUTCAR_#000.translate.xyz","0"],
["published_params.HN3.2+3+4b.Tersoff.special.offsets.txt","HN3.2gcc_3000K.OUTCAR_#000.xyz","0"],
["validated_params.TiO2.2+3b.Tersoff.txt","TiO2.unitcell_arbrot_#000.xyz","0"],
["test_params.CHON.txt","CHON.testfile.000.xyz","0"],
["published_params.liqCO.2+3b.cubic.txt","diam.64_#000.xyz","1"],
["published_params.liqCO.2+3b.cubic.txt","diam.16_#000.xyz","1"],
["published_params.liqCO.2+3b.cubic.txt","diam.8_#000.xyz","1"],
["published_params.liqCO.2+3b.cubic.txt","diam.2_#000.xyz","1"],
["published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt","CO.9GPa_2400K.OUTCAR_#000.xyz","0"],
["test_params.h2o_2bcheby.txt","H2O.1.50gcc_2000K_#000.xyz","0"],
["test_params.HN3.penalty.txt","HN3.2.04gcc_20000K_#000.xyz","0"]
]


exe='../../build/CPP-interface'


def find_energy(lines1, atom, xyz):
    found_e = 0
    for j in range(0, len(lines1)):
        line = lines1[j]
        line.rstrip()
        if found_e == 1:
            print("Energy is {:21.14e}".format(float(line)))
            energy = float(line)
            break 

        if "Energy (kcal/mol)" in line:
            found_e = 1

    if found_e == 0:
        print("Did not find energy in output")
        sys.exit(1)
        
    return energy

def find_force_all(lines1):
# Find all forces in the initial calculation.    
    found_force = 0
    forces=[]
    for j in range(0, len(lines1)):
        line = lines1[j]
        line.rstrip()

        if found_force == 1:
            found_force = 0
            for k in range(j, len(lines1) ):
                fields = lines1[k].split()
                if ( len(fields) == 3 ):
                    for l in range(0,3):
                        forces.append(float(fields[l]))
        
        if "Forces (kcal/mol/A)" in line:
            found_force = 1

    return forces

def run_chimes(args):
# Chimes calculation with atom position offset
#
           # Unpack args.
           cfgj, model, option, atom, xyz, stage = tuple(args)
           # print("Working on ", cfgj)
           filename = 'force_test.' + str(atom) + "." + str(xyz) + "." + str(stage) + '.log'
           command=exe + " " + "force_fields/" + model + " " + cfgj + " " + option + " > {fn}".format(fn=filename)

           if os.system(command) != 0:
               print(command, " failed")
               sys.exit(1)

def read_chimes_energy(atom, xyz, stage):
# Chimes calculation with atom position offset: return energy.
           filename = 'force_test.' + str(atom) + "." + str(xyz) + "." + str(stage) + '.log'
           try:
               file1 = open(filename, 'r')
           except IOError:
               print("Couldn't open ", filename)
               sys.exit(1)
               
           lines1 = file1.readlines()
           energy_tmp = find_energy(lines1, atom, xyz)
           file1.close()
           return energy_tmp

def init_chimes(cfgj, model, option):
# Initial chimes calculation with no position offset: return all forces.
           # print("Working on ", cfgj)
           filename = 'force_test0.log'
           command=exe + " " + "force_fields/" + model + " " + cfgj + " " + option + " > force_test0.log"

           print(command)
           if os.system(command) != 0:
               print(command, " failed")
               sys.exit(1)
               
           file1 = open(filename, 'r')
           lines1 = file1.readlines()
           force_tmp = find_force_all(lines1)
           file1.close()

           return force_tmp
       
def make_config_file(name, atom, xyz, delta, name2, lines):
## Create a single atom configuration file
    file2 = open(name2, "w")

    for j in range(0,atom+2):
        file2.write(lines[j])

    fields = lines[atom+2].split()

    pos = float(fields[xyz+1])
    if delta == 1:
        pos = pos + dx
        fields[xyz+1] = pos
        file2.write(str(fields[0]) + " " + str(fields[1]) + " " + str(fields[2]) + " " + str(fields[3]) + "\n")
    elif delta == 2:
        pos = pos - dx
        fields[xyz+1] = pos
        file2.write(str(fields[0]) + " " + str(fields[1]) + " " + str(fields[2]) + " " + str(fields[3]) + "\n")
    else:
        file2.write(lines[atom+2])

    for j in range(atom+3,len(lines)):
        file2.write(lines[j])

    file2.close
        
                
def make_input_files(cfg, dx, name):
# Create input configuration files with displaced coordinates for atoms.
# There is one file for each displacement.

    base_file="configurations/" + cfg
    if not exists(base_file):
        print(base_file, " does not exist.")
    file1 = open(base_file, "r")
    lines = file1.readlines()
    if ( len(lines) < 3 ):
        print("Less than 3 lines in ", base_file)
        exit(1)
    file1.close()

    # Clean up old files.
    os.system("rm {n}.[0-9]*.[0-9]*.[0-9]*.xyz".format(n=name))

    file_list = []
    # Write non-displaced input file
    name2 = name + "." + ".xyz"
    file_list.append(name2)
    make_config_file(name, 0, 0, 0.0, name2, lines)                
    
    # Write displaced input files.
    natoms = len(lines) - 2
    for atom in range(0,natoms):
        for xyz in range(0,3):
            for delta in range(1,3):
                name2 = name + "." + str(atom) + "." + str(xyz) + "." + str(delta) + ".xyz"
                make_config_file(name, atom, xyz, delta, name2, lines)
                file_list.append(name2)
                
    return file_list


#### MAIN PROGRAM #####
if os.path.isfile("force_test0.log"):
    os.remove("force_test0.log")    

    

for test in tests:
       model=test[0]
       cfg=test[1]
       option=test[2]

       print("Model = ", model)
       print("CFG   = ", cfg)
       print("option = ", option)

       
       os.system("rm -f force_test.*.log")
       os.system("rm -f force_test.*.xyz")
       
       file_list = make_input_files(cfg, dx, "force_test")
       
       energy=[]
       force=[]
       count = 0

       # Run chimes with no atom displacements.
       force_all = init_chimes(file_list[0], model, option)               

       chimes_args=[]

       # Package chimes arguments.
       for step in range(1,len(file_list)):
           atom = (step-1)//6
           xyz = (step-1)//2 % 3
           stage= (step-1) % 2
           arg = [file_list[step], model, option, atom, xyz, stage]
           chimes_args.append( arg )

        # Parallel execution of chimes calculations.
       pool = Pool(os.cpu_count())
       pool.map(run_chimes, chimes_args)

       # Read the output files.
       for step in range(1,len(file_list)):
           atom = (step-1)//6
           xyz = (step-1)//2 % 3
           stage= (step-1) % 2           
           energy_tmp = read_chimes_energy(atom, xyz, stage)
           energy.append(energy_tmp)           
           if ( (step - 1) % 2 == 1 and step > 1 ) :
               test_force = ( energy[step-1] - energy[step-2] ) / (2.0 * dx)

               print("Atom = ", str(atom))
               print("XYZ component = ", str(xyz))
               print("Finite difference force = ", test_force)
               print("Analytic force          = ", force_all[count])
               if abs((force_all[count]-test_force)) / (1.0 + abs(test_force)) < 1.0e-04:
                   print("Test passed\n")
               else:
                   print("Test failed\n")
                   sys.exit(1)
               count = count + 1


os.system("rm -f force_test.*.log")
os.system("rm -f force_test.*.xyz")
