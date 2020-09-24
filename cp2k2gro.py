#!/usr/bin/python
import sys
import os
import numpy as np
if len(sys.argv)==1:
    if os.path.exists("cp2k-pos-1.xyz"):
        xyz="cp2k-pos-1.xyz"
        if os.path.exists("cp2k-1.cell"):
            cell="cp2k-1.cell"
        else:
            print("No cell found.")
            sys.exit(1)
    else:
        print("No position found.")
        sys.exit(1)
elif len(sys.argv)==2:
    xyz=sys.argv[1]
    if not os.path.exists(xyz):
        print("No position found.")
        sys.exit(1)
    else:
        if os.path.exists("cp2k-1.cell"):
            cell="cp2k-1.cell"
        else:
            print("No cell found.")
            sys.exit(1)               
else:
    xyz=sys.argv[1]
    cell=sys.argv[2]
    if not os.path.exists(xyz):
        print("No position found.")
    if not os.path.exists(cell):
        print("No cell found.")

masses= [1.00794,4.002602,6.941,9.012182,10.811,12.011,14.00674,\
     15.9994,18.9984032,20.1797,22.989768,24.3050,26.981539,28.0855,\
     30.97362,32.066,35.4527,39.948,39.0983,40.078,44.955910,47.88,\
      50.9415,51.9961,54.93085,55.847,58.93320,58.69,63.546,65.39,\
      69.723,72.61,74.92159,78.96,79.904,83.80,85.4678,87.62,88.90585,91.224,\
      92.90638,95.94,98,101.07,102.90550,106.42,107.8682,112.411,\
      114.82,118.710,121.75,127.60,126.90447,131.29,132.90543,137.327,138.9055,140.115,\
      140.90765,144.24,145,150.36,151.965,157.25,158.92534,162.50,\
      164.93032,167.26,168.93421,173.04,174.967,178.49,180.9479,\
      183.85,186.207,190.2,192.22,195.08,196.96654,200.59,204.3833,\
      207.2,208.98037,209,210.0,222,223,226.025,227.028,232.0381,\
      231.03588,238.0289,237.048,244.,243.,247.,247.,251.,252.,\
      257.,258.,259.,260.0]


element=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',\
 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',\
  'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', \
  'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', \
  'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', \
  'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', \
  'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', \
  'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', \
  'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

def write_gro(time,atoms,newlattice):
    with open("system.gro",'a+') as writer:
        writer.write("t= %11.4f\n%d\n" %(time/1000.0,len(atoms))) #//ps
        ind=1
        for i in atoms:
            writer.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(ind,"MOL",i[0],ind,i[1][0]/10.0,i[1][1]/10.0,i[1][2]/10.0))
            ind+=1
        writer.write("%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n" %(newlattice[0][0]/10.0,newlattice[1][1]/10.0,
        newlattice[2][2]/10.0,0.0,0.0,newlattice[1][0]/10.0,0.0,newlattice[2][0]/10.0,newlattice[2][1]/10.0))
    #v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y), the last 6 values may be omitted (they will be set to zero). GROMACS only supports boxes with v1(y)=v1(z)=v2(z)=0. 

cells=np.loadtxt(cell)
if len(cells.shape) == 1:
    cells= cells[np.newaxis,:]
if os.path.exists("system.gro"):
    os.remove("system.gro")
with open(xyz,'r') as reader:
    atom_numbers=int(reader.readline())
    reader.seek(0)
    while(reader.tell() < os.path.getsize(xyz)):
        reader.readline()
        info=reader.readline().split()
        time=float(info[5].strip(","))
        lattice_infos=cells[cells[:,1]==time].tolist()[0]
        newlattice=[]
        newlattice.append(list(map(float,lattice_infos[2:5])))
        newlattice.append(list(map(float,lattice_infos[5:8])))
        newlattice.append(list(map(float,lattice_infos[8:11])))
        atoms=[]
        for i in range(atom_numbers):
            info=reader.readline().split()
            atom_label=info[0]
            atom_xyz=list(map(float,info[1:]))
            atoms.append((atom_label,atom_xyz))
        #print(atoms)
        write_gro(time,atoms,newlattice)

def elementmass(element_label):
    index=0
    for i in element:
        if i.upper()==element_label.upper():
            return masses[index]
        index+=1
             
all_element__=list(set([i[0] for i in atoms]))
tobe_writen=""
tobe_writen+="; WARNING, it cannot be used for a simulation.\n\n"
tobe_writen+="[ defaults ]\n; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n1 3 yes 0.5 0.5\n\n"
tobe_writen+="[ atomtypes ]\n; name bond_type mass charge ptype sigma epsilon\n"
for i in all_element__:
    tobe_writen+="%s C 1.0 0.0 A 0.0 0.0\n" %(i)
tobe_writen+="\n"
tobe_writen+="; i j func b0 kb\nC C 1 0.13 1000.0 ; totally bogus\n\n"
tobe_writen+="[ angletypes ]\n; i j k func th0 cth\n  C C C 1 109.500 100.0 ; totally bogus\n\n"
tobe_writen+="[ dihedraltypes ]\n; i j k l func coefficients\n  C C C C 1 0.0 3 10.0 ; totally bogus\n\n"
tobe_writen+="[ moleculetype ]\n; Name      nrexcl\nmolecule0     3\n\n"
tobe_writen+="[ atoms ]\n; nr  type  resnr residue atom cgnr charge  mass\n"
ii = 1
for i in atoms:
    tobe_writen+="%6d%12s      1      MOL%7s%6d     0.0000%11.4f\n" %(ii,i[0],i[0],ii,elementmass(i[0]))
    ii+=1
tobe_writen+="\n"
tobe_writen+="[ bonds ]\n; i  j  func\n\n[ system ]\n; Name\nMymolecule0\n\n"
tobe_writen+="[ molecules ]\n; Compound    #mols\nmolecule0    1\n"
with open("system.top",'w') as writer:
    writer.write(tobe_writen)
with open("system.mdp",'w') as writer:
  
    writer.write("; minim.mdp - used as input into grompp to generate em.tpr\n")
    writer.write("integrator  = steep     ; Algorithm (steep = steepest descent minimization)\n")
    writer.write("emtol       = 1000.0    ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm\n")
    writer.write("emstep      = 0.01      ; Energy step size\n")
    writer.write("nsteps      = 50000     ; Maximum number of (minimization) steps to perform\n")
    writer.write("\n")
    writer.write("; Parameters describing how to find the neighbors of each atom and how to calculate the interactions\n")
    writer.write("nstlist         = 1         ; Frequency to update the neighbor list and long range forces\n")
    writer.write("cutoff-scheme   = Verlet\n")
    writer.write("ns_type         = grid      ; Method to determine neighbor list (simple, grid)\n")
    writer.write("coulombtype     = cut-off       ; Treatment of long range electrostatic interactions\n")
    writer.write("rcoulomb        = 0.01       ; Short-range electrostatic cut-off\n")
    writer.write("rvdw            = 0.01       ; Short-range Van der Waals cut-off\n")
    writer.write("pbc             = xyz       ; Periodic Boundary Conditions (yes/no)\n\n")