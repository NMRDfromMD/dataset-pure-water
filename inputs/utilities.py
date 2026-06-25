import numpy as np
import random
import copy
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R

def neighborsearch(neighbor, molecule, x, y, z, Lx, Ly, Lz):
    """Search all neighbor to a molecule in a box and return the closest distance."""
    box = np.array([Lx, Ly, Lz])
    minr = 10
    for m in molecule:
        x0 = m[0] + x
        y0 = m[1] + y
        z0 = m[2] + z
        dxdydz = np.remainder(neighbor - np.array([x0,y0,z0]) + box/2., box) - box/2.
        minr = np.min([minr,np.min(norm(dxdydz,axis=1))])
    return minr

def randomlocation(Lx,Ly,Lz):
    """Choose a random location within a given box."""
    txlo, txhi = -Lx/2, Lx/2
    tylo, tyhi = -Ly/2, Ly/2
    tzlo, tzhi = -Lz/2, Lz/2    
    x = random.randint(1,1000)/1000*(txhi-txlo)
    y = random.randint(1,1000)/1000*(tyhi-tylo)
    z = random.randint(1,1000)/1000*(tzhi-tzlo)
    return x, y, z

def place_molecules(NH2O):

    attempt = 0
    cptH2O = 0

    while (cptH2O<NH2O):

        ## choose initial box dimensions
        Lx, Ly, Lz = 3.5+0.5*attempt, 3.5+0.5*attempt, 3.5+0.5*attempt # nm
        
        box = np.array([Lx, Ly, Lz])

        ## initialise matrix
        atoms = np.zeros((100000,7))
        cptatoms = 0
        cptres = 0
        resName = ["" for x in range(1000000)]
        atoName = ["" for x in range(1000000)]

        ### place water molecules

        rOH = 0.9572/10 
        rOM = 0.105/10 # tip4p/epsilon
        thetaHOH = 104.52    
        atomH2O = np.array([[1, 1, 0, 0, 0, 0], \
        [2, 2, 0.527, rOH*np.cos((thetaHOH/2)*np.pi/180),   rOH*np.sin((thetaHOH/2)*np.pi/180), 0.0], \
        [3, 2, 0.527, rOH*np.cos((thetaHOH/2)*np.pi/180),   -rOH*np.sin((thetaHOH/2)*np.pi/180),  0.0], \
        [4, 3, -1.054, rOM,  0.0, 0.0]])
        cptH2O = 0
        cptH = 0
        fail_attempt = 0
        while (cptH2O < NH2O) & (fail_attempt <1e6):
            x,y,z = randomlocation(Lx,Ly,Lz)
            pos = np.array([x,y,z])

            if cptH2O > 0:
                dxdydz = np.remainder((pos - atoms[:cptatoms].T[4:].T) + box/2., box) - box/2.
                d = np.min(norm(dxdydz,axis=1))
            else:
                d = 1

            if d > 0.24:
                for m in atomH2O:
                    atoms[cptatoms] = cptatoms+1, cptres+1, m[1], m[2], m[3]+x, m[4]+y, m[5]+z
                    resName[cptatoms] = 'SOL'
                    if m[1] == 1:
                        atoName[cptatoms] = 'OW'
                    elif m[1] == 2:
                        if cptH % 2 == 0:
                            atoName[cptatoms] = 'HW1'
                        else:
                            atoName[cptatoms] = 'HW2'
                        cptH += 1
                    elif m[1] == 3:
                        atoName[cptatoms] = 'MW'
                    cptatoms += 1    
                cptH2O += 1
                cptres += 1
            else:
                fail_attempt += 1
        attempt += 1
    
    print('box size = '+str(Lx)+' nm')

    return atoms, atoName, resName, Lx, Ly, Lz

def write_lammps(atoms, bonds, angles, Lx, Ly, Lz):

    f = open("initial.data", "w")
    f.write('# LAMMPS data file \n\n')
    f.write(str(len(atoms))+' atoms\n')
    f.write(str(len(bonds))+' bonds\n')
    f.write(str(len(angles))+' angles\n')
    f.write('\n')
    f.write('2 atom types\n')
    f.write('1 bond types\n')
    f.write('1 angle types\n')
    f.write('\n')
    f.write('0 '+str(Lx*10)+' xlo xhi\n')
    f.write('0 '+str(Ly*10)+' ylo yhi\n')
    f.write('0 '+str(Lz*10)+' zlo zhi\n')
    f.write('\n')
    f.write('Atoms\n')
    f.write('\n')
    for myatom in atoms:
        if myatom[2] < 3:
            if myatom[2] == 1:
                myatom[3] = -2*0.527
            myatom[4] *= 10
            myatom[5] *= 10
            myatom[6] *= 10
            for col in range(len(myatom)):
                if col < 3:
                    f.write(str(int(myatom[col]))+' ')
                else :
                    f.write(str(myatom[col])+' ')
            f.write('\n')
    f.write('\n')
    f.write('Bonds\n')
    f.write('\n')
    for cpt, mybond in enumerate(bonds):
        id1 = mybond[0]
        id2 = mybond[1]
        bond_types = 1
        myline = [cpt + 1, bond_types, id1, id2]
        for col in range(len(myline)):
            f.write(str(int(myline[col]))+' ')
        f.write('\n')
    f.write('\n')
    f.write('Angles\n')
    f.write('\n')   
    for cpt, myangle in enumerate(angles):
        id1 = np.int32(myangle[0])
        id2 = np.int32(myangle[1])
        id3 = np.int32(myangle[2])
        angle_types = 1
        if angle_types is not None:
            myline = [cpt + 1, angle_types, id1, id2, id3]
            for col in range(len(myline)):
                f.write(str(int(myline[col]))+' ')
            f.write('\n')
    f.close()

def prepare_lammps(atoms):
    # GROMACS to LAMMPS

    # suppress the extra oxygen
    atoms = atoms[atoms.T[2] < 3]

    # renumber the matrix
    cpt1 = 1
    new_atoms = []
    for atom in atoms:
        if np.sum(atom) == 0:
            continue
        newid = cpt1
        atom[0] = newid
        # replace 
        new_atoms.append(atom)
        cpt1 += 1
    atoms = np.array(new_atoms)

    # add bond for water
    new_bonds = []
    new_angles = []
    cpt = 1
    for atom in atoms:
        if atom[2] == 1:
            new_bonds.append([cpt, cpt+1])
            new_bonds.append([cpt, cpt+2])
            new_angles.append([cpt+1, cpt, cpt+2])
        cpt += 1
    bonds = np.array(new_bonds)
    angles = np.array(new_angles)

    return atoms, bonds, angles
