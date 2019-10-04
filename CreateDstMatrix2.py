import sys
import numpy as np
from math import sqrt
print (sys.argv[1])
if len(sys.argv)<4:
    print(len(sys.argv), "Too few arguments. There must be 4 parameters: matrix size (beads), path to lmp file, path to output matrix and path to output mol2 file.")
    sys.exit()
elif len(sys.argv)>4:
    print(len(sys.argv), "Too many arguments. There must be 4 parameters: matrix size (beads), path to lmp file, path to output matrix and path to output mol2 file.")
    sys.exit()
    
matr_size, path_to_lmp, path_to_matrix, path_to_mol2 = sys.argv

class bead:
    ''' class bead contains all info about bead in chain: global number, 
    the remaining valence, type of the bead and coordinates '''
    numbd=int()
    valence=int()
    typep=int()
    x=float()
    y=float()
    z=float()
    def __lt__(self, other):
        return self.numbd < other.numbd
    
class bond:
    '''class bond contains all info about bond: which beads connected by this bond'''
    first=int()
    last=int()
    
class chain:
    '''class chain has two lists of beads and bonds and general info about system such as total number of particles, density and box size along each axis'''
    def __init__(self):
        self.bd=[]
        self.bnd=[]
        self.number_of_beads=float()
        self.number_of_bonds=float()
        self.density=float()
        self.xbox=float()
        self.ybox=float()
        self.zbox=float()

def readLmpRst (f, polymer):
    one=bead()
    sb=bond()
    flag_bead=False
    flag_bond=False
    for i,line in enumerate(f):
        if i>1:
            if not line.strip():
                continue
            if 'Atoms # angle' in line:
                flag_bond=False
                flag_bead=True
                continue
            if 'Bonds' in line:
                flag_bond=True
                flag_bead=False
                continue
            if 'Velocities' in line:
                flag_bond=False
                flag_bead=False
                continue
            if 'atoms' in line:
                head,tail = line.split()
                polymer.number_of_beads=head
            elif 'bonds' in line:
                head,tail = line.split()
                polymer.number_of_bonds=head
            elif 'xlo xhi' in line:
                q,w,e,r=line.split()
                polymer.xbox=float(w)
            elif 'ylo yhi' in line:
                q,w,e,r=line.split()
                polymer.ybox=float(w)
            elif 'zlo zhi' in line:
                q,w,e,r=line.split()
                polymer.zbox=float(w)
            elif flag_bead==True:
                one=bead()
                numbd,valence,typep,x,y,z,t1,t2,t3 = line.split()
                one.numbd=int(numbd)
                one.valence=int(valence)
                one.typep=int(typep)
                one.x=float(x)
                one.y=float(y)
                one.z=float(z)
                polymer.bd.append(one)
            elif flag_bond==True:
                sb=bond()
                num, typeb, head, tail=line.split()
                sb.first=int(head)
                sb.last=int(tail)
                polymer.bnd.append(sb)        

def removePBC (polymer):
    '''revome periodic boundary conditions in case of single chain and numbers of beads corresspond to the global numbers and start from 1'''
    itx=0
    ity=0
    itz=0
    for i in range(len(polymer.bd)-1):
        if polymer.bd[i].x - itx * polymer.xbox - polymer.bd[i+1].x > polymer.xbox / 2:
            itx = itx + 1
            polymer.bd[i+1].x = polymer.bd[i+1].x + itx*polymer.xbox
        elif polymer.bd[i].x - itx * polymer.xbox - polymer.bd[i+1].x < -polymer.xbox / 2:
            itx = itx - 1
            polymer.bd[i+1].x = polymer.bd[i+1].x + itx*polymer.xbox
        else:
            polymer.bd[i+1].x = polymer.bd[i+1].x + itx*polymer.xbox
            
        if polymer.bd[i].y - ity * polymer.ybox - polymer.bd[i+1].y > polymer.ybox / 2:
            ity = ity + 1
            polymer.bd[i+1].y = polymer.bd[i+1].y + ity*polymer.ybox
        elif polymer.bd[i].y - ity * polymer.ybox - polymer.bd[i+1].y < -polymer.ybox / 2:
            ity = ity - 1
            polymer.bd[i+1].y = polymer.bd[i+1].y + ity*polymer.ybox
        else:
            polymer.bd[i+1].y = polymer.bd[i+1].y + ity*polymer.ybox
            
        if polymer.bd[i].z - itz * polymer.zbox - polymer.bd[i+1].z > polymer.zbox / 2:
            itz = itz + 1
            polymer.bd[i+1].z = polymer.bd[i+1].z + itz*polymer.zbox
        elif polymer.bd[i].z - itz * polymer.zbox - polymer.bd[i+1].z < -polymer.zbox / 2:
            itz = itz - 1
            polymer.bd[i+1].z = polymer.bd[i+1].z + itz*polymer.zbox
        else:
            polymer.bd[i+1].z = polymer.bd[i+1].z + itz * polymer.zbox

def distance (a, b):
    '''calculate distance between two beads'''
    return sqrt((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2)

def FillDistanceMatrix(poly,i,j,matr):
    
    coarse=int(i/len(matr))
    gap=len(matr)*coarse
    for iter1 in range(i,j):
        for iter2 in range(iter1+1,j):
            matr[iter1-gap][iter2-gap]+=distance(poly.bd[iter1],poly.bd[iter2])

def writeMol2 (polymer, path):
    bstr='1  ala'
    the_file=open(path, 'w')
    the_file.write('@<TRIPOS>MOLECULE\n')
    the_file.write('mol_name\n')
    the_file.write('\t %d \t %d \t %s \t %s \t %s \n' %(len(polymer.bd), len(polymer.bnd), '0', '0', '0'))
    the_file.write('SMALL\n')
    the_file.write('USER_CHARGES\n')
    the_file.write('@<TRIPOS>ATOM\n')
    for i in range(len(polymer.bd)):
        ty='O'
        if polymer.bd[i].typep==1:
            ty='C'
        the_file.write('%d \t %s \t %f \t %f \t %f \t %s \t %s \t %f \n' %(i+1, ty, polymer.bd[i].x, polymer.bd[i].y, polymer.bd[i].z, ty, bstr, float(i)))
    the_file.write('@<TRIPOS>BOND\n')
    for i in range(len(polymer.bnd)):
        the_file.write('%d \t %d \t %d \t %s \n' %(i+1, polymer.bnd[i].first, polymer.bnd[i].last, '1'))
    the_file.close()
    
matr=np.zeros((matr_size,matr_size))
poly=chain()
f=open(path_to_lmp)
readLmpRst(f, poly)
poly.bd.sort()
removePBC(poly)
numMaps=int(len(poly.bd)/matr_size)
print("Lenght of the polymer chain is ",len(poly.bd))
for i in range(0,len(poly.bd),matr_size):
    if i+matr_size<len(poly.bd):
        FillDistanceMatrix(poly,i,i+matr_size,matr)
for i in range(matr_size):
    for j in range(i+1,matr_size):
        matr[i][j]=matr[i][j]/numMaps
        matr[j][i]=matr[i][j]
np.savetxt(path_to_matrix,matr,fmt="%10.5f")
writeMol2(poly,path_to_mol2)