import sys
import copy
import os
import subprocess as sp

# Reads the Parameterfile from CGENFF
def read_ParameterFile(fname,out):
    ring=[]
    atoms={}
    bonds={}
    penalties={'ALL':[0,0,0],'ANGLES':[0,0,0],'CHARGES':[0,0,0],'DIHEDRALS':[0,0,0]}
    angle,di=False,False
    angles,dihedrals={},{}
    with open(fname,'r') as f:
        for line in f:
            l=line.split()
            if line[0:4]=='ATOM':
                atoms.update({l[1]:l[2:4]}))
                if 'R' in l[2]:
                    ring.append(l[1])
                p=l[5]
                if '!' in p:
                    p.remove('!')
                penalties['CHARGES'][whatPenalty(float(p))]+=1
                penalties['ALL'][whatPenalty(float(p))]+=1
            elif line[0:5]=='BOND ':
                if l[1] not in bonds:
                    bonds.update({l[1]:[l[2]]})
                elif l[2] not in bonds[l[1]]:
                    bonds[l[1]]+=[l[2]]
                if l[2] not in bonds:
                    bonds.update({l[2]:[l[1]]})
                elif l[1] not in bonds[l[2]]:
                    bonds[l[2]]+=[l[1]]
            elif line[0:6] =='ANGLES':
                angle=True
            elif line[0:9]=='DIHEDRALS':
                angle=False
                di=True
            elif len(l)<1:
                di=False
            elif angle:
                angles.update({tuple(l[:3]):"{0:.1f}".format(float(l[4]))})
                if l[-1].replace('.','').isdigit():
                    penalties['ANGLES'][whatPenalty(float(l[-1]))]+=1
                    penalties['ALL'][whatPenalty(float(l[-1]))]+=1
            elif di:
                if tuple(l[:4]) not in dihedrals:
                    dihedrals.update({tuple(l[:4]):"{0:.1f}".format(float(l[6]))})
                if l[-1].replace('.','').isdigit():
                    penalties['DIHEDRALS'][whatPenalty(float(l[-1]))]+=1
                    penalties['ALL'][whatPenalty(float(l[-1]))]+=1
    # Write Penalty File
    with open(out,'w') as o:
        names=['CHARGES','ANGLES','DIHEDRALS','ALL']
        for name in names:
            o.write(name+'\n')
            o.write(str(penalties[name])+'\n')

    return(atoms, bonds, ring,angles,dihedrals)

# Determines in which category the penalty belongs
def whatPenalty(n):
    if n<10.0:
        return 0
    elif n<50.0:
        return 1
    else:
        return 2

# reads Molfile to get the coordinates of each atom
def read_Molfile(fname,atoms):
    count_H=0
    doubles={}
    with open(fname,'r') as f:
        for line in f:
            l=line.split()
            if len(l)>1:
                if l[1] in atoms:
                    atoms[l[1]].append(l[2:5])
                if l[1]=='H':
                    count_H+=1
                    if 'H'+str(count_H) in atoms:
                        atoms['H'+str(count_H)].append(l[2:5])
                        atoms['H'+str(count_H)].append([])

# Brings the atom names in the correct order
def sort_names(atoms,bonds):
    names=[]
    keys=list(atoms.keys())
    keys.sort(key=lambda x: x[1:])
    for atom in keys:
        if (atom[0]!='H'):
            names.append(atom)
            for i in range(len(bonds[atom])):
                if bonds[atom][i][0]=='H':
                    names.append(bonds[atom][i])
    return names


# Returns bounded atoms in Dowser format
def boundedAtoms(atom,bonds):
    if 'H'==atom[0]:
        bonds[atom]=[bonds[atom][0],'NOT']
        return bonds
    else:
        bounded=bonds[atom]
        hlist=[x for x in bounded if x[0]=='H']
        bounded=[item for item in bounded if item not in hlist]
        bounded.sort(key=lambda x: x[1:])
        before,after='NOT','NOT'
        remove=[]
        for elem in bounded:
            if  (elem=='SAW') | (elem=='NOT'):
                continue
            elif (len(bonds[elem]) <3):
                continue
            elif (bonds[elem][0]==atom) & (bonds[elem][2]=='SAW'):
                after=elem
                remove.append(elem)
            elif (bonds[elem][1]==atom) & (bonds[elem][2]=='SAW'):
                before=elem
                remove.append(elem)
        bounded=[item for item in bounded if item not in remove]
        if (len(bounded)>0) & (before=='NOT'):
            before=bounded[0]
            remove.append(before)
        bounded=[item for item in bounded if item not in remove]
        if (len(bounded)>0) & (after=='NOT'):
            after=bounded[-1]
        elif (len(hlist)>0) & (after=='NOT'):
            after=hlist[0]
        bonds[atom]=[before,after,'SAW']
    return bonds

# Returns the corresponding residue name in atomparams.db
def residue(atom,atoms,bonds,ring):
    if atom[0]=='H':
        return 'H'
    # aromatic C
    elif (atom[0]=='C') & (atom in ring) & (len([x for x in bonds[atom] if x[0]=='H'])==1):
        return 'CHR'
    elif (atom[0]=='C')  & (len([x for x in bonds[atom] if x[0]=='H'])==1):
        return 'CH1'
    elif (atom[0]=='C')  & (len([x for x in bonds[atom] if x[0]=='H'])==2):
        return 'CH2'
    elif (atom[0]=='C') & (len([x for x in bonds[atom] if x[0]=='H'])==3):
        return 'CH3'
    elif (atom[0]=='C') & (len(bonds[atom])==3):
        return 'CR'
    elif atom[0]=='C':
        return 'CA'
    elif (atom[:2]=='NA') & (atoms[atom][0][:2]=='NA'):
        return 'NA'
    elif atom[0]=='N':
        return 'N'
    elif atom[0]=='S':
        return 'S'
    elif atom[0]=='O':
        return 'O'
    elif atom[:2]=='ZN':
        return 'ZN'
    elif atom[:3]=='CAL':
        return 'CAL'

# Makes sure everything has a length of given n
def addspace(word,n,begin=True):
    if begin:
        while len(word)<n:
            word+=' '
    else:
        while len(word)<n:
            word=' '+word
    return word

# Calculates the euclidean distance of two molecules
def distance(atom1,atom2):
    return "{0:.3f}".format(((sum([(x*x) for x in [float(atom1[2][i])-float(atom2[2][i]) for i in range(3) ]]))**0.5))

# Return Tuple of atoms of a backbone starting at given atom of the length n
def goBack(atom,atoms,bonds,names,n):
    state=()
    while len(state)<n:
        if (atom=='NOT') | (atom=='SAW'):
            return 0
        state=(atoms[atom][0],)+state
        atom=bonds[atom][0]
    return state

# Creates the Paramterfile for Dowser++
def make_DowserParamFile(output,atoms, bonds, ring,ligand,angles,dihedrals):
    names=sort_names(atoms,bonds)
    bond2=copy.copy(bonds)
    for atom in names:
        bonds=boundedAtoms(atom,bonds)
    with open(output, 'a') as f:
        f.write('RESIDUE '+ligand+'\n')
        for atom in names:
            partialCharge=atoms[atom][1]
            if float(partialCharge)>=0:
                partialCharge=' '+partialCharge
            dist='0.000'
            angle,dihedral='0.0','0.0'
            a=goBack(atom,atoms,bonds,names,3)
            b=goBack(atom,atoms,bonds,names,4)
            if (a!=0) & (a in angles):
                angle=angles[a]
            if (b!=0) & (b in dihedrals):
                dihedral=dihedrals[b]
            if bonds[atom][0]!='NOT':
                dist=distance(atoms[atom],atoms[bonds[atom][0]])
            f.write('ATOM '+ligand+'  '+addspace(atom,3)+'  '+addspace(bonds[atom][0],5) + addspace(bonds[atom][1],6)+dist+'  '+addspace(angle,5,False)+'  '+addspace(dihedral,5,False)+'  '+partialCharge+' '+residue(atom,atoms,bond2,ring)+'\n')

# Execute command on command line
def cmd(args):
    p=sp.Popen(args)
    p.wait()

# Returns a list of ligand names and a list with all ligands, including different chains and residues
def get_ligands(pdb,PATH):
    if not os.path.isfile(pdb):
        pdb=PATH+'/'+pdb
    ligands=[]
    all_ligands=[]
    with open(pdb,'r') as f:
        for line in f:
            if line[:6]=='HETATM':
                if line[17:20] not in ligands:
                    ligands.append(line[17:20])
                if (line[17:20],line[21],line[23:26]) not in all_ligands:
                    all_ligands.append((line[17:20],line[21],line[23:26]))
    with open(PATH+'Ligands.txt','w') as o:
        for ligand in ligands:
            o.write(ligand+'\n')
    return ligands, all_ligands

# Creates a seperate pdb file for a ligand, so it can be used for parameterization
def get_ligandpdb(pdb,lfile,ligand, PATH, resi='', chain='',start=True):
    if not os.path.isfile(pdb):
        pdb=PATH+'/'+pdb
    with open(pdb,'r') as f:
        with open(lfile,'w') as o:
            for line in f:
                if line[:6]=='HETATM':
                    if (line[17:20]==ligand) & (((line[23:26]==resi) & (chain==line[21]))| (start)):
                        if start:
                            resi=line[23:26]
                            chain=line[21]
                            start=False
                        o.write('ATOM  '+line[6:])

# Creates the Dowser and Dowser++ pdb, where water of the ligands are included
# Dowser++ has not HETATMs, only ATOMs
def makePdb(pdb,ligands,all_ligands,PATH):
    pdb2=copy.copy(pdb)
    if not os.path.isfile(pdb):
        pdb2=PATH+'/'+pdb
    with open(pdb2, 'r') as f:
        with open(PATH+'Water_'+pdb,'w') as o:
            for line in f:
                if (line[:6] not in ['HETATM','ANISOU']):
                    o.write(line)
                elif (line[:6] in ['HETATM','ANISOU']) & ((line[17:20],line[21],line[23:26]) in all_ligands):
                    ligand=line[17:20]+line[23:26]+line[21]
                    lfile=PATH+ligand+'.pdb'
                    wlfile=PATH+'Water_'+ligand+'.pdb'
                    molfile=PATH+ligand+'.mol2'
                    get_ligandpdb(pdb,lfile,line[17:20],PATH,line[23:26],line[21],False)
                    cmd(['babel','-i','pdb',lfile,'-o','mol2',molfile,'-h'])
                    cmd(['babel','-i','mol2',molfile,'-o','pdb',wlfile,'-h'])
                    all_ligands.remove((line[17:20],line[21],line[23:26]))
                    with open(wlfile,'r') as ff:
                        count=1
                        for line2 in ff:
                            if line2[:4]=='ATOM':
                                new=line2.replace(' A ',' '+line[21] +' ')
                                if line2[13]=='H':
                                    new=new.replace('H  ','H'+addspace(str(count),2),1)
                                    count+=1
                                o.write(new)
                    os.remove(lfile)
                    os.remove(wlfile)
                    os.remove(molfile)

    with open(PATH+'Water_'+pdb,'r') as f:
        with open(PATH+'Dowser_'+pdb,'w') as o:
            with open(PATH+'Dowser++_'+pdb,'w') as oo:
                count=0
                for line in f:
                    new=line
                    new2=line
                    if line[:4] in ['ATOM','TER ']:
                        count+=1
                    if line[17:20] in ligands:
                        new=line[:6]+addspace(str(count),5,False)+line[11:]
                        new2=new.replace('ATOM  ','HETATM')
                    o.write(new2)
                    oo.write(new)
    os.remove(PATH+'Water_'+pdb)


def main(pdb,outputPath):
    outputPath+='/'
    ligands,all_ligands=get_ligands(pdb,outputPath)
    for ligand in ligands:
        cgenff=outputPath+ligand+'_cgenff.txt'
        lfile=outputPath+ligand+'.pdb'
        wlfile=outputPath+'Water_'+ligand+'.pdb'
        molfile=outputPath+ligand+'.mol2'
        output=outputPath+ligand+'.db'
        get_ligandpdb(pdb,lfile,ligand,outputPath)
        cmd(['babel','-i','pdb',lfile,'-o','mol2',molfile,'-h'])
        cmd(['/home/student/Software/cgenff/cgenff',molfile, '-f',cgenff])
        atoms, bonds, ring,angles,dihedrals=read_ParameterFile(cgenff,outputPath+ligand+'_Penalty.txt')
        read_Molfile(molfile,atoms)
        make_DowserParamFile(output,atoms,bonds,ring,ligand,angles,dihedrals)
        cmd(['babel','-i','mol2',molfile,'-o','pdb',wlfile,'-h'])
    makePdb(pdb,ligands,all_ligands,outputPath)



if __name__ == "__main__":
    if len(sys.argv)!=3:
        print("Error! Too many or too less arguments! Two files are needed.")
    else:
        main(sys.argv[1],sys.argv[2])