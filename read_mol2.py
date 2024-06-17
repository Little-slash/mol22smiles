def read_mol2(filename):
    atom_list = []
    bond_list = []
    dika_list =[]
    with open(filename, 'r') as f:
        flag = [1]
        while True:
            line = f.readline()
            if not line:
                break
            if line.find("@<TRIPOS>BOND") != -1:
                flag[0] = 2
                continue
            if "ATOM" in line:
                flag[0] = 3
                continue
            line.strip()
            line.replace('\n', '')
            sp = line.split()
            if flag[0] == 2 and len(sp) != 0:
                if sp[3]=="Ar" or sp[3]=="ar":
                    sp[3] = 1.5
                if sp[3]=="am" or sp[3]=="Am":
                    sp[3] = 2
                bond_add = [int(sp[1]), int(sp[2]), int(sp[3])]
                bond_list.append(bond_add)
            elif flag[0] == 3 and len(sp) != 0:
                if sp[5][:1]=="H" or sp[5][:1]=="h":
                    continue
                atom_add = [int(sp[0]), sp[1], sp[5][:1]]
                atom_list.append(atom_add)
                dika_list.append([float(sp[2]),float(sp[3]),float(sp[4])])
        f.close()
    return atom_list, bond_list
    """, dika_list"""


if __name__ == '__main__':
    file = "./sp.mol2"
    atom, bond,dika = read_mol2(file)
    maxx1 = -10000
    maxy1 = -10000
    maxz1 = -10000
    minx1 = 10000
    miny1 = 10000
    minz1 = 10000
    for i in range(0,len(dika)):
        maxx1 = max(maxx1,dika[i][0])
        maxy1 = max(maxy1, dika[i][1])
        maxz1 = max(maxz1, dika[i][2])
        minx1 = min(minx1, dika[i][0])
        miny1 = min(miny1, dika[i][1])
        minz1 = min(minz1, dika[i][2])
    print(maxx1,minx1,maxy1,miny1,maxz1,minz1)
    print(maxx1-minx1,maxy1-miny1,maxz1-minz1)
    # print(atom)
    # print(bond)
