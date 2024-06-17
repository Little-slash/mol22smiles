from pysmiles import write_smiles, fill_valence
import matplotlib.pyplot as plt
import read_mol2
import networkx as nx


# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

def linearco(molecule1, lian):
    fp = molecule1.copy()
    for k in range(0,num-1):
        spk = molecule1.copy()
        mapping = {node: node + max(fp.nodes) for node in spk.nodes}
        spk = nx.relabel_nodes(spk, mapping)
        G_combined = nx.compose(fp, spk)
        G_combined.add_edge(lian[0], lian[1] + max(molecule1.nodes), order=1.0)
        fp = G_combined.copy()

    return fp


def hunco(molecules, lian):   # 方形
    num = 2  # 聚合度


def paint_molecule(G):
    element_colors = {
        'C': 'lightblue',
        'O': 'orange',
        'N': 'green'
    }
    node_colors = [element_colors[G.nodes[node]['element'][:1]] for node in G.nodes]

    # 获取节点标签（显示element属性）
    node_labels = {node: f"{node}: {G.nodes[node]['element']}" for node in G.nodes}

    # 绘制图
    pos = nx.spring_layout(G)  # 布局
    nx.draw(G, pos, labels=node_labels, with_labels=True, node_size=500, node_color=node_colors, font_weight='bold')
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    filename = "./sp.mol2"
    atom_list, bond_list = read_mol2.read_mol2(filename)
    atom_index = [atom_list[sp][0] for sp in range(0, len(atom_list))]
    max1 = 0
    for i in range(0, len(atom_index)):
        max1 = max(atom_index[i], max1)
    hbond = [0 for i in range(0, max1 + 1)]
    molecule1 = nx.Graph()
    molecule1.add_nodes_from(atom_index)
    for i in range(0, len(bond_list)):
        if atom_index.count(bond_list[i][0]) == 0:
            hbond[bond_list[i][1]] += 1
        elif atom_index.count(bond_list[i][1]) == 0:
            hbond[bond_list[i][0]] += 1
        else:
            molecule1.add_edge(bond_list[i][0], bond_list[i][1], order=bond_list[i][2])

    for i in range(0, len(atom_list)):
        molecule1.nodes[atom_list[i][0]]['element'] = atom_list[i][1]
        molecule1.nodes[atom_list[i][0]]['hcount'] = hbond[atom_list[i][0]]
    num = 2  # 聚合度

    paint_molecule(molecule1)

    lian = [41, 7]
    keywords = input()

    if keywords == "linear":
        G_co = linearco(molecule1, lian)
    elif keywords == "hun":
        G_co = hunco(molecule1,lian)


    print(write_smiles(molecule1))
    print(write_smiles(G_co))
    sp = write_smiles(G_co)
    filename_out = "try2.smi"
    with open(filename_out, 'w') as file:
        file.write(write_smiles(G_co))

