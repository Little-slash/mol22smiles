from pysmiles import write_smiles, fill_valence
import matplotlib.pyplot as plt
import read_mol2  ## tips: the py-files in this catalogue
import networkx as nx


def linearco(molecule1, lian):
    fp = molecule1.copy()
    num = 2
    for k in range(0, num-1):
        spk = molecule1.copy()
        mapping = {node: node + max(fp.nodes) for node in spk.nodes}
        spk = nx.relabel_nodes(spk, mapping)
        G_combined = nx.compose(fp, spk)
        G_combined.add_edge(lian[0], lian[1] + max(molecule1.nodes), order=1.0)
        fp = G_combined.copy()

    return fp


def hunco(molecules, lian):   # 方形
    num = 2  # 聚合度(degree of polymerization,dp)


def paint_molecule(G):
    element_colors = {
        'C': 'lightblue',
        'O': 'orange',
        'N': 'green'
    }
    node_colors = [element_colors[G.nodes[node]['element'][:1]] for node in G.nodes]

    # 获取节点标签（显示element属性）
    node_labels = {node: f"{node}: {G.nodes[node]['element']}" for node in G.nodes}

    # plot
    pos = nx.spring_layout(G)  # 布局
    nx.draw(G, pos, labels=node_labels, with_labels=True, node_size=500, node_color=node_colors, font_weight='bold')
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
    plt.show()


def paint1(atom_info, atom_coords):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Extract atomic sequence number, element type, and coordinates
    atom_indices = [info[0] for info in atom_info]
    atom_types = [info[1] for info in atom_info]
    xs = [coord[0] for coord in atom_coords]
    ys = [coord[1] for coord in atom_coords]
    zs = [coord[2] for coord in atom_coords]

    # 定义元素种类对应的颜色
    color_map = {'C': 'k', 'N': 'b', 'O': 'r'}
    colors = [color_map[atom_type] for atom_type in atom_types]

    # 绘制原子
    ax.scatter(xs, ys, zs, c=colors, marker='o')

    # 绘制连线
    # 假设我们需要绘制所有原子之间的连线（可以根据需求调整连线的策略）

    # 设置标签
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # 添加原子标签（编号）
    for i in range(len(atom_info)):
        label = f"{atom_indices[i]}"
        ax.text(xs[i], ys[i], zs[i], label, fontsize=10, color='black')

    plt.show()

if __name__ == '__main__':
    filename = "./result.mol2"
    atom_list, bond_list, dika_list = read_mol2.read_mol2(filename)
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

    paint1(atom_list,dika_list)
    # paint_molecule(molecule1)

    lian = [45, 2]
    keywords = input()

    if keywords == "linear":
        G_co = linearco(molecule1, lian)
    # elif keywords == "hun":
    #     G_co = hunco(molecule1,lian)

    print(write_smiles(molecule1))
    print(write_smiles(G_co))
    sp = write_smiles(G_co)
    filename_out = "try2.smi"
    with open(filename_out, 'w') as file:
        file.write(write_smiles(G_co))

