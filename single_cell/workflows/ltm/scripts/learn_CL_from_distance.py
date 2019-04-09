import networkx as nx


def learn_CL(cell_data):
    node_ids = cell_data.keys()
    G = nx.Graph()

    for i in range(len(node_ids)):
        G.add_node(node_ids[i])
        for j in range(i+1,len(node_ids)):
            G.add_edge(node_ids[i], node_ids[j], weight=calc_MI(cell_data[node_ids[i]],cell_data[node_ids[j]]))
       
    T = nx.minimum_spanning_tree(G)
    return T


def learn_CL_from_distance(d_files, tree_path):
    G = nx.Graph()

    for file in d_files:
        infile=open(file)
        for line in infile:
            tmp=line.strip().split(',')
            G.add_edge(tmp[0],tmp[1], weight=float(tmp[2]))
        infile.close()

    T = nx.minimum_spanning_tree(G)
    nx.write_gml(T, tree_path)
    
    
