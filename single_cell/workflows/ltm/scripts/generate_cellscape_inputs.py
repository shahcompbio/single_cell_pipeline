import os
import sys
import networkx as nx
import pandas as pd
import argparse


def root_tree(G, root_id):

    G2=nx.Graph()
    for ed in G.edges():
        G2.add_edge(ed[0],ed[1])
    T=nx.dfs_tree(G2, root_id)
    return T


def generate_edges_list(gml_tree, root_id, edges_list_path, rooted_gml_tree):
    G=nx.read_gml(gml_tree)
    outfile=open(edges_list_path,'w')
    outfile.write('source,target'+'\n')
    rooted_tree = root_tree(G, root_id)
    nx.write_gml(rooted_tree, rooted_gml_tree)
    for e in rooted_tree.edges():
        outfile.write(str(e[0])+','+str(e[1])+'\n')
    outfile.close()


def generate_cn_data(merged_matrix, cn_data_path):
    bins_file=open(merged_matrix)
    outfile=open(cn_data_path,'w')
    outfile.write('"chr","start","end","copy_number","single_cell_id"'+'\n')

    # reading the bins list
    bins_list=[]
    l=bins_file.readlines()
    for i in range(1,len(l)):
        tmp=l[i].strip().split(',')
        bins_list.append(tmp[0]+','+tmp[1]+','+tmp[2])
    bins_file.close()

    cn_pd=pd.read_csv(merged_matrix)
    cell_ids=list(cn_pd)[4:]
    for i in range(len(cell_ids)):
        current_cell=cell_ids[i]
        for b in range(len(bins_list)):
            outfile.write(bins_list[b]+','+str(cn_pd[current_cell][b])+','+current_cell+'\n')
    outfile.close()


def generate_annotations(merged_matrix, annotations_path):
    outfile=open(annotations_path,'w')
    outfile.write('single_cell_id' + ',' + 'genotype' + '\n')
    infile=open(merged_matrix)
    cells=infile.readline().strip().split(',')[4:]
    for cell_id in cells:
        outfile.write(cell_id + ',' + '0' + '\n') # All cells annotated with 0
    outfile.close()


def main_generate_all(merged_matrix, annotations_path, edges_list_path, cn_data_path, gml_tree, rooted_gml_tree, root_id):
    generate_edges_list(gml_tree, root_id, edges_list_path, rooted_gml_tree)
    generate_annotations(merged_matrix, annotations_path)
    generate_cn_data(merged_matrix, cn_data_path)


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--path_to_data" ,type=str, help='path to data')
    parser.add_argument("--path_to_annotations" ,type=str, help='path to output annotations csv')
    parser.add_argument("--path_to_edges_list" ,type=str, help='path to output edges list csv')
    parser.add_argument("--path_to_cn_data" ,type=str, help='path to output cn data csv')
    parser.add_argument("-t","--path_to_tree" ,type=str, help='path to tree in gml format')
    parser.add_argument("--path_to_rooted_tree" ,type=str, help='path to rooted tree in gml format')
    parser.add_argument("-r","--root_id" ,type=str, help='id of the root node')
    args = parser.parse_args()

    main_generate_all(args.path_to_data, args.path_to_annotations, args.path_to_edges_list, args.path_to_cn_data, args.path_to_tree, 
        args.path_to_rooted_tree, args.root_id)
