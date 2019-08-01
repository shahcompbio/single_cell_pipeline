import numpy as np
import networkx as nx
import timeit
from sklearn.metrics import mutual_info_score
import sklearn
import argparse
import string
import pandas as pd


def convert_n_alpha_cn(inar):
    num2alpha = dict(zip(range(0, 27), string.ascii_lowercase))
    
    s=''
    for n in inar:
        t=int(n)
        s+=num2alpha[t]
    return s


def breakpoint_conevert(bps_list, pre_cns):
    res_list=[]
    for i in range(len(pre_cns)):
        if i in bps_list:
            res_list.append(pre_cns[i])
    return res_list


def calc_MI(x, y):
    bins=len(x)/10
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return -mi


def calculate_distance(cell_data, infile_path, outfile_path):
    outfile=open(outfile_path,'w')
    pairs=set([])
    infile=open(infile_path)
    for line in infile:
        tmp=line.strip().split(',')
        t1=timeit.default_timer()
        weight = calc_MI(cell_data[tmp[0]],cell_data[tmp[1]])
        t2=timeit.default_timer()
        outfile.write(tmp[0]+','+tmp[1]+','+str(weight)+','+str(t2-t1)+'\n')
    outfile.close()
    infile.close()


def break_points_finder(infile_path):
    data = pd.read_hdf(data_path)

    # find breakpoint
    bps_set=set([])
    for j in range(data.shape[1]):
        for i in range(1,data.shape[0]):
            if data.iloc[i,j]!=data.iloc[i-1,j]:
                bps_set.add(i)
    bps_set.add(0)
    bps_list=list(bps_set)

    return bps_list


def with_break_point_correction_read_input_np(data_path, infile_path):
   
    bps_list= break_points_finder(data_path)

    cell_data={}
    
    ###### read list of filtered cells
    cells_set=set([])
    cells_file=open(infile_path)
    for line in cells_file:
        tmp=line.strip().split(',')
        cells_set.add(tmp[0])
        cells_set.add(tmp[1])
    
    #######
    
    data = pd.read_hdf(data_path)
    cell_ids = data.columns.values.tolist()[4:]

    for cell_id in cell_ids:
        #cell_data[i]=convert_n_alpha_cn(data[:,i])
        j=cell_ids.index(cell_id)
        if cell_id in cells_set:
            #cell_data[i]=breakpoint_conevert(bps_list, data[:,j])
            cell_data[cell_id]=breakpoint_conevert(bps_list, data.iloc[:,j+4])

    return cell_data


def read_input_np(data_path, infile_path):

    cell_data={}
    
    ###### read list of filtered cells
    cells_set=set([])
    cells_file=open(infile_path)
    for line in cells_file:
        tmp=line.strip().split(',')
        cells_set.add(tmp[0])
        cells_set.add(tmp[1])
    
    #######
    
    data = pd.read_csv(data_path)
    cell_ids = data.columns.values.tolist()[4:]

    for cell_id in cell_ids:
        j=cell_ids.index(cell_id)
        if cell_id in cells_set:
            cell_data[cell_id]=data.iloc[:,j+4]

    return cell_data


def read_hmm_cn_data(data_path, outfile_path, infile_path):
    
    cell_data = read_input_np(data_path, infile_path)

    CL_start = timeit.default_timer()
    calculate_distance(cell_data, infile_path, outfile_path)
    CL_end = timeit.default_timer()
    
    print ("total time "+str(CL_end-CL_start))
    

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_file", help="path to the edges lists folder.")
    parser.add_argument("-output_path", help="Path to the results.")
    parser.add_argument("-data_path", help="Path to the cn data.")

    args = parser.parse_args()
    input_file=args.input_file
    output_path=args.output_path
    data_path=args.data_path
    
    read_hmm_cn_data(data_path, output_path, input_file)
    

    