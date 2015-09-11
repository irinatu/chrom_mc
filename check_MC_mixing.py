#read pdb file and give the contact (Hi-C) matrix

import sys, optparse
import numpy as np


def pars_inp():
    ## read/validate command-line arguments
    optparser = optparse.OptionParser(usage = "%prog -p structure.pdb [<options>]")

    optparser.add_option('-f', type = "string",
        dest = "First_traj",
        default = '',
        help = "The first file with trajectory in pdb format")
    optparser.add_option('-s', type = "string",
        dest = "Second_traj",
        default = '',
        help = "The second file with trajectory in pdb format")
			
    (opts, args) = optparser.parse_args()

    print len(sys.argv)
    if len(sys.argv) < 2:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)
    return opts	
                 
                 
def extract_coord(file):
    coord_mtx_list = []
    inF = open(file, 'r')
    coord_mtx = []
    for line in inF:
        if "HEADER" in line:
            if len(coord_mtx) > 0:
                coord = np.array(coord_mtx)
                coord_mtx_list.append(coord)
                coord_mtx = []
            else: coord_mtx = []
        elif line[0:4] == "ATOM" and line[13] == "C":
            line_sp = line.split()
            coord_mtx.append([float(line_sp[6]), float(line_sp[6]), float(line_sp[7])])
    return coord_mtx_list

    
def subset_prep(coor_list, mid, pol_l):
    intra_dist_list = []
    inter_mat = np.zeros((plo_l, pol_l), float)
    steps = [mid *k/20 for k in range(10,21)]
    for s in steps:
        for i in range(pol_l):
            for j in range(i, pol_l):
                dis = np.sqrt(np.sum((coor_list[s][i] - coord_list[s][j])**2))
                inter_mat[i,j] = inter_mat[j,i] = dis
        intra_dist_list.append(inter_mat)                

def aver_dist(sub1, sub2):
    dis_list = []
    for su in sub1:
        for us in sub2:
            dis_mat = np.sqrt(np.sum((su-us)**2)) 
            dis_list.append(dis_mtx)
    return np.mean(dis_list)           
    
    
    
opts = pars_inp()
coord_from_traj1 = extract_coord(opts.First_traj)
coord_from_traj2 = extract_coord(opts.Second_traj)
#print len(coord_from_traj1), len(coord_from_traj2)

if coord_from_traj1[0].shape[0] != coord_from_traj2[0].shape[0]:
    print "ERROR: The polimers lenght in the first and the second trajectories are unequal!!! "
    sys.exit(1)
else: pol_len = coord_from_traj1[0].shape[0]
 
if len(coord_from_traj1) == len(coord_from_traj2):
    middle = len(coord_from_traj1)/2
elif len(coord_from_traj1) < len(coord_from_traj2): 
    middle = len(coord_from_traj1)/2
else: middle = len(coord_from_traj2)/2

subset_1 = subset_prep(coord_from_traj1, middle, pol_len)
subset_2 = subset_prep(coord_from_traj2, middle, pol_len)
av1 = aver_dist(subset_1, subset_1)
av12 = aver_dist(subset_1, subset_2)
print av1, av12 


 