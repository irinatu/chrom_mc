#read two trajectories written in pdb format and calculate mixing point based on  

import sys, optparse, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def pars_inp():
    ## read/validate command-line arguments
    optparser = optparse.OptionParser(usage = "%prog  [<options>]")

    optparser.add_option('-f', type = "string",
        dest = "First_traj",
        default = '',
        help = "The first file with trajectory in pdb format")
    optparser.add_option('-s', type = "string",
        dest = "Second_traj",
        default = '',
        help = "The second file with trajectory in pdb format")
    optparser.add_option('-b',
        dest = "Bound_atm",
        action = 'store_true',
        default = False,
        help = "Only BOU and LAM atoms will be considered")
    optparser.add_option('-c', type = "string",
        dest = "Chain",
        default = '',
        help = "The chain for mixing point")
			
    (opts, args) = optparser.parse_args()

#    print len(sys.argv)
    if len(sys.argv) < 2:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)
    return opts	
                 
                 
def extract_coord(files, bound, ch):
    print files
    coord_mtx_list = []
    for file in files.split():
        print file
        inF = open(file, 'r')
        coord_mtx = []
        for line in inF:
            if "HEADER" in line:
                if len(coord_mtx) > 0:
                    coord = np.array(coord_mtx)
                    coord_mtx_list.append(coord)
                    coord_mtx = []
                else: coord_mtx = []
            elif bound:
                #print line[21], ch, line
                if line[0:4] == "ATOM" and line[13] == "C" and line[17:20] != "UNB" and line[21] == ch:
                    #line_sp = line.split()
                    coord_mtx.append([float(line[30:37]), float(line[38:45]), float(line[46:53])])
                else: 
                    pass 
                    #print line[17:20], "tak"
            else:
                if line[0:4] == "ATOM" and line[13] == "C" and line[21] == ch:
                    #line_sp = line.split()
                    coord_mtx.append([float(line[30:37]), float(line[38:45]), float(line[46:53])])
        #print len(coord_mtx)
    return coord_mtx_list

    
def subset_prep(coor_list, mid, pol_l):
    intra_dist_list = []
    steps = [mid *k/20 for k in range(10,21)]
    for s in steps:
        inter_mat = np.zeros((pol_l, pol_l), float)
        for i in range(pol_l):
            for j in range(i, pol_l):
                dis = np.sqrt(np.sum((coor_list[s][i] - coor_list[s][j])**2))
                #print "DISTANCE", dis
                inter_mat[i,j] = inter_mat[j,i] = dis
        intra_dist_list.append(inter_mat)
    return intra_dist_list                

def aver_dist(sub1, sub2):
    dis_list = []
    for su in sub1:
        for us in sub2:
            dis_mat = np.sqrt(np.sum(np.power((su-us),2))) 
            #print dis_mat
            dis_list.append(dis_mat)
    return np.mean(dis_list)           
    
    
    
opts = pars_inp()
bound_at = opts.Bound_atm
coord_from_traj1 = extract_coord(opts.First_traj, bound_at, opts.Chain)
coord_from_traj2 = extract_coord(opts.Second_traj, bound_at, opts.Chain)
print "klatki", len(coord_from_traj1), len(coord_from_traj2)

if coord_from_traj1[0].shape[0] != coord_from_traj2[0].shape[0]:
    print "ERROR: The polimers lenght in the first and the second trajectories are unequal!!! ", coord_from_traj1[0].shape[0], coord_from_traj2[0].shape[0] 
    sys.exit(1)
else: pol_len = coord_from_traj1[0].shape[0]
 
if len(coord_from_traj1) == len(coord_from_traj2):
    middle = len(coord_from_traj1)/2
elif len(coord_from_traj1) < len(coord_from_traj2): 
    middle = len(coord_from_traj1)/2
else: middle = len(coord_from_traj2)/2
#incr = (middle/2) -1
incr = middle/5
whole = middle*2
#middle = 100000

print len(coord_from_traj1[0]), len(coord_from_traj2[0])
l_av1 = []
l_av12 = []
while middle < whole:
    subset_1 = subset_prep(coord_from_traj1, middle, pol_len)
    subset_2 = subset_prep(coord_from_traj2, middle, pol_len)
    print pol_len, middle, len(subset_1), len(subset_2)
    av1 = aver_dist(subset_1, subset_1)
    av12 = aver_dist(subset_1, subset_2)
    l_av1.append(av1)
    l_av12.append(av12)
    per = (av12-av1)*100.0/av1 
    print av1, av12, per, '%'
    if av12 > av1-av1*0.1 and av12 < av1+av1*0.1:
        print middle, "is enought steps", av1, av12
        break
    else:
        middle = middle + incr

fig = plt.figure()
p_one = plt.plot( l_av1, "b-", linewidth =0.5, label="Inside")
p_two = plt.plot( l_av12, "r-", linewidth =0.5, label="Between")
plt.axis([0,len(l_av1), 0, max(l_av12)])
plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
#plt.show()
head1, tail1 = os.path.split(opts.First_traj)
head2, tail2 = os.path.split(opts.Second_traj)
fig.savefig(tail1.split('.')[0]+ "_"+tail2.split('.')[0]+"_dist.png")
