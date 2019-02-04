#read two trajectories written in pdb format and calculate mixing point based on  

import sys, optparse, os, re
import numpy as np
import matplotlib
from collections import Counter
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
        default = ' ',
        help = "The chain for mixing point, default ' '")
			
    (opts, args) = optparser.parse_args()


    if len(sys.argv) < 2:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)
    return opts	
                 
                 
def extract_coord(files, bound, ch, steps_subset):
    #print files
    coord_mtx_list = []
    str_nr = 0
    for file in files.split():
        #print file
        inF = open(file, 'r')
        coord_mtx = []
        for line in inF:
            if "HEADER" in line:
                str_nr = str_nr +1
                continue
                #if str_nr in steps_subset:
                #    print "STR_NR", str_nr
                   
                #else: continue
            elif str_nr in steps_subset:
                if bound:
                    #print line[0:4], line[13], line[17:20], line[21], ch
                    if line[0:4] == "ATOM" and line[13] == "C" and line[17:20] != "UNB" and line[21] == ch:
                        coord_mtx.append([float(line[30:37]), float(line[38:45]), float(line[46:53])])
                    else: 
                        pass 

                else:
                    if line[0:4] == "ATOM" and line[13] == "C" and line[21] == ch:
                        #line_sp = line.split()
                        coord_mtx.append([float(line[30:37]), float(line[38:45]), float(line[46:53])])
            elif len(coord_mtx) > 0:
                coord = np.array(coord_mtx)
                coord_mtx_list.append(coord)
                coord_mtx = []
            else:continue
        #print len(coord_mtx)
    return coord_mtx_list

    
def subset_prep(coor_list, pol_l):
    intra_dist_list = []
    for macie in coor_list:
        inter_mat = np.zeros((pol_l, pol_l), float)
        for i in range(pol_l):
            for j in range(i, pol_l):
                dis = np.sqrt(np.sum((macie[i] - macie[j])**2))
                #print "DISTANCE", dis
                inter_mat[i,j] = inter_mat[j,i] = dis
        intra_dist_list.append(inter_mat)
    return intra_dist_list                              # list of martices with distances between atoms             

def aver_dist(sub1, sub2):
    dis_list = []
    for su in sub1:
        for us in sub2:
            dis_mat = np.sqrt(np.sum(np.power((su-us),2))) 
            #print dis_mat
            dis_list.append(dis_mat)
    return np.mean(dis_list)                  # mean of distances between structures in two subsets
    


def calc_frames(files):
    print files
    whole = 0
    for one in files.split():
        f = open(one)
        passage = f.read()
        words = re.findall(r'HEADER', passage)
        #word_counts = Counter(words)
        whole = whole + len(words)
    print whole
    return whole

    
opts = pars_inp()
bound_at = opts.Bound_atm
whole_nr1 = calc_frames(opts.First_traj)
whole_nr2 = calc_frames(opts.Second_traj)


middle = min(whole_nr1, whole_nr2)/2
whole = middle*2
incr = min(middle/5, 3000)

l_av1 = []
l_av12 = []
while middle < whole:
    steps = [middle *k/20 for k in range(10,21)]
    print "STEPS", steps
    coord_from_traj1 = extract_coord(opts.First_traj, bound_at, opts.Chain, steps) #list of matrices of coordinates
    coord_from_traj2 = extract_coord(opts.Second_traj, bound_at, opts.Chain, steps)
    
    if coord_from_traj1[0].shape[0] != coord_from_traj2[0].shape[0]:
        print "ERROR: The polimers lenght in the first and the second trajectories are unequal!!! ", coord_from_traj1[0].shape[0], coord_from_traj2[0].shape[0] 
        sys.exit(1)
    else: pol_len = coord_from_traj1[0].shape[0]
    
    subset_1 = subset_prep(coord_from_traj1, pol_len)
    subset_2 = subset_prep(coord_from_traj2, pol_len)
    print pol_len, middle, whole, len(subset_1), len(subset_2)
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
