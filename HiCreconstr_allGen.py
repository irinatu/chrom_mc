import sys, optparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns
from matplotlib.colors import LogNorm, SymLogNorm

DISTANCE = 5.0

def pars_inp():
    ## read/validate command-line arguments
    optparser = optparse.OptionParser(usage = "%prog  [<options>]")

    optparser.add_option('-f', type = "string",
        dest = "First_traj",
        default = '',
        help = "The first file with trajectory in pdb format")
    #optparser.add_option('-s', type = "string",
     #   dest = "Second_traj",
      #  default = '',
       # help = "The second file with trajectory in pdb format")
    optparser.add_option('-b',
        dest = "Bound_atm",
        action = 'store_true',
        default = False,
        help = "Only BOU and LAM atoms will be considered")
    optparser.add_option('-s',
        dest = "Start",
        type = "int",
        default = 1,
        help = "The first step to make a map")
    optparser.add_option('-l',
        dest = "End",
        type = "int",
        default = 10,
        help = "The last step to make a map")
    optparser.add_option('-t',
        dest = "Step",
        type = "int",
        default = 1,
        help = "The step size to collect structures and make a HiC map")
    
			
    (opts, args) = optparser.parse_args()

#    print len(sys.argv)
    if len(sys.argv) < 2:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)
    return opts	

def calculate_dist_mtx(coor_m, dis_mtx):
    #print coor_m, dis_mtx
    for i in range(coor_m.shape[0]):
        for j in range(i, coor_m.shape[0]):
            dis = np.sqrt(np.sum((coor_m[i] - coor_m[j])**2))
            #print "DISTANCE", dis
          #  if dis <= DISTANCE * 3.0:
          #      dis = dis/3.0
          #      if dis == 0.0: 
          #          dis = 1.0 
          #          if i != j: print i,j, dis, coor_m[i], coor_m[j]
          #      dis_mtx[i,j] = dis_mtx[j,i] = dis_mtx[j,i] + 1.0
          #      #print i, j, dis
          #  else: 
          #      pass
                #print i,j, dis
            dis = dis/3.0
            if dis == 0.0: 
                dis = 1.0 
                if i != j: print i,j, dis, coor_m[i], coor_m[j]
            dis_mtx[i,j] = dis_mtx[j,i] = dis_mtx[j,i] + 1.0/dis # contact calculation
                #print i, j, dis
    return dis_mtx
    
def extract_contacts(files, bou, start, end, step):
    #print files
    coord_mtx_list = []
    count_str = 0
    steps = start

    for file in files.split():
        print file
        inF = open(file, 'r')
        chains = {}  # dictionary with dist_mtx for each chain {'A':dist_mtx, 'B':dist_mtx...}
        coord_mtx = [] # list of [x,y,z] for one chain [[x1,y1,z1],[x2,y2,z2]...]
        for line in inF:
            #print line
            if "HEADER" in line:
                #print line
                count_str = count_str + 1
                COOR = False
                if count_str >=  start and count_str <= end and count_str == steps:
                    COOR = True
                    #print start, end, steps, count_str
                    print line, steps
                    steps += step
                elif count_str >= end: 
                    break 
                else: coord_mtx = []
            elif line[0:3] == "TER":
                if len(coord_mtx) > 0:
                    coord = np.array(coord_mtx)
                    coord_mtx = []
                    #print coord.shape, coord
                    try:
                        chains[chain]
                        chains[chain] = calculate_dist_mtx(coord, chains[chain])
                    except KeyError:
                        chains[chain] = np.zeros((coord.shape[0], coord.shape[0]), float)
                        chains[chain] = calculate_dist_mtx(coord, chains[chain])
                else: coord_mtx = []
                
            elif bou and COOR:
                #print 'DODAJE'
                if line[0:4] == "ATOM" and line[13] == "C" and line[17:20] != "UNB":
                    #line_sp = line.split()
                    chain = line[21]
                    #print line
                    coord_mtx.append([float(line[30:37]), float(line[38:45]), float(line[46:53])])
                else: 
                    pass 
                    #print line[17:20], "tak"
            elif COOR:
                #print 'dodaje'
                if line[0:4] == "ATOM" and line[13] == "C":
                    chain = line[21]
                    line_sp = line.split()
                    #print line
                    coord_mtx.append([float(line[30:37]), float(line[38:45]), float(line[46:53])])
    return chains
    
def plot_hic(mt, cha):
    fig = plt.figure()
    plt.imshow(mt,origin='lower',norm=LogNorm(), interpolation='nearest')
    plt.colorbar()
    plt.show()

    fig.savefig("HiC_%s.png" %(cha))
    
opts = pars_inp()
bound = opts.Bound_atm
interactions_for_chains = extract_contacts(opts.First_traj, bound, opts.Start, opts.End, opts.Step )
for ch in interactions_for_chains:
    plot_hic(interactions_for_chains[ch], ch)

