import sys, optparse
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
from matplotlib.colors import LogNorm, SymLogNorm

DISTANCE = 1.0

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
        default = 5,
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
    for i in range(coor_m.shape[0]):
        for j in range(i, coor_m.shape[0]):
            dis = np.sqrt(np.sum((coor_m[i] - coor_m[j])**2))
            #print "DISTANCE", dis
            if dis <= DISTANCE * 3.0:
                dis_mtx[i,j] = dis_mtx[j,i] = dis_mtx[j,i] + dis
                print i, j, dis
            else: 
                pass
                #print i,j, dis
    return dis_mtx
    
def extract_contacts(files, bou, start, end, step):
    #print files
    coord_mtx_list = []

    for file in files.split():
        print file
        count_str = 0
        inF = open(file, 'r')
        coord_mtx = []
        for line in inF:
            #print line
            if "HEADER" in line:
                count_str = count_str + step
                if count_str >=  start:
                    if len(coord_mtx) > 0:
                        coord = np.array(coord_mtx)
                        try:
                            inter_mat
                            inter_mat = calculate_dist_mtx(coord, inter_mat)
                        except NameError:
                            inter_mat = np.zeros((coord.shape[0], coord.shape[0]), float)
                            inter_mat = calculate_dist_mtx(coord, inter_mat)
                        coord_mtx = []
                    else: coord_mtx = []
            elif bou:
                if line[0:4] == "ATOM" and line[13] == "C" and line[17:20] != "UNB":
                    line_sp = line.split()
                    coord_mtx.append([float(line_sp[6]), float(line_sp[6]), float(line_sp[7])])
                else: 
                    pass 
                    #print line[17:20], "tak"
            else:
                if line[0:4] == "ATOM" and line[13] == "C":
                    line_sp = line.split()
                    coord_mtx.append([float(line_sp[6]), float(line_sp[6]), float(line_sp[7])])
    return inter_mat
    
def plot_hic(mt):
    fig = plt.figure()
    plt.imshow(mt,origin='lower',norm=LogNorm(), interpolation='nearest')
    plt.show()
    
opts = pars_inp()
bound = opts.Bound_atm
interactions = extract_contacts(opts.First_traj, bound, opts.Start, opts.End, opts.Step )
plot_hic(interactions)

