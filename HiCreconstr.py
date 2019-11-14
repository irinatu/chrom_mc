import sys, optparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib import  cm
from util import clip_and_blur

DISTANCE = 5.0

def pars_inp():
    ## read/validate command-line arguments
    optparser = optparse.OptionParser(usage = "%prog  [<options>]")

    optparser.add_option('-p', type = "string",
        dest = "First_traj",
        default = '',
        help = "The trajectory in pdb format")
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
    optparser.add_option('-m',
        dest = "Mode",
        type = "int",
        default = 1,
        help = "The mode for work (1 or 2): for 1/dist use 1, for cut_off - use 2. 1 is default")
    optparser.add_option('-c',
        dest = "CutOff",
        type = "float",
        default = 3.5,
        help = "TCutOff for calculation of contac in 'cut_aff mode', default is 3.5")
    
			
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

def calculate_dist_mtx_cutoff(coor_m, dis_mtx):
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
            if dis <= opts.CutOff and i != j: # calculate only contacts within 3.5 cutoff
                dis_mtx[i,j] = dis_mtx[j,i] = dis_mtx[i,j] + 1.0 # contact calculation
                #print i, j, dis
        dis_mtx[i,i] = dis_mtx[i,i] +1.0
    #print dis_mtx[5:25, 5:25], dis_mtx.min(), dis_mtx.max()
    return dis_mtx
    
def extract_contacts(files, bou, start, end, step):
    #print files
    coord_mtx_list = []
    count_str = 0
    steps = start
    file_nr = 0.0

    for file in files.split():
        print file
        file_nr = file_nr + 1.0
        inF = open(file, 'r')
        coord_mtx = []
        for line in inF:
            #print line
            if "MODEL" in line:
                #print line
                count_str = count_str + 1
                COOR = False
                #print start, end, steps
                if count_str >=  start and count_str <= end and count_str == steps:
                    COOR = True
                    #print start, end, steps, count_str
                    print line, steps
                    steps += step
                elif count_str >= end: 
                    break 
                else: coord_mtx = []
            elif line[0:6] == "ENDMDL":
                if len(coord_mtx) > 0:
                    coord = np.array(coord_mtx)
                    coord_mtx = []
                    #print coord.shape, coord
                    try:
                        inter_mat
                        if opts.Mode ==1:
                            inter_mat = calculate_dist_mtx(coord, inter_mat)
                        elif opts.Mode ==2:
                            inter_mat = calculate_dist_mtx_cutoff(coord, inter_mat)
                    except NameError:
                        inter_mat = np.zeros((coord.shape[0], coord.shape[0]), float)
                        if opts.Mode ==1:
                            inter_mat = calculate_dist_mtx(coord, inter_mat)
                        elif opts.Mode ==2:
                            inter_mat = calculate_dist_mtx_cutoff(coord, inter_mat)
                else: coord_mtx = []
                
            elif bou and COOR:
                #print 'DODAJE'
                if line[0:4] == "ATOM" and line[13] == "C" and line[17:20] != "UNB":
                    line_sp = line.split()
                    #print line
                    coord_mtx.append([float(line_sp[5]), float(line_sp[6]), float(line_sp[7])])
                else: 
                    pass 
                    #print line[17:20], "tak"
            elif COOR:
                #print 'dodaje'
                if line[0:4] == "ATOM" and line[13] == "C":
                    line_sp = line.split()
                    #print line
                    coord_mtx.append([float(line_sp[5]), float(line_sp[6]), float(line_sp[7])])
    #inter_mat = inter_mat/ file_nr
    return inter_mat
    
def plot_hic(mt):
    fig = plt.figure()
    #print mt[5:25, 5:25]
    if opts.Mode ==1:
        mt = clip_and_blur(mt)
    plt.imshow(mt,origin='lower',norm=LogNorm(), cmap=cm.jet, interpolation='nearest')
    #plt.imshow(mt,origin='lower', cmap=cm.jet)
    plt.colorbar()
    #plt.show()

    fig.savefig("HiC.png", dpi=1000)
    
opts = pars_inp()
bound = opts.Bound_atm
interactions = extract_contacts(opts.First_traj, bound, opts.Start, opts.End, opts.Step )
np.save("HiC.npy", interactions)
plot_hic(interactions)

