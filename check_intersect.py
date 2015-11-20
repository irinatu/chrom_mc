#check if there is intersect of polimer in the pdb file

import sys, optparse
import numpy as np


def pars_inp():
    ## read/validate command-line arguments
    optparser = optparse.OptionParser(usage = "%prog  [<options>]")

    optparser.add_option('-f', type = "string",
        dest = "Pdb_file",
        default = '',
        help = "The pdb file of trajectory")
			
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
            coord_mtx.append([float(line_sp[5]), float(line_sp[6]), float(line_sp[7])])
    inF.close()
    return coord_mtx_list

def distance (p1, p2):
    return np.sqrt(np.sum((p1-p2)**2))  
    
def cross(cha, po1, po2):
    pier = np.where((cha == po1).all(axis=1))
    dru = np.where((cha == po2).all(axis=1))
    #print "pir-dru ", len(pier), len(dru), pier, dru, pier[0][0]
    if len(pier) == 1 and len(dru) == 1:
        if pier[0] == dru[0]+1 or pier[0] == dru[0]-1:
            #print 'TRUE cross'
            return True
        else:
            #print 'FALSE cross'
            return False
    else: print "Two atoms of the chain have the same coordinates. Please check ", po1, po2, pier[0], dru[0]
    
def intersect(new_p, next, atC): #coordinates of two next polimer points
    #print dist(next, new_p), distance(next, new_p)
    if  distance(next, new_p) == 3: 
        return False
    elif distance(next, new_p) > 3:
        differ = new_p - next
        if differ[0] == 0:
            pos1 = np.array(([next[0], next[1], new_p[2]]))
            pos2 = np.array(([next[0], new_p[1], next[2]]))
            return cross(atC, pos1, pos2)
         
        elif differ[1] == 0:
            pos1 = np.array(([next[0], next[1], new_p[2]]))
            pos2 = np.array(([new_p[0], next[1], next[2]]))
            return cross(atC, pos1, pos2)

        elif differ[2] == 0:
            pos1 = np.array(([next[0], new_p[1], next[2]]))
            pos2 = np.array(([new_p[0], next[1], next[2]]))
            return cross(atC, pos1, pos2)

        else: print new_p, next, "The distance between these two positions are not = sqrt(18)"
        
#def check_pdb_steps(traj):
#    for pdb in traj:
#        pdb = numpy.array(pdb)
#        for atom in pdb:
      
    
opts = pars_inp()
coord_from_traj = extract_coord(opts.Pdb_file)
for pdb in coord_from_traj:
    for i in range (0, len(pdb)):
        if i != len(pdb) -1 and intersect(pdb[i], pdb[i+1], pdb):
            print "Intersect in the structure ", coord_from_traj.index(pdb), " between ", i, i+1
        else: pass
