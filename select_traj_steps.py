## for giving trajectory/ies chose every defined number of frames
import numpy as np
import os, argparse, sys

def pars_inp():
    ## read/validate command-line arguments
    argparser = argparse.ArgumentParser( usage = "%prog  [<options>]")

    argparser.add_argument('-i',
        dest = "Inp",
        type=str,
        nargs='+',
        help = "input file/s with pdb trajectories")
    argparser.add_argument('-o', type = str,
        dest = "Out",
        default = '',
        help = "A output file with pdb trajectory")
    argparser.add_argument('-s', type = int,
        dest = "Skipp",
        default = 50,
        help = "Skipping frames in trajectory/ies. The smallest number is 2 - no skipping")
        
    opts = argparser.parse_args()
    #if len(args) < 2:
    #print opts
     #   print optparser.format_help() #prints help if no arguments
     #   sys.exit(1)
    return opts	

def find_lines_in_model(my_file, out): #find first HEADER and calculate number of positions for one model
    find = 0
    while 1:
        line = my_file.readline()
        if "HEADER" in line:
            pos_in_line = len(line)
            find += 1
            if find ==2: break
        out.write(line)
    pos_for_model = my_file.tell()-pos_in_line
    return pos_for_model
    
def write_model(m_file, out_f):
    find = 0
    while 1:
        line = m_file.readline()
        if "HEADER" in line:
            find += 1
            if find == 2: break
        out_f.write(line)

opts = pars_inp()
#print opts.Skipp
N = opts.Skipp #number of skipping models    
traj = open(opts.Inp[0])
#print opts.Inp
out_traj = open(opts.Out,"w")
lines_for_model  = find_lines_in_model(traj, out_traj)
print lines_for_model
for plik in opts.Inp:
    tra =open(plik)
    file_size = os.fstat(tra.fileno()).st_size
    while 1:
        tra.seek((N-2)*lines_for_model, 1)
        print tra.tell(), file_size
        if (tra.tell() + (N-2)*lines_for_model) >= file_size: #zeby za daleko nie skoczyc
            tra.close()
            break
        else: write_model(tra, out_traj)
out_traj.close()
 

