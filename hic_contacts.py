#read pdb file and give the contact (Hi-C) matrix

import sys


def pars_inp():
    ## read/validate command-line arguments
    optparser = optparse.OptionParser(usage = "%prog -p structure.pdb [<options>]")

    optparser.add_option('-p', type = "string",
        dest = "Int_pdb",
        default = '',
        help = "An input filename with 3D structure in pdb format")
			
    (opts, args) = optparser.parse_args()

    if len(args) < 1:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)
    return opts	
    
    
    
opts = pars_inp()