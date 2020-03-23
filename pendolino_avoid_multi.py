import numpy,math,time, optparse, sys, pickle, json
import random as random

#accepted move vectors
#MOVES = numpy.array([[0,0,1], [0,1,0], [1,0,0], [0,0,-1], [0,-1,0], [-1,0,0], [1,0,1], [0,1,1], [1,1,0], [-1,0,-1], [0,-1,-1], [-1,-1,0], [1,0,-1], [0,1,-1], [1,-1,0], [-1,0,1], [0,-1,1], [-1,1,0], [-1,1,-1], [-1,-1,-1], [1,1,1], [1,-1,1], [1,-1,-1], [-1,-1,1], [-1,1,1], [1,1,-1]])
MOVES = numpy.array([[0,0,1], [0,1,0], [1,0,0], [0,0,-1], [0,-1,0], [-1,0,0], [1,0,1], [0,1,1], [1,1,0], [-1,0,-1], [0,-1,-1], [-1,-1,0], [1,0,-1], [0,1,-1], [1,-1,0], [-1,0,1], [0,-1,1], [-1,1,0]])
#accepted matching positions of binding sites
BMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]])
#BMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1], [1,1,0], [-1,-1,0], [-1,1,0], [1,-1,0], [1,0,1], [-1,0,1], [1,0,-1], [-1, 0, -1], [0,1,-1], [0,-1,-1], [0,-1,1], [0,1,1]])

# radius of the nucleus
R = 20
# 2 x radius + a fringe, because lamin barrier has to be hermetic
BOUND = 2 * R + 2
#ran_seed = 2
#random.seed(ran_seed)
#print ran_seed

EMPTY = 0
LAMIN = 1
BSITE_L = 2
REGDNA = 3
BSITE_R = []
BINDER = []

DIST = 3

# max distance of good neighbors
GOOD_NEIGH = 3

def pars_inp():
    ## read/validate command-line arguments
    optparser = optparse.OptionParser(usage = "%prog regular_bsites.txt lamin_bsites.txt [<options>]")

    optparser.add_option('-p', type = "string",
        dest = "Out_str",
        default = '',
        help = "An output filename with MC trajectory")
    optparser.add_option('-i', type = "string",
        dest = "In_str",
        default = '',
        help = "An input state (pickled file)")
    optparser.add_option('-s', type = "int",
        dest = "Steps",
        default = 100000,
        help = "Number of steps for simulation (default 100000)")
    optparser.add_option('-l', type = "int",
        dest = "Ch_lenght",
        default = 512,
        help = "Lenght of the chain (default 512)")
    optparser.add_option('-b', type = "string",
        dest = "Binders",
        default = 256,
        help = "Number of binders (default 256)")
    optparser.add_option('-e', type = "int",
        dest = "Save",
        default = 100,
        help = "Save every .... accepted step (default 100")
    optparser.add_option('-n', type = "string",
        dest = "Regular_bsites",
        default = "",
        help = "Full path to the file names with regular binding sites, separated by comma. Mandatory!")
    optparser.add_option('-a', type = "string",
        dest = "Lamin_bsites",
        default = "",
        help = "Full path to the file names with lamin binding sites, separated by comma. Mandatory!")
			
    (opts, args) = optparser.parse_args()

    #if len(args) < 2 and opts.In_str =='':
    #    print optparser.format_help() #prints help if no arguments
    #    sys.exit(1)
    return opts	

def init_dist_matrix(max_d = GOOD_NEIGH + 1):
    dist_matrix = numpy.zeros((max_d, max_d, max_d), dtype = numpy.float32)
    for i in range(max_d):
        for j in range(max_d):
            for k in range(max_d):
                dist_matrix[i][j][k] = math.sqrt(i**2 + j**2 + k**2)
    return dist_matrix
DIST_MATRIX = init_dist_matrix()

def initialize_import(f):
    list_ob = pickle.load(open(f))
    ch = list_ob[0]
    b = list_ob[1]
    a = list_ob[2]
    state = list_ob[3]
    b_nr = list_ob[4]
    b_sites = list_ob[5]
    binder_l = list_ob[6]
    bsite_binder = list_ob[7]
    return ch, b, a, state, b_nr, b_sites, binder_l, bsite_binder
    
def initialize_import_json(f):
    list_ob = json.load(open(f))
    ch = numpy.asarray(list_ob[0])
    b = numpy.asarray(list_ob[1])
    a = list_ob[2]
    state = numpy.asarray(list_ob[3])
    b_nr = list_ob[4]
    b_sites = list_ob[5]
    binder_l = list_ob[6]
    bsite_binder = list_ob[7]
    print bsite_binder
    bsite_binder = {int(k):v for k,v in bsite_binder.items()}
    print bsite_binder
    return ch, b, a, state, b_nr, b_sites, binder_l, bsite_binder

def dist_from_mi(x, y, z, mi):
        return math.sqrt((x - mi)**2 + (y - mi)**2 + (z - mi)**2)

def getStateWithLamins(bound, f):

    state = numpy.zeros((bound, bound, bound), dtype = numpy.int)
    MIDDLE = bound / 2
    lam_name = f.split('.pdb')[0] + '_lamin.pdb'
    save_lam = open(lam_name, "w")
    save_lam.write("MODEL LAMINA")
    at_nr = 1

    

    for x in range(BOUND):
        for y in range(BOUND):
            for z in range(BOUND):
                border = abs(dist_from_mi(x, y, z, MIDDLE) - MIDDLE + 1)
                if border <= 2:
                    state[x, y, z] = LAMIN
                    if border == 1:
                        #print border
                        line = "\nATOM  " + str(at_nr).rjust(5) + " " + "P".center(4) + " " + "LAM" + "  " + str(at_nr).rjust(4) + "    " + str(round(x * DIST, 3)).rjust(8) + str(round(y * DIST, 3)).rjust(8) + str(round(z * DIST, 3)).rjust(8) + "  0.00 00.00"
                        at_nr += 1
                        save_lam.write(line)
    save_lam.write("\nENDMDL")
    save_lam.close()
    return state


def dist(r1,r2):
    return DIST_MATRIX[tuple(abs(r1 - r2))]
    
#def distance (p1, p2):
 #   return numpy.sqrt(numpy.sum((p1-p2)**2))

def cross(cha, po1, po2):
    pier = numpy.where((cha == po1).all(axis=1))
    dru = numpy.where((cha == po2).all(axis=1))
    #print "pir-dru ", len(pier), len(dru), pier, dru, pier[0][0]
    if len(pier) == 1 and len(dru) == 1:
        if pier[0] == dru[0]+1 or pier[0] == dru[0]-1:
            #print 'TRUE cross'
            return True
        else:
            #print 'FALSE cross'
            return False
    else: print "Two atoms of the chain have the same coordinates. Please check ", po1, po2, pier[0], dru[0]
    
def intersect(new_p, next, sta, ch): #coordinates of two next polimer points
    #import functools, operator # to reduce list of list to the flat list
    #print dist(next, new_p), distance(next, new_p)
    def flatten(lst):
	    return sum( ([x] if not isinstance(x, list) else flatten(x) for x in lst), [] )
    
    #BINDER_flat = list(set(flatten(BINDER_dic.values()))) # uniqe values of BINDER 
    if  dist(next, new_p) == 1: 
        return False	    
    elif dist(next, new_p) > 1:
        differ = new_p - next
        #print "TO", BINDER, flatten(BINDER)
        if differ[0] == 0:
            pos1 = numpy.array(([next[0], next[1], new_p[2]]))
            pos2 = numpy.array(([next[0], new_p[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN] + BINDER  and sta[tuple(pos2)] not in [EMPTY, LAMIN] + BINDER :
                return cross(ch, pos1, pos2)
            else: return False         
        elif differ[1] == 0:
            pos1 = numpy.array(([next[0], next[1], new_p[2]]))
            pos2 = numpy.array(([new_p[0], next[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN] + BINDER  and sta[tuple(pos2)] not in [EMPTY, LAMIN] + BINDER :
                return cross(ch, pos1, pos2)
            else: return False
        elif differ[2] == 0:
            pos1 = numpy.array(([next[0], new_p[1], next[2]]))
            pos2 = numpy.array(([new_p[0], next[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN] + BINDER  and sta[tuple(pos2)] not in [EMPTY, LAMIN] + BINDER :
                return cross(ch, pos1, pos2)
            else: return False
        else: print new_p, next, "The distance between these two positions are not = sqrt(2)"

def initialize_random(n, m, fa, bound = BOUND):
    chain = numpy.zeros((n, 3), dtype = numpy.int)
    binders = numpy.zeros((sum(m), 3), dtype = numpy.int)
    state = getStateWithLamins(bound, fa)
    attached_to_lamins = []

    chain[0] = [bound / 2] * 3

    def get_site_type_list(fpath, arra, site_type):
        #print "lent", length_list, fpath
        for l in open(fpath):
            arra[int(l) -1][site_type] = 1
        return arra

    
    regular_bsites = numpy.zeros((n, len(opts.Regular_bsites.split(",")))) ## 3D array, columns - different binding types, each array for different binding sites
    for porz, file_regular in enumerate(opts.Regular_bsites.split(",")):
        regular_bsites = get_site_type_list(file_regular, regular_bsites, porz)
    lamin_bsites = numpy.zeros((n, 1))
    if opts.Lamin_bsites != "":
        lamin_bsites   = get_site_type_list(opts.Lamin_bsites, lamin_bsites, 0)
    else: pass
    
    #regular_bsites = get_site_type_list(sys.argv[1], n)
    #lamin_bsites   = get_site_type_list(sys.argv[2], n)

    def get_site_type(i, regular_bsites, lamin_bsites): # BSITE_R interacts with binders whereas BSITE_L interacts both with lamins and binders
        #b = [el[i] for el in regular_bsites] # list with values of ith bin for each binders, len of list = numer of binders
        #indices = [k for k, x in enumerate(b) if x == 1] # indexes of b list with 1 values
        indices =  numpy.where(regular_bsites[i] == 1) # regular_bsites - 3D array with layer nr = nr of binding_sites_types, rows - atom_number, columns - chromosomes
        if len(indices[0]) > 1: 
            print "Warning!!! The  %i atom is assign to the different binders. " %(i+1)
            ke = tuple([ BSITE_R[index] for index in indices[0]]) # key of the reverse dictionary, to find number for state item for  the atom that can bind more than one binders
            val = BSITE_dic_re[ke]
            #print indices, indices[0], ke, val
            #sys.exit(1)
            return val
            #return BSITE_R[BINDER_dic.keys().index(val)]
            #return BSITE_R[random.choice(indices[0])]
        if lamin_bsites[i,0] == 1:
            return BSITE_L
        elif len(indices[0]) == 1:
            return BSITE_R[indices[0][0]]
        else:
            return REGDNA


    cur = chain[0]
    #if opts.Lamin_bsites != "":
    state[tuple(cur)] = get_site_type(0, regular_bsites, lamin_bsites)
    #else:
    #    state[tuple(cur)] = get_site_type(0, regular_bsites)
        
    for i in range(1, n):
        mov = random.choice(MOVES)
        tries = 0
        while tries < 100 and (not (no_collisions(tuple(cur + mov), state)) or intersect(cur, cur+mov, state, chain)):
            mov = random.choice(MOVES)
            tries += 1    
            
        assert tries != 100, "unable to find initialization"
        chain[i] = cur + mov
        #if opts.Lamin_bsites != "":
        state[tuple(chain[i])] = get_site_type(i, regular_bsites, lamin_bsites)
        #else:
        #    state[tuple(cur)] = get_site_type(0, regular_bsites)

        if state[tuple(chain[i])] == BSITE_L and count_bonds(chain[i], [LAMIN], state) > 0:
            attached_to_lamins.append(tuple(chain[i]))

        cur = chain[i]
        
    mid = bound/2
    at_bin = -1
    for bind, b_nr in zip(m, BINDER):
        for i in range(bind):
            at_bin += 1
            x = random.randint(0, bound)
            y = random.randint(0, bound)
            z = random.randint(0, bound)
            tries = 0
            distance = dist_from_mi(x, y, z, mid)
            while distance > mid-3 or (not (no_collisions((x, y, z), state)) and tries < 100):
                x = random.randint(0, bound)
                y = random.randint(0, bound)
                z = random.randint(0, bound)
                tries += 1
                distance = dist_from_mi(x, y, z, mid)
            binders[at_bin] = [x, y, z]
            state[tuple(binders[at_bin])] = b_nr
    
    

    return chain, binders, attached_to_lamins, state


def good_neighbors(x, i, chain):
    #check neighbors (chain)
    if i > 0:
        d1 = dist(chain[i - 1], x)
        if d1 > 0 and d1 < GOOD_NEIGH:
            pass
        else:
            return False
    if i < len(chain) - 1:
        d2 = dist(chain[i + 1], x)
        if d2 > 0 and d2 < GOOD_NEIGH:
            pass
        else:
            return False
    return True
    
def no_collisions(x, state):
    if state[tuple(x)] != EMPTY:
        return False
    else:
        return True

def bonds(chain, stat):

    bonds = 0
    for j in range(chain.shape[0]):
        molecule_pos = chain[j]
        molecule = stat[tuple(molecule_pos)]
        #print j
        if molecule == REGDNA:
            continue
 
        elif molecule in BSITE_R:
            #print "BSITE_R", molecule, BSITE_R, BSITE_BINDER, BSITE_BINDER.keys(), type(BSITE_BINDER.keys()[0])
            #binding = ([BSITE_BINDER[molecule]] if not isinstance(BSITE_BINDER[molecule], list) else BSITE_BINDER[molecule])
            binding = BSITE_BINDER[molecule]
        elif molecule == BSITE_L:
            binding = [LAMIN] 
        else: print "INNY STAN!!!", j, molecule
        #print j, molecule
        
        one_ch_at = 0
        for bmove in BMOVES: 
            new = molecule_pos + bmove
            enc = stat[tuple(new)]
            #print enc, binding
            if enc in binding:
                bonds += 1
                #one_ch_at += 1
                #if one_ch_at > 6: break # one atom can have only 6 binding partners 

    return bonds

def modify(chain, binders, state, bound = BOUND):
    #move binders
    if random.randint(0, 1):
        i = random.randint(0, len(binders) - 1)
        move = random.choice(MOVES)
        new = move + binders[i]

        if no_collisions(tuple(new), state):
            return False, i, move
    #move residues
    else:
        i = random.randint(0, len(chain) - 1)
        move = random.choice(MOVES)
        new = move + chain[i]
        if good_neighbors(new, i, chain) and no_collisions(tuple(new), state):  # test if there is no collisions (the same place by different atoms) and no intersect of bonds
            if i != len(chain) - 1:
                if dist(chain[numpy.absolute(i-1)], new) <= numpy.sqrt(2) and dist(chain[numpy.absolute(i+1)], new) <= numpy.sqrt(2) and not intersect(new, chain[numpy.absolute(i-1)], state, chain) and not intersect(new, chain[numpy.absolute(i+1)], state, chain):
                    #print "Nie przecin", i, chain[i], chain[i+1], chain[i-1] 
                    return True, i, move
                else: pass
                         
            elif dist(chain[numpy.absolute(i-1)], new) <= numpy.sqrt(2) and not intersect(new, chain[numpy.absolute(i-1)], state, chain):
                #print "Nie przecin", i
                return True, i, move
                
            else:
                pass
                #print "Za duza odleglosc", i
        else:
            pass 
            #print i, "No movement"
    return None

DIST = 3
def write_as_pdb(chain, binders, attached_to_lamins, state, f, nb, metr_step, name = "chromosome and binders"):

    l = chain.shape[0]
    n = binders.shape[0]
    f.write("MODEL %i %d step %i\nTITLE %s" % (nb, l + n, metr_step, name))
    at_nr = 0

    def pdb_line(at_name, at_nr, desc, pos):
        return "\nATOM  " + str(at_nr).rjust(5) + " " + at_name.center(4) + " " + desc + "  " + str(at_nr).rjust(4) + "    " + str(round(pos[0] * DIST, 3)).rjust(8) + str(round(pos[1] * DIST, 3)).rjust(8) + str(round(pos[2] * DIST, 3)).rjust(8) + "  0.00 00.00"

    for i in range(l):
        at_n = 'C'
        cur_chain = chain[i]
        if state[tuple(cur_chain)] == REGDNA:
            r = "UNB"
        elif state[tuple(cur_chain)] in BSITE_R:
            if len(str(state[tuple(cur_chain)])) == 1:
                r = "SI"+str(state[tuple(cur_chain)])
            else: 
                r = "S"+str(state[tuple(cur_chain)])
        else: 
            r = "LAM" #BSITE_L
        #else: # BSITE_L
        #    #print type(tuple(cur_chain)), type(attached_to_lamins), tuple(cur_chain), attached_to_lamins
        #    if tuple(cur_chain) in attached_to_lamins:
        #        r = "LAS"
        #    else:
        #        r = "NLA"
        at_nr += 1
        f.write(pdb_line(at_n, at_nr, r, chain[i]))


    chain_at = at_nr

    for i in range(n):
        at_nr += 1
        ap_binder = binders[i]
        if len(str(state[tuple(ap_binder)])) == 1:
            r = "BI" + str(state[tuple(ap_binder)])
        else:
            r = "B" + str(state[tuple(ap_binder)])
        #r = "BIN"
        at_n = 'O'
        f.write(pdb_line(at_n, at_nr, r, ap_binder))

    #for hl in highlighted_lamins:
    #    at_nr += 1
    #    r = "HLA" 
    #    at_n = 'P'
    #    f.write(pdb_line(at_n, at_nr, r, hl))

    for i in range(1, chain_at - 1):
        if i == 1:
            line = "\nCONECT" + str(i).rjust(5) +  str(i + 1).rjust(5)
        else: line = "\nCONECT" + str(i).rjust(5) + str(i - 1).rjust(5) + str(i + 1).rjust(5)
        f.write(line)
    f.write("\nENDMDL   \n")

def count_bonds(pos, accepted, state):
    bonds = 0
    for bmove in BMOVES:
        if state[tuple(pos + bmove)] in accepted:
             bonds += 1
      #       if bonds > 6: break
    #print bonds
    return bonds

def radius_gyr(chai):
    length = chai.shape[0]
    r = 0.0
    for at in range(length):
        for ato in range(at):
            r = r + math.sqrt(numpy.sum((chai[at] - chai[ato])**2))
    r_gyr = 2*r/(length*(length-1))
    return r_gyr


DELTA=2
GYRATION = True
CHECK_E = False

def metropolis(chain, binders, attached_to_lamins, state, out_fname, name = "chromosome", n = 100):

    if GYRATION:
        out_file_out.write("iter, step, Energy, R_gyr\n")
    else:
        out_file_out.write("iter, step, Energy\n")

    def put_as_pickle(p_out,  p_chain, p_binders, p_attached_to_lamins, p_state, p_bindNR, p_bsites, p_bindersState, p_bsiteBinder):
        # dump the last state to the pickle
        l_obj = [p_chain, p_binders, p_attached_to_lamins, p_state, p_bindNR, p_bsites, p_bindersState, p_bsiteBinder]
        pickle_fname = p_out.split('.pdb')[0] + ".pick"
        pickle_file = open(pickle_fname, 'w')
        pickle.dump(l_obj, pickle_file)
    
    def put_as_json(p_out,  p_chain, p_binders, p_attached_to_lamins, p_state, p_bindNR, p_bsites, p_bindersState, p_bsiteBinder):
        p_chain = p_chain.tolist()
        p_binders = p_binders.tolist()
        p_state = p_state.tolist()
        l_obj = [p_chain, p_binders, p_attached_to_lamins, p_state, p_bindNR, p_bsites, p_bindersState, p_bsiteBinder]
        #print type(p_chain), type(p_binders), type(p_attached_to_lamins), type(p_state), type(p_bindNR), type(p_bsites), type(p_bindersState)
        json_fname = p_out.split('.pdb')[0] + ".json"
        json_file = open(json_fname, 'w')
        json.dump(l_obj, json_file)

    out_file = open(out_fname, "w")
    st_nr = 0
    pick_step = opts.Save
    E = bonds(chain, state)
    out_file_out.write("Starting energy: %f\n" %E)
    write_as_pdb(chain, binders, attached_to_lamins, state, out_file, st_nr, 0, name + ";bonds=" + str(E))

    for step in xrange(n):

        resp = modify(chain, binders, state)

        ch = numpy.array(chain, copy = True)
        b = numpy.array(binders, copy = True)
        if CHECK_E:
            st = numpy.array(state,copy=True)
        else: pass
        if resp:
            in_chain, i, move = resp
            if in_chain:
                #print i, "change"
                old = numpy.copy(ch[i])
                ch[i] = ch[i] + move
                if CHECK_E:
                    st[tuple(ch[i])] = st[tuple(old)]
                    st[tuple(old)] = EMPTY
                else: pass
                #print state[tuple(old)]
                if state[tuple(old)] in BSITE_R:
                    #ac = ([BSITE_BINDER[state[tuple(old)]]] if not isinstance(BSITE_BINDER[state[tuple(old)]], list) else BSITE_BINDER[state[tuple(old)]])
                    #print ([BSITE_BINDER[state[tuple(old)]]] if not isinstance(BSITE_BINDER[state[tuple(old)]], list) else BSITE_BINDER[state[tuple(old)]])
                    Enew = E + count_bonds(ch[i], BSITE_BINDER[state[tuple(old)]], state) - count_bonds(old, BSITE_BINDER[state[tuple(old)]], state)
                    if CHECK_E:
                        Eslow = bonds(ch, st)
                        if Enew != Eslow:
                            print "R", 'Enew', Enew, 'Eslow', Eslow
                    else: pass
                elif state[tuple(old)] == BSITE_L:
                    Enew = E + count_bonds(ch[i], [LAMIN], state) - count_bonds(old, [LAMIN], state)
                    if CHECK_E:
                        Eslow = bonds(ch, st)
                        if Enew != Eslow: 
                            print "L", 'Enew', Enew, 'Eslow', Eslow
                    else: pass
                    if tuple(ch[i]) not in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) > 0:
                        #print "NOT i > 0", tuple(ch[i]),  attached_to_lamins, count_bonds(ch[i], [LAMIN], state), ch[i]
                        attached_to_lamins.append(tuple(ch[i]))
                    elif tuple(ch[i]) in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) == 0:
                        #print "IN i ==0", tuple(ch[i]),  attached_to_lamins, count_bonds(ch[i], [LAMIN], state), ch[i]
                        attached_to_lamins.remove(ch[i])
                    #elif tuple(ch[i]) in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) > 0:
                    #    print "JEST, ale ma lamine!", ch[i], count_bonds(ch[i], [LAMIN], state)
                    else:
                        pass 
                        #print "WYJEATEK ", tuple(ch[i]), attached_to_lamins, count_bonds(ch[i], [LAMIN], state)

                else: # REGDNA
                    Enew = E
            else:
                old = numpy.copy(b[i])
                b[i] = b[i] + move
                if CHECK_E:
                    st[tuple(b[i])] = st[tuple(old)]
                    st[tuple(old)] = EMPTY
                else: pass
                bsite_r_accep = [ke for ke, va in BSITE_BINDER.items() if state[tuple(old)] in va] # find keys based on values of dictionary
                #print b[i],old, tuple(b[i]), state[tuple(old)], BINDER.index(state[tuple(old)]), [BSITE_R[BINDER.index(state[tuple(old)])]]
                #Enew = E + count_bonds(b[i],  [BSITE_R[BINDER.index(state[tuple(old)])]], state) - count_bonds(old,  [BSITE_R[BINDER.index(state[tuple(old)])]], state) 
                Enew = E + count_bonds(b[i],  bsite_r_accep, state) - count_bonds(old,  bsite_r_accep, state) 
                if CHECK_E:
                    Eslow = bonds(ch, st)
                    if Enew != Eslow:
                        print  "B, Enew", Enew, 'Eslow', Eslow
                else: pass
        else:
            Enew = E

        if Enew >= E or random.uniform(0.0, 1.0) < math.exp((Enew - E) * DELTA): #accept

            chain = ch
            binders = b
            #write_as_pdb(chain, binders, attached_to_lamins, state, out_file, st_nr, name + ";bonds=" + str(E))
            if resp:
                state[tuple(old + move)] = state[tuple(old)]
                state[tuple(old)] = EMPTY

            if E != Enew:

                E = Enew
            st_nr += 1
            if (st_nr%opts.Save)==0 or st_nr == opts.Steps: ###ZAPISUJE KAZDE 100 KROKOW!!!!
                if GYRATION:
                    out_file_out.write ("%i, %i, %i, %f\n" %(step, st_nr, E, radius_gyr(chain)))
                else:
                    out_file_out.write ("%i, %i, %f\n" %(step, st_nr, E))
                #print st_nr
                write_as_pdb(chain, binders, attached_to_lamins, state, out_file, st_nr, step, name + ";bonds=" + str(E))
                #print "WRITE!!!"
            if st_nr == pick_step or st_nr == opts.Steps:
                pick_step += opts.Save
                #put_as_pickle(out_fname, chain, binders, attached_to_lamins, state, M,  BSITE_R, BINDER, BSITE_BINDER)
                put_as_json(out_fname, chain, binders, attached_to_lamins, state, M,  BSITE_R, BINDER, BSITE_BINDER)

    #put_as_pickle(out_fname, chain, binders, attached_to_lamins, state, M,  BSITE_R, BINDER)
    put_as_json(out_fname, chain, binders, attached_to_lamins, state, M,  BSITE_R, BINDER, BSITE_BINDER)
    out_file.close()


def output_name(ou, m, n):
    if ou == '':
        f_n = "MC_traj_%ibin_%ichain.pdb" % (m[0], n)
    elif "." in ou:
        f_n = ou.split('.')[0] + ".pdb"
    else:
        f_n = ou + ".pdb"
    return f_n



opts = pars_inp()


if opts.In_str == '':
    rand_init = True
else: 
    rand_init = False


t1 = time.time()
if rand_init: 
    N = opts.Ch_lenght # length of the chain
    M = [int(e) for e in opts.Binders.split(",")] # nr of binders
    fn = output_name(opts.Out_str, M, N)
    BSITE_pre = []
    for r in range(len(opts.Regular_bsites.split(","))):
        BSITE_pre.append(4+r)
        #BINDER_pre.append(4+len(opts.Regular_bsites.split(","))+r)
    
    import itertools # add all possible combinations of binders
    BSITE_R = list(BSITE_pre)
    for l in range(2, len(BSITE_R)+1):
        for subset in itertools.combinations(BSITE_R, l):
            BSITE_pre.append(list(subset))
    BSITE_R = range(4,4+len(BSITE_pre))
    BSITE_dic = dict(zip(BSITE_R, BSITE_pre)) # dictionary for more than one type of connected binder
    BSITE_dic_re = dict(zip(map(lambda x: tuple(x) if isinstance(x, list) else x, BSITE_pre), BSITE_R))
    
    for r in range(len(opts.Regular_bsites.split(","))):
        BINDER.append(4+len(BSITE_R)+r)
    
    BINDER_pre = list(BINDER)
    BINDER_pre = list(map(lambda el:[el], BINDER_pre))    
    for l in range(2, len(BINDER)+1):
        for subset in itertools.combinations(BINDER, l):
            BINDER_pre.append(list(subset))
    
    BSITE_BINDER = dict(zip(BSITE_R, BINDER_pre))
            
    print BSITE_R, BINDER, BSITE_dic, BSITE_dic_re, BINDER_pre, BSITE_BINDER
    
    c, b, a, state = initialize_random(N, M, fn)
else: 
    c, b, a, state, M, BSITE_R,  BINDER, BSITE_BINDER = initialize_import_json(opts.In_str)
    N = c.shape[0]
    #M = b.shape[0]
    fn = output_name(opts.Out_str, M, N)


outputn_out = fn.split('.pdb')[0] + ".out"
out_file_out = open(outputn_out, 'w')

out_file_out.write( "The lenght of the chain is %i,  the number of binders is  %s,  the nucleus radius is %i, the number of steps is %i, bsites %s, lamin %s\n" %(N, M, R, opts.Steps, opts.Regular_bsites, opts.Lamin_bsites))

t2 = time.time()
out_file_out.write( "initialization: %f\n" %(t2 - t1))
BOUND = numpy.max(c)

a = []
metropolis(c, b, a, state, fn, n = opts.Steps)
t1 = t2
t2 = time.time()
out_file_out.write( "metropolis: %f\n" %(t2 - t1))
