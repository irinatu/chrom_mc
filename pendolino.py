import numpy,math,time, optparse, sys, pickle
import random as random

#accepted move vectors
MOVES = numpy.array([[0,0,1], [0,1,0], [1,0,0], [0,0,-1], [0,-1,0], [-1,0,0], [1,0,1], [0,1,1], [1,1,0], [-1,0,-1], [0,-1,-1], [-1,-1,0], [1,0,-1], [0,1,-1], [1,-1,0], [-1,0,1], [0,-1,1], [-1,1,0], [-1,1,-1], [-1,-1,-1], [1,1,1], [1,-1,1], [1,-1,-1], [-1,-1,1], [-1,1,1], [1,1,-1]])
#accepted matching positions of binding sites
BMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]])

# radius of the nucleus
R = 100
# 2 x radius + a fringe, because lamin barrier has to be hermetic
BOUND = 2 * R + 2
random.seed(1)

EMPTY = 0
BINDER = 1
LAMIN = 2
BSITE_R = 3
BSITE_L = 4
REGDNA = 5

DIST = 3

# max distance of good neighbors
GOOD_NEIGH = 3

def pars_inp():
    ## read/validate command-line arguments
    optparser = optparse.OptionParser(usage = "%prog regular_dsites.txt lamin_bsites.txt [<options>]")

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
    optparser.add_option('-b', type = "int",
        dest = "Binders",
        default = 256,
        help = "Number of binders (default 256)")
			
    (opts, args) = optparser.parse_args()

    if len(args) < 2:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)
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
    return ch, b, a, state

def getStateWithLamins(bound, f):

    state = numpy.zeros((bound, bound, bound), dtype = numpy.int)
    MIDDLE = bound / 2
    lam_name = f.split('.pdb')[0] + '_lamin.pdb'
    save_lam = open(lam_name, "w")
    save_lam.write("HEADER LAMINA")
    at_nr = 1

    def dist(x, y, z):
        return math.sqrt((x - MIDDLE)**2 + (y - MIDDLE)**2 + (z - MIDDLE)**2)

    for x in range(BOUND):
        for y in range(BOUND):
            for z in range(BOUND):
                border = abs(dist(x, y, z) - MIDDLE + 1)
                if border <= 2:
                    state[x, y, z] = LAMIN
                    if border == 1:
                        #print border
                        line = "\nATOM  " + str(at_nr).rjust(5) + " " + "O".center(4) + " " + "LAM" + "  " + str(at_nr).rjust(4) + "    " + str(round(x * DIST, 3)).rjust(8) + str(round(y * DIST, 3)).rjust(8) + str(round(z * DIST, 3)).rjust(8) + "  0.00 00.00"
                        at_nr += 1
                        save_lam.write(line)
                    
    save_lam.close()
    return state

def initialize_random(n, m, fa, bound = BOUND):

    chain = numpy.zeros((n, 3), dtype = numpy.int)
    binders = numpy.zeros((m, 3), dtype = numpy.int)
    state = getStateWithLamins(bound, fa)
    attached_to_lamins = []

    chain[0] = [bound / 2] * 3

    def get_site_type_list(fpath, length):
        positions = [0] * length
        for l in open(fpath):
            positions[int(l)] = 1
        return positions

    regular_bsites = get_site_type_list(sys.argv[1], n)
    lamin_bsites   = get_site_type_list(sys.argv[2], n)

    def get_site_type(i, regular_bsites, lamin_bsites): # BSITE_R interacts with binders whereas BSITE_L interacts both with lamins and binders
        if regular_bsites[i] == 1:
            return BSITE_R
        elif lamin_bsites[i] == 1:
            return BSITE_L
        return REGDNA

    cur = chain[0]
    state[tuple(cur)] = get_site_type(0, regular_bsites, lamin_bsites)
    for i in range(1, n):
        mov = random.choice(MOVES)
        tries = 0
        while tries < 100 and not (no_collisions(tuple(cur + mov), state)):
            mov = random.choice(MOVES)
            tries += 1
        assert tries != 100, "unable to find initialization"
        chain[i] = cur + mov
        state[tuple(chain[i])] = get_site_type(i, regular_bsites, lamin_bsites)

        if state[tuple(chain[i])] == BSITE_L and count_bonds(chain[i], [LAMIN], state) > 0:
            attached_to_lamins.append(tuple(chain[i]))

        cur = chain[i]

    for i in range(m):
        cur = chain[i * n / m]
        mov = 2 * random.choice(MOVES)
        tries = 0
        while tries < 100 and not (no_collisions(tuple(cur + mov), state)):
            mov = 2 * random.choice(MOVES)
            tries += 1
        binders[i] = cur + mov
        state[tuple(binders[i])] = BINDER

    return chain, binders, attached_to_lamins, state

def dist(r1,r2):
    return DIST_MATRIX[tuple(abs(r1 - r2))]

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
    return True

def bonds(chain, state):

    bonds = 0
    for j in range(chain.shape[0]):
        molecule_pos = chain[j]
        molecule = state[tuple(molecule_pos)]

        if molecule == REGDNA:
            continue
 
        if molecule == BSITE_R:
            binding = [BINDER]
        elif molecule == BSITE_L:
            binding = [LAMIN, BINDER]

        for bmove in BMOVES: 
            new = molecule_pos + bmove
            enc = state[tuple(new)]
            if enc in binding:
                bonds += 1 

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
        if good_neighbors(new, i, chain) and no_collisions(tuple(new), state):
            return True, i, move
    return None

DIST = 3
def write_as_pdb(chain, binders, attached_to_lamins, state, f, name = "chromosome and binders"):

    l = chain.shape[0]
    n = binders.shape[0]
    f.write("HEADER %d\nTITLE %s" % (l + n, name))
    at_nr = 0

    def pdb_line(at_nr, desc, pos):
        return "\nATOM  " + str(at_nr).rjust(5) + " " + "C".center(4) + " " + desc + "  " + str(at_nr).rjust(4) + "    " + str(round(pos[0] * DIST, 3)).rjust(8) + str(round(pos[1] * DIST, 3)).rjust(8) + str(round(pos[2] * DIST, 3)).rjust(8) + "  0.00 00.00"

    for i in range(l):
        cur_chain = chain[i]
        if state[tuple(cur_chain)] == REGDNA:
            r = "UNB"
        elif state[tuple(cur_chain)] == BSITE_R:
            r = "BOU"
        else: # BSITE_L
            if tuple(cur_chain) in attached_to_lamins:
                r = "LAS"
            else:
                r = "NLA"
        at_nr += 1
        f.write(pdb_line(at_nr, r, chain[i]))

    def neighborhood(pos, type_to_search, size = 2):
        res = []
        for i in range(-size, size + 1):
            for j in range(-size, size + 1):
                for k in range(-size, size + 1):
                    p = pos + numpy.array([i, j, k])
                    if state[tuple(p)] == type_to_search:
                        res.append(p)
        return res

    highlighted_lamins = []
    for pos in attached_to_lamins:
        # find lamins nearby
        neigh = neighborhood(pos, LAMIN)
        # if they are not in the list, add them
        for nei in neigh:
            if not tuple(nei) in highlighted_lamins:
                highlighted_lamins.append(tuple(nei))

    chain_at = at_nr

    for i in range(n):
        at_nr += 1
        r = "BIN"
        f.write(pdb_line(at_nr, r, binders[i]))

    for hl in highlighted_lamins:
        at_nr += 1
        r = "HLA" 
        f.write(pdb_line(at_nr, r, hl))

    for i in range(1, chain_at - 1):
        line = "\nCONECT" + str(i).rjust(5) +  str(i + 1).rjust(5)
        f.write(line)
    f.write("\nEND   \n")

def count_bonds(pos, accepted, state):
    bonds = 0
    for bmove in BMOVES:
        if state[tuple(pos + bmove)] in accepted:
             bonds += 1
    return bonds

DELTA=2
def metropolis(chain, binders, attached_to_lamins, state, fn, name = "chromosome", n = 100):

    E = bonds(chain, state)

    for step in range(n):

        resp = modify(chain, binders, state)

        ch = numpy.array(chain, copy = True)
        b = numpy.array(binders, copy = True)
        if resp:
            in_chain, i, move = resp
            if in_chain:
                old = numpy.copy(ch[i])
                ch[i] = ch[i] + move
                if state[tuple(old)] == BSITE_R:
                    Enew = E + count_bonds(ch[i], [BINDER], state) - count_bonds(old, [BINDER], state)
                elif state[tuple(old)] == BSITE_L:
                    Enew = E + count_bonds(ch[i], [LAMIN, BINDER], state) - count_bonds(old, [LAMIN, BINDER], state)
                    if tuple(ch[i]) not in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) > 0:
                        attached_to_lamins.append(ch[i])
                    if tuple(ch[i]) in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) == 0:
                        attached_to_lamins.remove(ch[i])

                else: # REGDNA
                    Enew = E
            else:
                old = numpy.copy(b[i])
                b[i] = b[i] + move
                Enew = E + count_bonds(b[i], [BSITE_R], state) - count_bonds(old, [BSITE_R], state)
        else:
            Enew = E

        if Enew > E or random.uniform(0.0, 1.0) < math.exp((Enew - E) * DELTA): #accept

            if E != Enew:

                traj.append((ch, b, Enew))

                E = Enew
                f=open(fn, "a")
                print "iter", step, "energy:", E, "accepted:", len(traj)
                write_as_pdb(chain, binders, attached_to_lamins, state, f, name + ";frame=" + str(len(traj)) + ";bonds=" + str(E))
                f.close()

            chain = ch
            binders = b
            if resp:
                state[tuple(old + move)] = state[tuple(old)]
                state[tuple(old)] = EMPTY

    # dump the last state to the pickle
    l_obj = [ch, b, attached_to_lamins, state]
    fn = fn.split('.pdb')[0] + ".pick"
    file = open(fn, 'w')
    pickle.dump(l_obj, file)
    return traj

opts = pars_inp()

if opts.In_str == '':
    rand_init = True
else: 
    rand_init = False

N = opts.Ch_lenght # length of the chain
M = opts.Binders # nr of binders

if opts.Out_str == '':
    fn = "MC_traj_%ibin_%ichain.pdb" % (M, N)
elif "." in opts.Out_str:
    fn = opts.Out_str.split('.')[0] + ".pdb"
else:
    fn = opts.Out_str + ".pdb"

t1 = time.time()
if rand_init: 
    c, b, a, state = initialize_random(N, M, fn)
else: 
    c, b, a, state = initialize_import(opts.In_str)
t2 = time.time()
print "initialization: ", t2 - t1
BOUND = numpy.max(c)

t = metropolis(c, b, a, state, fn, n = opts.Steps)
t1 = t2
t2 = time.time()
print "metropolis: ", t2 - t1
