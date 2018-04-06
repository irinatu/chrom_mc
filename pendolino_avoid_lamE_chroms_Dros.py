import numpy,math,time, optparse, sys, pickle
import random as random

#accepted move vectors
#MOVES = numpy.array([[0,0,1], [0,1,0], [1,0,0], [0,0,-1], [0,-1,0], [-1,0,0], [1,0,1], [0,1,1], [1,1,0], [-1,0,-1], [0,-1,-1], [-1,-1,0], [1,0,-1], [0,1,-1], [1,-1,0], [-1,0,1], [0,-1,1], [-1,1,0], [-1,1,-1], [-1,-1,-1], [1,1,1], [1,-1,1], [1,-1,-1], [-1,-1,1], [-1,1,1], [1,1,-1]])
MOVES = numpy.array([[0,0,1], [0,1,0], [1,0,0], [0,0,-1], [0,-1,0], [-1,0,0], [1,0,1], [0,1,1], [1,1,0], [-1,0,-1], [0,-1,-1], [-1,-1,0], [1,0,-1], [0,1,-1], [1,-1,0], [-1,0,1], [0,-1,1], [-1,1,0]])
#accepted matching positions of binding sites
#BMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]])
BMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1], [1,1,0], [-1,-1,0], [-1,1,0], [1,-1,0], [1,0,1], [-1,0,1], [1,0,-1], [-1, 0, -1], [0,1,-1], [0,-1,-1], [0,-1,1], [0,1,1]])
#TERMOVES = numpy.array([[-2, -2, -2], [-2, -2, -1], [-2, -2, 0], [-2, -2, 1], [-2, -2, 2], [-2, -1, -2], [-2, -1, -1], [-2, -1, 0], [-2, -1, 1], [-2, -1, 2], [-2, 0, -2], [-2, 0, -1], [-2, 0, 0], [-2, 0, 1], [-2, 0, 2], [-2, 1, -2], [-2, 1, -1], [-2, 1, 0], [-2, 1, 1], [-2, 1, 2], [-2, 2, -2], [-2, 2, -1], [-2, 2, 0], [-2, 2, 1], [-2, 2, 2], [-1, -2, -2], [-1, -2, -1], [-1, -2, 0], [-1, -2, 1], [-1, -2, 2], [-1, -1, -2], [-1, -1, -1], [-1, -1, 0], [-1, -1, 1], [-1, -1, 2], [-1, 0, -2], [-1, 0, -1], [-1, 0, 0], [-1, 0, 1], [-1, 0, 2], [-1, 1, -2], [-1, 1, -1], [-1, 1, 0], [-1, 1, 1], [-1, 1, 2], [-1, 2, -2], [-1, 2, -1], [-1, 2, 0], [-1, 2, 1], [-1, 2, 2], [0, -2, -2], [0, -2, -1], [0, -2, 0], [0, -2, 1], [0, -2, 2], [0, -1, -2], [0, -1, -1], [0, -1, 0], [0, -1, 1], [0, -1, 2], [0, 0, -2], [0, 0, -1], [0, 0, 1], [0, 0, 2], [0, 1, -2], [0, 1, -1], [0, 1, 0], [0, 1, 1], [0, 1, 2], [0, 2, -2], [0, 2, -1], [0, 2, 0], [0, 2, 1], [0, 2, 2], [1, -2, -2], [1, -2, -1], [1, -2, 0], [1, -2, 1], [1, -2, 2], [1, -1, -2], [1, -1, -1], [1, -1, 0], [1, -1, 1], [1, -1, 2], [1, 0, -2], [1, 0, -1], [1, 0, 0], [1, 0, 1], [1, 0, 2], [1, 1, -2], [1, 1, -1], [1, 1, 0], [1, 1, 1], [1, 1, 2], [1, 2, -2], [1, 2, -1], [1, 2, 0], [1, 2, 1], [1, 2, 2], [2, -2, -2], [2, -2, -1], [2, -2, 0], [2, -2, 1], [2, -2, 2], [2, -1, -2], [2, -1, -1], [2, -1, 0], [2, -1, 1], [2, -1, 2], [2, 0, -2], [2, 0, -1], [2, 0, 0], [2, 0, 1], [2, 0, 2], [2, 1, -2], [2, 1, -1], [2, 1, 0], [2, 1, 1], [2, 1, 2], [2, 2, -2], [2, 2, -1], [2, 2, 0], [2, 2, 1], [2, 2, 2]])

TERMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1], [1,1,0], [-1,-1,0], [-1,1,0], [1,-1,0], [1,0,1], [-1,0,1], [1,0,-1], [-1, 0, -1], [0,1,-1], [0,-1,-1], [0,-1,1], [0,1,1]])

# radius of the nucleus
R = 20
# 2 x radius + a fringe, because lamin barrier has to be hermetic
BOUND = 2 * R + 2
#ran_seed = 5
#random.seed(ran_seed)
#print ran_seed

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
    optparser.add_option('-l', type = "string",
        dest = "Ch_lenght",
        default = 512,
        help = "Lenght of the chains, separated by comma (two chains: 234, 456) (default one chain with 512 bins lenght)")
    optparser.add_option('-r', type = "string",
        dest = "Revers",
        default = 1058,
        help = "Lenght of the chains, where tail should be near the center (2R, 3R), separated by comma (two chains: 234, 456) (default one chain with 512 bins lenght)")
    optparser.add_option('-b', type = "int",
        dest = "Binders",
        default = 256,
        help = "Number of binders (default 256)")
			
    (opts, args) = optparser.parse_args()

    if len(args) < 2 and opts.In_str =='':
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
    print len(list_ob)
    ch = list_ob[0]
    b = list_ob[1]
    a = list_ob[2]
    state = list_ob[3]
    ch_n = list_ob[4]
    re = list_ob[5]
    return ch, b, a, state, ch_n, re

def dist_from_mi(x, y, z, mi):
        return math.sqrt((x - mi)**2 + (y - mi)**2 + (z - mi)**2)

def getStateWithLamins(bound, f):

    state = numpy.zeros((bound, bound, bound), dtype = numpy.int)
    MIDDLE = bound / 2
    lam_name = f.split('.pdb')[0] + '_lamin.pdb'
    save_lam = open(lam_name, "w")
    save_lam.write("HEADER LAMINA")
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
                    
    save_lam.close()
    return state

def getTerritories(angles):
    # angles shoud be in the format [(0,20),(20,70),(70,200), (200,360)]
    stateT = numpy.zeros((BOUND, BOUND, BOUND), dtype = numpy.int)
    
    radius_cut = round(MIDDLE/5.0)
    start = 0.
    radius_cut_list = []
    for ra in range(5):
        start = start+radius_cut
        radius_cut_list.append(start)
        
    
    for katy in angles: 
        kat0 =  katy[0]
        kat1 = katy[1]

        for x in range(BOUND):
            for y in range(BOUND):
                for z in range(BOUND):
                    border = dist_from_mi(x-MIDDLE, y-MIDDLE, z-MIDDLE, 0)
                    if border <= MIDDLE and border!= 0:
                        kat_tet = math.atan2((y-MIDDLE),(x-MIDDLE))
                        if  (kat0 < 180 and kat1 < 180 and math.radians(kat0) <= kat_tet < math.radians(kat1)) or (kat0 < 180 and kat1 > 180 and ((math.radians(kat0) <= kat_tet <= math.pi) or (-math.pi < kat_tet < (-math.pi+math.radians(kat1-180))))) or (kat0 > 180 and kat1 > 180 and (-math.pi+math.radians(kat0-180) < kat_tet < -math.pi+math.radians(kat1-180))):
                            bor = dist_from_mi(x-MIDDLE, y-MIDDLE, MIDDLE-MIDDLE, 0)
                            if 0. <= bor < radius_cut_list[0] : stateT[x, y, z] = angles.index(katy)+1
                            elif radius_cut_list[0] <= bor < radius_cut_list[1]: stateT[x, y, z] = str(angles.index(katy)+1) + '1'
                            elif radius_cut_list[1] <= bor < radius_cut_list[2]: stateT[x, y, z] = str(angles.index(katy)+1) + '2'
                            elif radius_cut_list[2] <= bor < radius_cut_list[3]: stateT[x, y, z] = str(angles.index(katy)+1) + '3'
                            elif radius_cut_list[3] <= bor < MIDDLE: stateT[x, y, z] = str(angles.index(katy)+1) + '4'

    return stateT


def dist(r1,r2):
    return DIST_MATRIX[tuple(abs(r1 - r2))]
    
#def distance (p1, p2):
 #   return numpy.sqrt(numpy.sum((p1-p2)**2))

def cross(cha, po1, po2):
    pier = numpy.where((cha == po1).all(axis=1))
    dru = numpy.where((cha == po2).all(axis=1))
    #print "pir-dru ", po1, po2, len(pier), len(dru), pier, dru, pier[0][0]
    if len(pier) == 1 and len(dru) == 1:
        if pier[0][0] == dru[0][0]+1 or pier[0][0] == dru[0][0]-1:
            #print 'TRUE cross'
            return True
        else:
            #print 'FALSE cross'
            return False
    else: print "Two atoms of the chain have the same coordinates. Please check ", po1, po2, pier[0], dru[0]
    
def intersect(new_p, next, sta, ch): #coordinates of two next polimer points
    #print dist(next, new_p), distance(next, new_p)
    #print "INTER"
    if  dist(next, new_p) == 1: 
        return False
    elif dist(next, new_p) > 1:
        differ = new_p - next
        if differ[0] == 0:
            pos1 = numpy.array(([next[0], next[1], new_p[2]]))
            pos2 = numpy.array(([next[0], new_p[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN, BINDER] and sta[tuple(pos2)] not in [EMPTY, LAMIN, BINDER]:
                return cross(ch, pos1, pos2)
            else: return False         
        elif differ[1] == 0:
            pos1 = numpy.array(([next[0], next[1], new_p[2]]))
            pos2 = numpy.array(([new_p[0], next[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN, BINDER] and sta[tuple(pos2)] not in [EMPTY, LAMIN, BINDER]:
                return cross(ch, pos1, pos2)
            else: return False
        elif differ[2] == 0:
            pos1 = numpy.array(([next[0], new_p[1], next[2]]))
            pos2 = numpy.array(([new_p[0], next[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN, BINDER] and sta[tuple(pos2)] not in [EMPTY, LAMIN, BINDER]:
                return cross(ch, pos1, pos2)
            else: return False
        else: print new_p, next, "The distance between these two positions are not = sqrt(2)"

def no_collisions(x, state):
    #print "COLI", x, state.shape, state[tuple(x)]
    if state[tuple(x)] != EMPTY:
        return False
    else:
        return True
        
def get_angels(at_nr):
    stos = [round(float(a)/min(at_nr)) for a in at_nr]
    ang_st = 360./sum(stos)
    ang = [round(s * ang_st) for s in stos]
    mi = min(ang)
    if mi < 15:
        rest = 15-mi
        ma = max(ang)
        ang[ang.index(mi)] = 15
        ang[ang.index(ma)] = ang[ang.index(ma)] - rest
    an0 = ang[0]
    angels = [(0, ang[0])]
    for an in ang[1:-1]:
         an1 = an0+an
         angels.append((an0, an1))
         an0 = an1
    angels.append((an0, 360))
    return angels
    
        
    

MIDDLE = BOUND / 2
def initialize_random(n, m, fa, rever, bound = BOUND): # n - list with chains atoms numbers, m - number of binders, fa - file name of lamins
    def fibonacci_sphere(samples=1,randomize=True): # distribute points on a sphere from web
        rnd = 1.
        if randomize:
            rnd = random.random() * samples
        points = []
        offset = 2./samples
        increment = math.pi * (3. - math.sqrt(5.));
        for i in range(samples):
            y = round(((i * offset) - 1) + (offset / 2))
            r = math.sqrt(9 - pow(y,2)) # R ~ 3
            #print "R", r
            phi = ((i + rnd) % samples) * increment
            y = y + MIDDLE
            x = round(math.cos(phi) * r) + MIDDLE
            z = round(math.sin(phi) * r) + MIDDLE
            points.append([x,y,z])
            #print points
        return points

    def get_site_type_list(fpath, length_list):
        positions = []
        chrom = -1
        #print "lent", length_list
        for length in length_list:
            pos = [0] * length
            positions.append(pos)
        for l in open(fpath):
                if "chr" in l:
                    #print "tak", l
                    chrom +=1
                    continue
                #print l, chrom
                positions[chrom][int(l) -1] = 1
        return positions
        
    def get_site_type(i, regular_bsites, lamin_bsites): # BSITE_R interacts with binders whereas BSITE_L interacts both with lamins and binders
        if lamin_bsites[i] == 1:
            return BSITE_L
        elif regular_bsites[i] == 1:
            return BSITE_R
        else:
            return REGDNA
    
    def check_gyr(atoms, nr_a, chai):
        r = 0.0
        for at in range(atoms-nr_a, atoms+1):
            for ato in range(at, atoms+1):
                r = r + math.sqrt(numpy.sum((chai[at] - chai[ato])**2))
                #print at, ato, r, chai[at], chai[ato], (chai[at] - chai[ato])**2,numpy.sum((chai[at] - chai[ato])**2) 
        ch_gyr = 2*r/float(nr_a)
        #print "NUMERY", atoms-nr_a, atoms+1, nr_a, ch_gyr
        return ch_gyr
        
    def rand_next(cu, st, ch, terri, nun, ti, at =1, ii=1):
        mov = random.choice(MOVES)
        tries = 0
        #print "dfgf"
        ch_cop = numpy.copy(ch)
        ch_cop[at]=cu+mov
        #while (tries < 100) and (not (no_collisions(tuple(cu + mov), st)) or intersect(cu, cu+mov, st, ch)  or terri[tuple(cu+mov)] != 0 or check_gyr(at, ii, ch_cop)>nun/5.):
        while (tries < 100) and (not (no_collisions(tuple(cu + mov), st)) or intersect(cu, cu+mov, st, ch)  or terri[tuple(cu+mov)] != ti):
        #while tries < 100 and (not (no_collisions(tuple(cu + mov), st)) or intersect(cu, cu+mov, st, ch)):
            #print "TRTRT"
            mov = random.choice(MOVES)
            tries += 1
            try:
                terri[tuple(cu+mov)]
            except IndexError: 
                print "Error"
                mov = random.choice(MOVES)
            ch_cop[at]=cu+mov
        if tries == 100: 
            return mov, True
        return mov, False
        
    def fill_one_terr(one_t, pos, t):
        for tm in TERMOVES:
            try:
                if one_t[tuple(tm+pos)] == t:
                    one_t[tuple(tm+pos)] = 0
            except IndexError: pass
        return one_t
    
    chain = numpy.zeros((sum(n), 3), dtype = numpy.int)
    binders = numpy.zeros((m, 3), dtype = numpy.int)
    state = getStateWithLamins(bound, fa)
    ter_ang = get_angels(n)
    territor = getTerritories(ter_ang) # state with territories for each chain defined as number form 1 to len(n)
    
    attached_to_lamins = []
    for x in range(bound):
        for y in range(bound):
            for z in range(bound):
                dist_m = dist_from_mi(x, y, z, MIDDLE)
                if dist_m <= 3:
                    territor[x, y, z] = 0   
    

    regular_bsites = get_site_type_list(sys.argv[1], n)
    lamin_bsites   = get_site_type_list(sys.argv[2], n)
    #print regular_bsites  
    
    
    at_nr = -1
    
    #points_on_sph = fibonacci_sphere(60)
    #print get_site_type(0, regular_bsites[0], lamin_bsites[0])
    for nu, re, la in zip(n, regular_bsites, lamin_bsites):
         
        cur0 = [bound / 2] * 3
        te = n.index(nu)+1
        if cur0 == [MIDDLE,MIDDLE,MIDDLE]:
            ter_point = numpy.where(territor == te)
            po_dis_min = 20.0 # big enough to find smaller dis - closest to the middle
            for px, py, pz in zip(ter_point[0], ter_point[1], ter_point[2]):
                po_dis = dist_from_mi(px, py, pz, MIDDLE)
                if po_dis_min > po_dis:
                    po_dis_min = po_dis
                    cur0 = [px,py,pz]
        
        at_nr += 1
        chain[at_nr] = numpy.array(cur0)
        print chain[0], chain[at_nr]
        if nu not in rever:
            state[tuple(chain[at_nr])] = get_site_type(0, re, la)
        else:
            state[tuple(chain[at_nr])] = get_site_type(-1, re, la)
        cur = chain[at_nr]
        print at_nr, chain[at_nr], nu
        
        ter_zero = numpy.where(territor == te)
        ter_one = numpy.where(territor == int(str(te) + '1'))
        ter_two = numpy.where(territor == int(str(te) + '2'))
        ter_three = numpy.where(territor == int(str(te) + '3'))
        ter_four = numpy.where(territor == int(str(te) + '4'))
        state_lam = numpy.where(state == LAMIN)
        try_chain = 0 
        succ = False  
        
        if nu not in rever:
            while try_chain < 100 and not succ:
                print try_chain, succ
                for i in range(1, nu):
                    chain_l_dev = round(nu/125.)
                    at_nr += 1
                    if chain_l_dev <= i < chain_l_dev * 10. and len(numpy.where(territor == int(str(te) + '1'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '1'
                        open_part = numpy.where(territor == int(str(te) + '1'))
                        territor[open_part]= te
                    elif chain_l_dev * 10. <= i < chain_l_dev * 30. and len(numpy.where(territor == int(str(te) + '2'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '2'
                        open_part = numpy.where(territor == int(str(te) + '2'))
                        territor[open_part]= te
                    elif chain_l_dev * 30. <= i < chain_l_dev * 70. and len(numpy.where(territor == int(str(te) + '3'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '3'
                        open_part = numpy.where(territor == int(str(te) + '3'))
                        territor[open_part]= te
                    elif chain_l_dev * 70. <= i < chain_l_dev * 125. and len(numpy.where(territor == int(str(te) + '4'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '4'
                        open_part = numpy.where(territor == int(str(te) + '4'))
                        territor[open_part]= te 
            
                    #print "AT_nr!!!", at_nr, i,  cur, type(cur)
                    mo, tang = rand_next(cur, state, chain, territor, nu, te, at_nr, i)
                    if tang:
                        try_chain +=1
                        assert try_chain != 100, "Chain is not possible to initialize"
                        at_nr = at_nr - i
                        old_state0 = state[tuple(chain[at_nr])]
                        territor[ter_zero]= te
                        state[ter_zero] = 0.
                        territor[ter_one]= int(str(te) + '1')
                        state[ter_one] = 0.
                        territor[ter_two]= int(str(te) + '2')
                        state[ter_two] = 0.
                        territor[ter_three]= int(str(te) + '3')
                        state[ter_three] = 0.
                        territor[ter_four]= int(str(te) + '4')
                        state[ter_four] = 0.
                        state[tuple(chain[at_nr])] = old_state0
                        state[state_lam] = LAMIN
                        cur = numpy.array(cur0)
                        print "Start one more time", try_chain, cur, at_nr
                        break
                    #print "MO", mo 
                    chain[at_nr] = cur + mo
                    print "Rosne", at_nr, chain[at_nr], nu, i, te
                    state[tuple(chain[at_nr])] = get_site_type(i, re, la)

                    cur = chain[at_nr]
                    if i == nu-1: # the last residue of the chain, so no need to remodelling the chain 
                        succ = True
        else:
            while try_chain < 100 and not succ:
                print try_chain, succ
                for i in range(nu-1, 0, -1):
                    chain_l_dev = round(nu/125.)
                    at_nr += 1
                    if chain_l_dev * 70. <= i < chain_l_dev * 125.  and len(numpy.where(territor == int(str(te) + '1'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '1'
                        open_part = numpy.where(territor == int(str(te) + '1'))
                        territor[open_part]= te
                    elif chain_l_dev * 30. <= i < chain_l_dev * 70. and len(numpy.where(territor == int(str(te) + '2'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '2'
                        open_part = numpy.where(territor == int(str(te) + '2'))
                        territor[open_part]= te
                    elif chain_l_dev * 10. <= i < chain_l_dev * 30. and len(numpy.where(territor == int(str(te) + '3'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '3'
                        open_part = numpy.where(territor == int(str(te) + '3'))
                        territor[open_part]= te
                    elif chain_l_dev <= i < chain_l_dev * 10. and len(numpy.where(territor == int(str(te) + '4'))[0]) != 0:
                        print "ZMIENIAM", str(te) + '4'
                        open_part = numpy.where(territor == int(str(te) + '4'))
                        territor[open_part]= te 
        
                    #print "AT_nr!!!", at_nr, i,  cur, type(cur)
                    mo, tang = rand_next(cur, state, chain, territor, nu, te, at_nr, i)
                    if tang:
                        try_chain +=1
                        assert try_chain != 100, "Chain is not possible to initialize"
                        at_nr = at_nr - nu + i
                        old_state0 = state[tuple(chain[at_nr])]
                        territor[ter_zero]= te
                        state[ter_zero] = 0.
                        territor[ter_one]= int(str(te) + '1')
                        state[ter_one] = 0.
                        territor[ter_two]= int(str(te) + '2')
                        state[ter_two] = 0.
                        territor[ter_three]= int(str(te) + '3')
                        state[ter_three] = 0.
                        territor[ter_four]= int(str(te) + '4')
                        state[ter_four] = 0.
                        state[tuple(chain[at_nr])] = old_state0
                        state[state_lam] = LAMIN
                        cur = numpy.array(cur0)
                        print "Start one more time", try_chain, cur, at_nr
                        break
                    #print "MO", mo 
                    chain[at_nr] = cur + mo
                    print at_nr, chain[at_nr], nu, i, te
                    state[tuple(chain[at_nr])] = get_site_type(i-1, re, la)

                    cur = chain[at_nr]
                    if i == 1: # the last residue of the chain, so no need to remodelling the chain 
                        succ = True
        
    mid = bound/2
    for i in range(m):
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
        binders[i] = [x, y, z]
        state[tuple(binders[i])] = BINDER

    return chain, binders, attached_to_lamins, state


def good_neighbors(x, i, chain, s_pos_chain, l_pos_chain):
    #check neighbors (chain)
    if i > 0 and (i not in s_pos_chain):
        #print "CHECK",s_pos_chain, i
        d1 = dist(chain[i - 1], x)
        if d1 > 0 and d1 < GOOD_NEIGH:
            pass
        else:
            return False
    if i not in l_pos_chain:
        d2 = dist(chain[i + 1], x)
        if d2 > 0 and d2 < GOOD_NEIGH:
            pass
        else:
            return False
    return True   


def bonds(chain, stat):

    bonds = 0
    for j in range(chain.shape[0]):
        molecule_pos = chain[j]
        molecule = stat[tuple(molecule_pos)]
        #print j
        if molecule == REGDNA:
            continue
 
        elif molecule == BSITE_R:
            binding = [BINDER]
        elif molecule == BSITE_L:
            binding = [LAMIN, BINDER]
        else: print "INNY STAN!!!", j, molecule
        #print j, molecule
        
        one_ch_at = 0
        for bmove in BMOVES: 
            new = molecule_pos + bmove
            enc = stat[tuple(new)]
            if enc in binding:
                bonds += 1
                #one_ch_at += 1
                #if one_ch_at > 6: break # one atom can have only 6 binding partners 

    return bonds

def modify(sta_pos_chain, la_pos_chain, chain, binders, state, bound = BOUND):
    #move binders
    
        
    if random.randint(0, 1):
        i = random.randint(0, len(binders) - 1)
        move = random.choice(MOVES)
        new = move + binders[i]
        #print "OLD_BIN", binders[i]
        if no_collisions(tuple(new), state):
            return False, i, move
    #move residues
    else:
        i = random.randint(0, len(chain) - 1)
        move = random.choice(MOVES)
        new = move + chain[i]
        #print "modify", i, move, chain[i], new

        #print "OLD", chain[i]
        if good_neighbors(new, i, chain, sta_pos_chain, la_pos_chain) and no_collisions(tuple(new), state):  # test if there is no collisions (the same place by different atoms) and no intersect of bonds
            if i not in sta_pos_chain and i not in la_pos_chain:
                if dist(chain[abs(i-1)], new) <= numpy.sqrt(2) and dist(chain[i+1], new) <= numpy.sqrt(2) and not intersect(new, chain[abs(i-1)], state, chain) and not intersect(new, chain[i+1], state, chain):
                    #print "Nie przecin", i
                    return True, i, move
                else: pass
                         
            elif i in la_pos_chain:
                #print "Last", i, la_pos_chain 
                if dist(chain[abs(i-1)], new) <= numpy.sqrt(2) and not intersect(new, chain[abs(i-1)], state, chain):
                #print "Nie przecin", i
                    return True, i, move
            elif i in sta_pos_chain:
                #print "FIRST", i, sta_pos_chain
                if dist(chain[i+1], new) <= numpy.sqrt(2) and not intersect(new, chain[i+1], state, chain):
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
def write_as_pdb(ch_nr, rev, chain, binders, attached_to_lamins, state, f, nb, metr_step, name = "chromosome and binders"):

    l = chain.shape[0]
    n = binders.shape[0]
    f.write("HEADER %i %d step %i\nTITLE %s" % (nb, l + n, metr_step, name))
    at_nr = 0

    def pdb_line(at_name, at_nr, res_nr, desc, pos, chain_n):
        return "\nATOM  " + str(at_nr).rjust(5) + " " + at_name.center(4) + " " + desc + " " + chain_n + str(res_nr).rjust(4) + "    " + str(round(pos[0] * DIST, 3)).rjust(8) + str(round(pos[1] * DIST, 3)).rjust(8) + str(round(pos[2] * DIST, 3)).rjust(8) + "  0.00 00.00"
        
    chain_list = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T"]
    for l, ch_n in zip(ch_nr, chain_list):
        if l not in rev:
            res_nu = 0
        else: res_nu = l+1
        for i in range(l):
            at_n = 'C'
            cur_chain = chain[at_nr]
            if state[tuple(cur_chain)] == REGDNA:
                r = "UNB"
            elif state[tuple(cur_chain)] == BSITE_R:
                r = "BOU"
            else: 
                r = "LAM" #BSITE_L
            at_nr += 1
            if l not in rev:
                res_nu += 1
            else: res_nu -= 1
            f.write(pdb_line(at_n, at_nr, res_nu, r, chain[at_nr-1], ch_n))
        f.write("\nTER")
            


    chain_at = at_nr
    res_nu = 0
    for i in range(n):
        res_nu += 1
        at_nr += 1
        r = "BIN"
        at_n = 'O'
        f.write(pdb_line(at_n, at_nr, res_nu, r, binders[i], "0"))

 
    ind = 0
    sum = ch_nr[ind]
    for i in range(1, chain_at):
        #print i, sum
        if i == sum:
            ind +=1
            sum = sum + ch_nr[ind]
            continue
        else:
            line = "\nCONECT" + str(i).rjust(5) +  str(i + 1).rjust(5)
        f.write(line) 
    f.write("\nEND   \n")

def count_bonds(pos, accepted, state):
    bonds = 0
    for bmove in BMOVES:
        if state[tuple(pos + bmove)] in accepted:
            if state[tuple(pos + bmove)] == LAMIN:
                bonds += 2
            else: bonds += 1
      #       if bonds > 6: break
    #print bonds
    return bonds

def radius_gyr(chai, last_pos):
    length = chai.shape[0]
    gyr_chains = []
    for la in last_pos:
        r = 0.0
        ind = last_pos.index(la)
        if ind == 0:
            for at in range(la):
                for ato in range(at, la):
                    r = r + math.sqrt(numpy.sum((chai[at] - chai[ato])**2))
            r_gyr = 2*r/(length*(length-1))
            gyr_chains.append(r_gyr)
        else:
            for at in range(last_pos[ind-1],la):
                for ato in range(at, la):
                    r = r + math.sqrt(numpy.sum((chai[at] - chai[ato])**2))
            r_gyr = 2*r/(length*(length-1))
            gyr_chains.append(r_gyr)
        
    return gyr_chains


DELTA=2
GYRATION = False
CHECK_E = False
def metropolis(nr_chrom, revers, chain, binders, attached_to_lamins, state, out_fname, name = "chromosome", n = 100):
    
    start_pos_chain = []
    po = 0
    for pos in nr_chrom:
        po = po+pos
        start_pos_chain.append(po)
        
    last_pos_chain = [s-1 for s in start_pos_chain]

    out_file = open(out_fname, "w")
    st_nr = 0
    E = bonds(chain, state)
    print "Starting energy:", E
    write_as_pdb(nr_chrom, revers, chain, binders, attached_to_lamins, state, out_file, st_nr, 0, name + ";bonds=" + str(E))

    for step in xrange(n):

        resp = modify(start_pos_chain, last_pos_chain, chain, binders, state)

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
                if state[tuple(old)] == BSITE_R:
                    Enew = E + count_bonds(ch[i], [BINDER], state) - count_bonds(old, [BINDER], state)
                    if CHECK_E:
                        Eslow = bonds(ch, st)
                        if Enew != Eslow:
                            print "R", 'Enew', Enew, 'Eslow', Eslow
                    else: pass
                elif state[tuple(old)] == BSITE_L:
                    Enew = E + count_bonds(ch[i], [LAMIN, BINDER], state) - count_bonds(old, [LAMIN, BINDER], state)
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
                Enew = E + count_bonds(b[i], [BSITE_R, BSITE_L], state) - count_bonds(old, [BSITE_R, BSITE_L], state) 
                if CHECK_E:
                    Eslow = bonds(ch, st)
                    if Enew != Eslow:
                        print "B", 'Enew', Enew, 'Eslow', Eslow
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
                if GYRATION:
                     print "iter", step, "step", st_nr, "energy:", E, "R_gyr ", radius_gyr(chain, last_pos_chain)
                else:
                    print "iter", step, "step", st_nr, "energy:", E
                
                write_as_pdb(nr_chrom, revers, chain, binders, attached_to_lamins, state, out_file, st_nr, step, name + ";bonds=" + str(E))
                #print "WRITE!!!"

    # dump the last state to the pickle
    l_obj = [chain, binders, attached_to_lamins, state, nr_chrom, revers]
    pickle_fname = out_fname.split('.pdb')[0] + ".pick"
    pickle_file = open(pickle_fname, 'w')
    pickle.dump(l_obj, pickle_file)

    out_file.close()


def output_name(ou, m, n):
    if ou == '':
        if type(m) == list:
            f_n = "MC_traj_%ibin_%ichain.pdb" % (len(m), len(n))
        else:
            f_n = "MC_traj_%ibin_%ichain.pdb" % (m, n)
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
    N = opts.Ch_lenght.split(',') # length of the chains
    N = map(int, N)
    R = opts.Revers.split(',')
    R = map(int, R)
    M = opts.Binders # nr of binders
    fn = output_name(opts.Out_str, M, N)
    c, b, a, state = initialize_random(N, M, fn, R)
else: 
    c, b, a, state, N, R = initialize_import(opts.In_str)
    #N = c.shape[0]
    M = b.shape[0]
    fn = output_name(opts.Out_str, M, N)

print "The lenght of the chain is ", N, ", the number of binders is ", M, ", the nucleus radius is ", R, ", the number of steps is ", opts.Steps, "random_seed ", #ran_seed

t2 = time.time()
print "initialization: ", t2 - t1
BOUND = numpy.max(c)

a = []
metropolis(N, R, c, b, a, state, fn, n = opts.Steps)
t1 = t2
t2 = time.time()
print "metropolis: ", t2 - t1
