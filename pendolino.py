import numpy,math,time, optparse, sys, pickle
import random as random

#accepted move vectors
MOVES=numpy.array([[0,0,1],[0,1,0],[1,0,0],[0,0,-1],[0,-1,0],[-1,0,0],[1,0,1],[0,1,1],[1,1,0],[-1,0,-1],[0,-1,-1],[-1,-1,0], [1,0,-1],[0,1,-1],[1,-1,0], [-1,0,1],[0,-1,1],[-1,1,0],[-1,1,-1],[-1,-1,-1],[1,1,1],[1,-1,1], [1,-1,-1], [-1, -1, 1], [-1, 1, 1], [1, 1, -1]])
#accepted matching positions of binding sites
BMOVES=numpy.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])

# radius of the nucleus
R=100
# 2 x radius + a fringe, because lamin barrier has to be hermetic
BOUND = 2 * R + 2
random.seed(1)

EMPTY=0
BINDER=1
LAMIN=2
BSITE_R=3
BSITE_L=4
REGDNA=5

# max distance of good neighbors
GOOD_NEIGH=3

def pars_inp():
    ## read/validate command-line arguments
    optparser=optparse.OptionParser(usage="%prog regular_dsites.txt lamin_bsites.txt [<options>]")

    optparser.add_option('-p',type="string",
        dest="Out_str",
        default='',
        help="An output filename with MC trajectory")
    optparser.add_option('-i',type="string",
        dest="In_str",
        default='',
        help="An input state (pickled file)")
    optparser.add_option('-s',type="int",
        dest="Steps",
        default=100000,
        help="Number of steps for simulation (default 100000)")
    optparser.add_option('-l',type="int",
        dest="Ch_lenght",
        default=512,
        help="Lenght of the chain (default 512)")
    optparser.add_option('-b',type="int",
        dest="Binders",
        default=256,
        help="Number of binders (default 256)")
			
    (opts, args)=optparser.parse_args()

    if len(sys.argv) <= 2:
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)
    return opts	

def init_dist_matrix(max_d=GOOD_NEIGH+1):
    dist_matrix=numpy.zeros((max_d, max_d, max_d),dtype=numpy.float32)
    for i in range(max_d):
        for j in range(max_d):
            for k in range(max_d):
                dist_matrix[i][j][k] = math.sqrt(i**2+j**2+k**2)
    return dist_matrix
DIST_MATRIX = init_dist_matrix()

def initialize_import(f):
    list_ob = pickle.load(open(f))
    ch = list_ob[0]
    b = list_ob[1]
    state = list_ob[2]
    return ch,b,state

def getStateWithLamins(bound):

    state = numpy.zeros((bound, bound, bound), dtype=numpy.int)
    MIDDLE = bound / 2

    def dist(x, y, z):
        return math.sqrt((x - MIDDLE) ** 2 + (y - MIDDLE) ** 2 + (z - MIDDLE) ** 2)

    for x in range(0, BOUND):
        for y in range(0, BOUND):
            for z in range(0, BOUND):
                if abs(dist(x, y, z) - MIDDLE + 1) <= 1:
                    state[x, y, z] = LAMIN

    return state

def initialize_random(n,m,bound=BOUND):
    chain=numpy.zeros((n,3), dtype=numpy.int)
    binders=numpy.zeros((m,3), dtype=numpy.int)
    state=getStateWithLamins(bound)

    chain[0, 0] = bound / 2
    chain[0, 1] = bound / 2
    chain[0, 2] = bound / 2

    def get_site_type_list(fpath, length):
        positions = [0] * length
        for l in open(fpath):
            positions[int(l)] = 1
        return positions

    regular_bsites = get_site_type_list(sys.argv[1], n)
    lamin_bsites   = get_site_type_list(sys.argv[2], n)

    def get_site_type(i, regular_bsites, lamin_bsites): # we should remember that it is possibility, that one site can interact with lamin and binders simultaneously
        if regular_bsites[i] == 1:
            site_type = BSITE_R
        elif lamin_bsites[i] == 1:
            site_type = BSITE_L
        else: site_type = REGDNA
        return site_type

    cur=chain[0]
    state[tuple(cur)] = get_site_type(0, regular_bsites, lamin_bsites)
    for i in range(1,n):
        mov=random.choice(MOVES)
        tries=0
        while tries<100 and not (no_collisions(tuple(cur+mov), state)):
            mov=random.choice(MOVES)
            tries+=1
        assert tries != 100, "unable to find initialization"
        chain[i]=cur+mov

        state[tuple(chain[i])] = get_site_type(i, regular_bsites, lamin_bsites)
        cur=chain[i]

    for i in range(m):
        cur=chain[i*n/m]
        mov=2*random.choice(MOVES)
        tries=0
        while tries<100 and not (no_collisions(tuple(cur+mov), state)):
            mov=2*random.choice(MOVES)
            tries+=1
        binders[i]=cur+mov
        state[tuple(binders[i])] = BINDER

    return chain,binders,state

def dist(r1,r2):
    return DIST_MATRIX[tuple(abs(r1-r2))]

#def rg(chain):
#   n=chain.shape[0]
#    res=0.0
#    for i in range(n):
#        for j in range(i):
#            res+=old_dist(chain[i],chain[j])
#    return 2*res/(n*(n-1))

def good_neighbors(x,i,chain):
    #check neighbors (chain)
    if i>0:
        d1=dist(chain[i-1],x)
        if d1>0 and d1<GOOD_NEIGH:
            pass
        else:
            return False
    if i< len(chain)-1:
        d2=dist(chain[i+1],x)
        if d2>0 and d2<GOOD_NEIGH:
            pass
        else:
            return False
    return True
    
def no_collisions(x, state):
    if state[tuple(x)] != EMPTY:
        return False
    return True

def bonds(chain,state):

    bonds=0
    for j in range(chain.shape[0]):
        molecule_pos=chain[j]
        molecule=state[tuple(molecule_pos)]

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

def modify(chain,binders,state,bound=BOUND):
    #move binders
    if random.randint(0,1):
        i=random.randint(0,len(binders)-1)
        move=random.choice(MOVES)
        new=move+binders[i]

        if no_collisions(tuple(new), state):
            return False, i, move
    #move residues
    else:
        i = random.randint(0,len(chain)-1)
        move=random.choice(MOVES)
        new=move+chain[i]
        if good_neighbors(new,i,chain) and no_collisions(tuple(new), state):
            return True, i, move
    return None

DIST=3
def write_as_xyz(chain,binders,f,name="chromosome and binders"):
    l=chain.shape[0]
    n=binders.shape[0]
    f.write("HEADER %d\nTITLE %s"%(l+n,name))
    at_nr = 0

    for i in range(l):
        if i%2:
            a="N"
            r = "BOU"
        else:
            a="C"
            r = "UNB"
        at_nr += 1
        line = "\nATOM  " + str(at_nr).rjust(5) + " " + a.center(4) + " " + r + "  " + str(at_nr).rjust(4) + "    " + str(round(chain[i,0]*DIST, 3)).rjust(8) + str(round(chain[i,1]*DIST, 3)).rjust(8) + str(round(chain[i,2]*DIST, 3)).rjust(8) + "  0.00 00.00"
        #f.write("\n%s %f %f %f"%(a,chain[i,0]*DIST,chain[i,1]*DIST,chain[i,2]*DIST))
        f.write(line)

    chain_at = at_nr
    for i in range(n):
        at_nr += 1
        a="O"
        r = "BIN"
        line = "\nATOM  " + str(at_nr).rjust(5) + " " + a.center(4) + " " + r + "  " + str(at_nr).rjust(4) + "    " + str(round(binders[i,0]*DIST, 3)).rjust(8) + str(round(binders[i,1]*DIST, 3)).rjust(8) + str(round(binders[i,2]*DIST, 3)).rjust(8) + "  0.00 00.00"
        f.write(line)
        #f.write("\n%s %f %f %f"%(a,binders[i,0]*DIST,binders[i,1]*DIST,binders[i,2]*DIST))

    ### output lamins

    for i in range(1, chain_at-1):
        line = "\nCONECT" + str(i).rjust(5) +  str(i+1).rjust(5)
        f.write(line)
    f.write("\nEND   \n")

def count_bonds(pos, accepted, state):
    bonds = 0
    for bmove in BMOVES:
        if state[tuple(pos+bmove)] in accepted:
             bonds += 1
    return bonds

DELTA=2
def metropolis(chain,binders,state,fn,name="chromosome",n=100):

    E=bonds(chain,state)
    traj=[(chain,binders,E)]

    f=open(fn,"w")
    write_as_xyz(chain,binders,f,name+";frame="+str(len(traj))+";bonds="+str(E))
    f.close()

    for step in range(n):

        resp=modify(chain,binders,state)

        ch=numpy.array(chain,copy=True)
        b=numpy.array(binders,copy=True)
        if resp:
            in_chain, i, move = resp
            if in_chain:
                old = numpy.copy(ch[i])
                ch[i] = ch[i] + move
                if state[tuple(old)] == BSITE_R:
                    Enew = E + count_bonds(ch[i], [BINDER], state) - count_bonds(old, [BINDER], state)
                elif state[tuple(old)] == BSITE_L:
                    Enew = E + count_bonds(ch[i], [LAMIN, BINDER], state) - count_bonds(old, [LAMIN, BINDER], state)
                else: # REGDNA
                    Enew = E
            else:
                old = numpy.copy(b[i])
                b[i] = b[i] + move
                Enew = E + count_bonds(b[i], BSITE_R, state) - count_bonds(old, BSITE_R, state)
        else:
            Enew = E

        if Enew>E or random.uniform(0.0,1.0)<math.exp((Enew-E)*DELTA): #accept

            if E != Enew:

                traj.append((ch,b,Enew))

                E=Enew
                f=open(fn,"a")
                print "iter",step, "energy:",E, "accepted:",len(traj)#,"Rg:",rg(chain)
                write_as_xyz(chain,binders,f,name+";frame="+str(len(traj))+";bonds="+str(E))
                f.close()

            chain=ch
            binders=b
            if resp:
                state[tuple(old + move)] = state[tuple(old)]
                state[tuple(old)] = EMPTY
            
    # load the last state to the pickle
    l_obj=[ch,b,state]
    fn = fn.split('.pdb')[0] + ".pick"
    file = open (fn, 'w')
    pickle.dump(l_obj, file)
    return traj

opts = pars_inp()

if opts.In_str == '':
    rand_init = True
else: 
    rand_init = False


N=opts.Ch_lenght # lenght of the chian
M=opts.Binders # nr of binders

if opts.Out_str == '':
    fn = "MC_traj_%ibin_%ichain.pdb" %(M,N)
elif "." in opts.Out_str:
    fn = opts.Out_str.split('.')[0] + ".pdb"
else: fn = opts.Out_str + ".pdb"


t1 = time.time()
if rand_init: 
    c,b,state=initialize_random(N, M)
else: 
    c,b,state = initialize_import(opts.In_str)
t2 = time.time()
print "initialization: ", t2 - t1
BOUND=numpy.max(c)

#print type(state), state.shape
#for i in state:
#    for j in i:
#        for k in j:
#            if k != 0:
#                print k

#fn="test_run_2-512-d3-rnd-prof_pendolino_test.pdb"
t=metropolis(c,b,state,fn,n=opts.Steps)
t1 = t2
t2 = time.time()
print "metropolis: ", t2 - t1
