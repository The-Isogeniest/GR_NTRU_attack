import time
from fpylll import BKZ as BKZ_FPYLLL, GSO, IntegerMatrix, FPLLL
from fpylll.tools.quality import basis_quality
from fpylll.algorithms.bkz2 import BKZReduction
from random import randint
import json
import DiTRU
from DiTRU import DiTRU
from utils import (get_key_norm, get_norm, is_it_zero, is_it_ternary, is_it_pm_2, add_vectors_with_centerlifting,
                   substract_vectors_with_centerlifting, divide_by_2, run_all, parse_args, dump_blocksize_for_layer,
                   get_next, dump_seed)

FPLLL.set_precision(120)
"""
Implementation of the lattice reduction for GR-NTRU based on the dihedral group for one and two layers of the attack.


How to:


Examples:

for dihderal group: 
    python attack.py 14  -q=128 --verbose=True --group="dihedral" --h="[115, 42, 117, 108, 73, 3, 53, 29, 108, 34, 72, 5, 36, 101]" --layer=1

"""


class attack:
    """
    Extension of Gentry’s attack to noncommutative group-ring NTRU over dihedral group
    We find the first block size that find the key for equivalent instances
    It includes the block size needed to find a {ternary key, non_ternary}

    Examples:
            - One-layer-attack: python attack.py 86 --layer=1  --bkz_betas=3:50 --seed=15721658425189707788
            - Two layer-attack: python attack.py 128 --layer=2  --bkz_betas=3:50 --verbose=True --dump=True --t=100
    """

    def __init__(self, params):

        #print("inside intialization")
        self.n = params['n']  # half the order of the dihedral group
        self.order = 2*self.n # the oder of the dihderal group
        self.q = params['q']  # The used modulo
        self.seed = params['seed']  # NTRU key seed
        self.layer = params['layer']  # one or two layers
        #self.nsamples = params['nsamples']  # number of NTRU samples
        self.blocksizes = params['blocksizes']  # range [a,b] with starting blocksize a, last blocksize b-1
        self.ntours = params['tours']  # number of bkz tours
        self.nthreads = params['threads']  # number of threads
        self.verbose = params['verbose']  # verbose mode
        self.filename = params['filename']  ##file name
        self.dump = params['dump']
        self.keynorm = get_key_norm(self.order)
        keyGenSuccess = False
        if params['h'] != None:
            self.h = json.loads(params['h'])  # for user-input h from NTRU Challenges
        else:
            self.h = None
        self.generator = DiTRU(self.n, self.q, self.seed, self.h,layer=self.layer)

        while not keyGenSuccess:
            try:
                self.keyseed = self.generator.newSeed()
                # print("seed: ", self.seed)
                # print(self.generator.get_key(self.keyseed))
                f,g,h,_,_= self.generator.get_key(self.keyseed)
                #print("f: ",f)
                #print("h: ", h)
                self.f = f
                self.g = g
                self.h = h
                #print("inside loop")
                self.lattices = self.generator.get_lattices(self.h)
                #print("after lattices")
                # self.lattice contains a tuple: as (lattice, plus_lattice, minus_lattice)
                # in the case of the cyclic group both of plus_lattice and minus_lattice are None
                #print("self.order: ", self.order)
                self.m = self.generator.get_random_message(8*self.order)
                #print("message: ", self.m)
                self.encrypted_message = self.generator.encrypt(self.m, self.h)
                #print("encrypted message: ", self.encrypted_message)
                keyGenSuccess = True
            except:
                ## exception happens if the key is not invertible
                self.seed = self.seed + 1
                self.generator.update_seed(self.seed)

        ### At this point we have the lattice of the element either in Z_qC_n or Z_qD_n
        #print("constructor: ", self.layer)
        if self.layer == 1:
            self.dim = self.order   ## the dimension where the lattice reducton algorithm is applied.
            self.threshold = 2*self.keynorm
            # threshold is two times of the key norm, (check in two smaller
            # lattices) and then pull back to the original lattice, will
            # get almost 4*||key||
        else:
            self.dim = self.n ##for two layers the reduction is applied on lattices with half the order of the dihedral
            self.threshold = self.keynorm
            self.upper_bound = 5 ##maximum number of vectors to check in the two layers attack
                                 ##for two layers attack the total number of comibinations to check will be upper_bound^4

        if self.dim <= 178:
            self.float_type = "long double"
        else:
            self.float_type = "mpfr"

        self.basis1 = None
        self.basis2 = None
        self.basis3 = None
        self.basis4 = None

        self.M1 = None
        self.M2 = None
        self.M3 = None
        self.M4 = None


        self.basis1 = IntegerMatrix.from_matrix(self.lattices[1], int_type="long")
        self.M1 = GSO.Mat(self.basis1, float_type=self.float_type,
                                U=IntegerMatrix.identity(self.basis1.nrows, int_type=self.basis1.int_type),
                                UinvT=IntegerMatrix.identity(self.basis1.nrows,
                                                            int_type=self.basis1.int_type))

        self.basis2 = IntegerMatrix.from_matrix(self.lattices[2], int_type="long")
        self.M2 = GSO.Mat(self.basis2, float_type=self.float_type,
                          U=IntegerMatrix.identity(self.basis2.nrows, int_type=self.basis2.int_type),
                          UinvT=IntegerMatrix.identity(self.basis2.nrows,
                                                       int_type=self.basis2.int_type))

        if self.layer==2:
            self.basis3 = IntegerMatrix.from_matrix(self.lattices[3], int_type="long")
            self.M3 = GSO.Mat(self.basis3, float_type=self.float_type,
                              U=IntegerMatrix.identity(self.basis3.nrows, int_type=self.basis3.int_type),
                              UinvT=IntegerMatrix.identity(self.basis3.nrows,
                                                           int_type=self.basis3.int_type))

            self.basis4 = IntegerMatrix.from_matrix(self.lattices[4], int_type="long")
            self.M4 = GSO.Mat(self.basis4, float_type=self.float_type,
                              U=IntegerMatrix.identity(self.basis4.nrows, int_type=self.basis4.int_type),
                              UinvT=IntegerMatrix.identity(self.basis4.nrows,
                                                           int_type=self.basis4.int_type))


    def __call__(self):
        self.progressive_search()  # call the function that retrieves the key

    def find_key(self):
        """
        For reduced basis, check if  (key/ternary key) exists
        The non-ternary key is accepted if its norm is smaller than the
        threshold.

        The ternary key is always accepted if exists.
        Apply one layer and two layers attack.

        """
        if self.layer==1:
            key_tuple, aux = self.check_for_one_layer()
        elif self.layer ==2:
            key_tuple = self.check_for_two_layers()

            #key_tuple = self.check_for_dihedral()

        return key_tuple

    def check_for_no_layer(self, keys_found_tuple):
        """
        This function is applicable on the original lattice, i.e, the lattice without
        applying any reduction.
        Upon a reduced basis of a lattice for GR-NTRU based on the dihedral group,
        check for the ternary/non-ternary key.
        keys_found_tuple: a tuple refers if the (non-ternary-found, ternary-found).
        Output: a tuple (k1, k2) where ki itself is a tuple as ([f,g], norm).
        if no key is returned, returns "failure".
        """

        key1_found = keys_found_tuple[0]
        key2_found = keys_found_tuple[1]
        k1 = None
        k2 = None
        norms = {}
        # print(self.M.B)
        B = self.M.B
        for i in range(self.dim):
            norms[i] = B[i].norm()
        sorted_norms = sorted(norms.items(), key=lambda x: x[1])
        for i in range(self.dim):
            if sorted_norms[i][1] > self.threshold:
                if key1_found or key2_found:
                    return (k1, k2)
                return "failure"

            fg = list(B[sorted_norms[i][0]])
            f = fg[self.n:]
            g = fg[:self.n]
            if not is_it_zero(g):
                if not key1_found and self.generator.is_invertible_R_p(f):
                    k1 = (fg, sorted_norms[i][1])  # (key, its norm)
                    key1_found = True

            if not key2_found and is_it_ternary(fg):
                if self.generator.is_invertible_R_p(f):
                    k2 = (fg, sorted_norms[i][1])  # (key, its norm)
                    key2_found = True

            if key1_found and key2_found:
                return (k1, k2)

        return "failure"

    def check_for_one_layer(self, keys_found_tuple):
        """
                Upon a reduced basis of a lattice for GR-NTRU based on the dihedral group,
                check for the ternary/non-ternary key.
                The reduced basis are two matrices of dimension n= 2d, where n is the order
                of the dihedral group.
                keys_found_tuple: a tuple refers of the (non_ternary found, ternary found)
                Output: a tuple (k1, k2) where ki itself is a tuple as ([f,g], norm).
                if no key is returned, returns "failure".
                along with (k1_aux k2_aux): the tuples of the smaller keys in the smaller lattices that helped
                retrieving the key in the larger lattice
        """
        #print("inside the one layer function: ")
        #print("key norm: ", self.keynorm)
        #print("threshold: ", self.threshold)

        key1_found = keys_found_tuple[0]
        key2_found = keys_found_tuple[1]
        k1 = (None, None)
        k2 = (None, None)
        k1_aux = ((None, None), (None,None))
        k2_aux = ((None, None), (None, None))


        d = int(self.n/2) ## half n( n is composite here n = 2d, the order of the dihedral = 2n)

        norms_plus = {}
        norms_minus = {}
        B_plus = self.M1.B  # reduced basis for the plus mat
        #print("B_plus[0]: ", list(B_plus[0]))
        B_minus = self.M2.B # reduced basis for the minus mat
        #print("B_minus[0]: ", list(B_minus[0]))
        for i in range(self.dim):
            norms_plus[i] = B_plus[i].norm()
        sorted_norms_plus = sorted(norms_plus.items(), key=lambda x: x[1])

        for i in range(self.dim):
            norms_minus[i] = B_minus[i].norm()
        sorted_norms_minus = sorted(norms_minus.items(), key=lambda x: x[1])

        for i in range(self.dim):
            if sorted_norms_plus[i][1] > self.threshold:
                if key1_found or key2_found:
                    return ((k1, k2), (k1_aux, k2_aux))
                return "failure"
            t1 = list(B_plus[sorted_norms_plus[i][0]])
            # print("N: ", N)
            g00 = t1[0:d]
            g01 = t1[d:2*d]
            f00 = t1[2*d:3*d]
            f01 = t1[3*d:4*d]
            if not is_it_zero(g00+g01):
                for j in range(self.dim):
                    if sorted_norms_minus[j][1] > self.threshold:
                        break
                    t2 = list(B_minus[sorted_norms_minus[j][0]])
                    g10 = t2[0:d]
                    g11 = t2[d:2*d]
                    f10 = t2[2*d:3*d]
                    f11 = t2[3*d:4*d]

                    if not is_it_zero(g10+g11):

                        fp0 = add_vectors_with_centerlifting(f00, f10, d, self.q)
                        fp1 = substract_vectors_with_centerlifting(f00, f10, d, self.q)
                        fp2 = add_vectors_with_centerlifting(f01, f11, d, self.q)
                        fp3 = substract_vectors_with_centerlifting(f01, f11, d, self.q)

                        gp0 = add_vectors_with_centerlifting(g00, g10, d, self.q)
                        gp1 = substract_vectors_with_centerlifting(g00, g10, d, self.q)
                        gp2 = add_vectors_with_centerlifting(g01, g11, d, self.q)
                        gp3 = substract_vectors_with_centerlifting(g01, g11, d, self.q)
                        # print("fp0: ", fp0)
                        # print("gp0: ", gp0)
                        F = fp0 + fp1+ fp2+fp3  # concatenating (fp0, fp1, fp2, fp3)
                        G = gp0 + gp1 +gp2+gp3 # concatenating  (gp0, gp1, gp2, gp3)

                        if not key1_found and self.generator.is_invertible_R_p(F):
                            # print("(f0,g0): ", f0+g0)
                            # print("first vecor norm: ", get_norm(g0 + f0))
                            # print("(f1,g1): ", f1+g1)
                            # print("second vecor norm: ", get_norm(g1 + f1))
                            k1 = (F + G, get_norm(F + G))
                            k1_aux = ((f00+f01+g00+g01, get_norm(t1)), (f10+f11+g10+g11, get_norm(t2)))
                            key1_found = True

                        if not key2_found and is_it_pm_2(F + G):
                            F = divide_by_2(F)
                            G = divide_by_2(G)

                            if self.generator.is_invertible_R_p(F):
                                k2 = (F + G, get_norm((F + G)))
                                k2_aux = ((f00 + f01 + g00 + g01, get_norm(t1)), (f10 + f11 + g10 + g11, get_norm(t2)))
                                key2_found = True

                            if key1_found and key2_found:
                                return ((k1, k2), [ k1_aux, k2_aux])
        # print("reached here")
        return "failure"

    def build_big_vector1(self, list_of_vectors):
       """
               Input: list_of_vectors for two layers attack, they are four vectors
               each vector from a reduced lattice in the second layer.
               In the last layer, we are building the higher lattices according to DCC and
               for the first layer, we add accordig to ACNS.
               Output: the big vectors built from the four smaller lattices and the auxiliary vectors that helped
                building them
        """



       index = 2
       aux = []
       result = []
       while len(list_of_vectors) > 1 or len(result)>1:
           #print("list of vectors: ", list_of_vectors)
           #print("index: ", index)
           #print("list of vecotrs: ", list_of_vectors)
           #print("0 ", list_of_vectors[0])
           #print("1", list_of_vectors[1])

           #d = int(len(t1) / 4)
           #n = d * 2
           if index == 0 or index==2:  ##index=2 for the second later, the positvie side ##0 for the first layer
               #print("inisde inedx: ", 0)
               if index ==2:
                   t1 = list_of_vectors.pop(0)
                   t2 = list_of_vectors.pop(0)
                   index = index-1
               else:
                   t1 = result.pop(0)
                   t2 = result.pop(0)
               d = int(len(t1)/4)
               g00 = t1[0:d]
               g01 = t1[d:2 * d]
               f00 = t1[2 * d:3 * d]
               f01 = t1[3 * d:4 * d]

               g10 = t2[0:d]
               g11 = t2[d:2 * d]
               f10 = t2[2 * d:3 * d]
               f11 = t2[3 * d:4 * d]

               fp0 = add_vectors_with_centerlifting(f00, f10, d, self.q)
               fp1 = substract_vectors_with_centerlifting(f00, f10, d, self.q)
               fp2 = add_vectors_with_centerlifting(f01, f11, d, self.q)
               fp3 = substract_vectors_with_centerlifting(f01, f11, d, self.q)

               gp0 = add_vectors_with_centerlifting(g00, g10, d, self.q)
               gp1 = substract_vectors_with_centerlifting(g00, g10, d, self.q)
               gp2 = add_vectors_with_centerlifting(g01, g11, d, self.q)
               gp3 = substract_vectors_with_centerlifting(g01, g11, d, self.q)
               # print("fp0: ", fp0)
               # print("gp0: ", gp0)
               F = fp0 + fp1 + fp2 + fp3  # concatenating (fp0, fp1, fp2, fp3)
               G = gp0 + gp1 + gp2 + gp3  # concatenating  (gp0, gp1, gp2, gp3)

               result.append(G+F)

               ### auxilarily
               aux.append(t1[2*d:]+t1[:2*d]) ## f,g
               aux.append(t2[2*d:]+t2[:2*d]) ## f,g


           else:
               t1 = list_of_vectors.pop(0)
               t2 = list_of_vectors.pop(0)
               #print("list of vectors after pop: ", list_of_vectors)
               #print("result: ", result)
               #print("t1: ", t1)
               #print("t2: ", t2)
               n = int(len(t1)/2)
               g0 = t1[:n]
               f0 = t1[n:]
               g1 = t2[:n]
               f1 = t2[n:]
               F0 = add_vectors_with_centerlifting(f0, f1, n, self.q)
               F1 = substract_vectors_with_centerlifting(f0, f1, n, self.q)
               G0 = add_vectors_with_centerlifting(g0, g1, n, self.q)
               G1 = substract_vectors_with_centerlifting(g0, g1, n, self.q)

               F = F0 + F1
               G = G0 + G1
               index = index-1
               result.append(G + F)

               aux.append(t1[n:]+t1[:n]) ##f,g

               aux.append(t2[n:]+t2[:n])  ##f,g




       # print("len inside the function: ", len(result[0]))
       #print("final: ", result)
       #print("result[0]: ", result[0])
       #x = get_x_vector(result[0], self.lattices[0])
       #print("to lay all should be integers: ", x)

       swapped = result[0][self.order:]+result[0][:self.order]
       #print("swapped: ", swapped)
       return swapped, aux


    def build_big_vector(self, list_of_vectors):
        """
        Input: list_of_vectors for two layers attack, they are four vectors
        each vector from a reduced lattice in the second layer.
        """

        #print("list of vecotrs: ", list_of_vectors)
        result = []
        index = 0
        while len(result)!=1:
            result = []
            while len(list_of_vectors) != 0:
                t1 = list_of_vectors.pop(0)
                t2 = list_of_vectors.pop(0)
                d = int(len(t1)/4)
                n = d*2
                if index == 0:
                    g00 = t1[0:d]
                    g01 = t1[d:2 * d]
                    f00 = t1[2 * d:3 * d]
                    f01 = t1[3 * d:4 * d]

                    g10 = t2[0:d]
                    g11 = t2[d:2 * d]
                    f10 = t2[2 * d:3 * d]
                    f11 = t2[3 * d:4 * d]

                    fp0 = add_vectors_with_centerlifting(f00, f10, d, self.q)
                    fp1 = substract_vectors_with_centerlifting(f00, f10, d, self.q)
                    fp2 = add_vectors_with_centerlifting(f01, f11, d, self.q)
                    fp3 = substract_vectors_with_centerlifting(f01, f11, d, self.q)

                    gp0 = add_vectors_with_centerlifting(g00, g10, d, self.q)
                    gp1 = substract_vectors_with_centerlifting(g00, g10, d, self.q)
                    gp2 = add_vectors_with_centerlifting(g01, g11, d, self.q)
                    gp3 = substract_vectors_with_centerlifting(g01, g11, d, self.q)
                    # print("fp0: ", fp0)
                    # print("gp0: ", gp0)
                    F = fp0 + fp1 + fp2 + fp3  # concatenating (fp0, fp1, fp2, fp3)
                    G = gp0 + gp1 + gp2 + gp3  # concatenating  (gp0, gp1, gp2, gp3)
                    index = 1
                elif index==1:
                    g0 = t1[:n]
                    f0 = t1[n:]
                    g1 = t2[:n]
                    f1 = t2[n:]
                    F0 = add_vectors_with_centerlifting(f0, g0,n,self.q)
                    F1 = substract_vectors_with_centerlifting(f0, f1, n, self.q)
                    G0 = add_vectors_with_centerlifting(g0, g1, n , self.q)
                    G1 = substract_vectors_with_centerlifting(g0, g1, n, self.q)

                    F = F0+F1
                    G = G0+G1
                result.append(F+G)
            list_of_vectors = result
        #print("len inside the function: ", len(result[0]))
        return result[0]

    def check_for_two_layers(self, keys_found_tuple):
        """
                Upon a reduced basis of a lattice for GR-NTRU based on the dihedral group,
                check for the ternary/non-ternary key.
                The reduced basis are four matrices of dimension d, where n= 2d is the order
                of the dihedral group.
                keys_found_tuple: a tuple refers of the (non_ternary found, ternary found)
                Output: a tuple (k1, k2) where ki itself is a tuple as ([f,g], norm).
                if no key is returned, returns "failure".
                along with (aux, aux2): aux1: (keypp, norm), (keypm, norm),
        """
        #print("inside the two layers function")
        key1_found = keys_found_tuple[0]
        key2_found = keys_found_tuple[1]
        k1 = (None, None) ##(key, norm)
        k2 = (None, None) ##(key, norm)
        #aux1 = ((None, None), (None, None)) ##layer's two pp, pm ((key, norm), (key,norm))
        #aux2 = ((None, None), (None, None)) ##layer's two mp, mm ((key, norm), (key,norm))
        #aux3 = ((None, None), (None, None))  ##layer 1 pulled back vectors from layer 2 (key, norm) for pulled back pp, pm, (key, norm ) for pulled back mp, mm
        #aux = [aux1, aux2, aux3]
        aux = []
        l = []
        n = 4##will be four for two layers attack.
        #print("key norm", self.keynorm)
        #print(list(self.M1.B[0]))
        #print(list(self.M2.B[0]))
        #print(list(self.M3.B[0]))
        #print(list(self.M4.B[0]))
        while True:
            list_of_vectors = []
            l = get_next(l,n,self.upper_bound)
            if l == "failure":
                if key1_found or key2_found:
                    return ((k1,k2), aux)
                return "failure"
            linv = l[::-1] ##inverse of l
            print(linv)


            list_of_vectors.append(list(self.M1.B[linv[0]]))
            list_of_vectors.append(list(self.M2.B[linv[1]]))
            list_of_vectors.append(list(self.M3.B[linv[2]]))
            list_of_vectors.append(list(self.M4.B[linv[3]]))

            big_vector, aux = self.build_big_vector1(list_of_vectors)
            #print("big_vector: ", big_vector)

            vector_norm = get_norm(big_vector)
            #print("vector norm: ", vector_norm)

            if vector_norm<=4*self.keynorm and not (key1_found):
                #print("here reached: ")
                #print("F: ", big_vector[:self.order])
                #print("is it invertible: ", self.generator.is_invertible_R_p(big_vector[:self.order]))
                if self.generator.is_invertible_R_p(big_vector[:self.order]):
                    aux_t =[]
                    for i in range(0,5,2):
                        aux_t.append( ( (aux[i], get_norm(aux[i])), (aux[i+1], get_norm(aux[i+1]))) )
                    k1 = (big_vector, vector_norm)
                    return ((k1, k2), aux_t)

            ####We have deleted the part that looks for ternary key, we couldn't find a reasonable way
            ####We have deleted the part that looks for ternary key, we couldn't find a reasonable way
            ## to retrieve the ternary key for two layers.
    def progressive_search(self):
        """
         Apply reduction algorithm with increased block sizes and return the block size
         that retrieves both a non-ternary and ternary keys
        """
        #print("progressive bkz")
        key1_found = False  # The non ternary key.
        key2_found = False  # The ternary key.
        beta = [0] * 2  ## block size needed to retrieve the (non-ternary key, ternary key)
        if self.layer==1:
            key_tuple = [(None, None), (None, None)]  ##[(non ternrary key, its norm), (ternary key, its norm)]
            aux1 = ((None, None), (None,None)) ## auxilary for non ternary, ternary and their respective norms
            aux2 = ((None, None), (None, None)) ## auxilary for ternary and their respective norms
            aux = [aux1, aux2]
        elif self.layer ==2:
            key_tuple = [(None, None), (None, None), (None, None)] ##for layer 2 two pairs, for layer one one pair (key, norm)
            aux1 = ((None, None), (None, None))  ## layer 2, for the positive side pair
            aux2 = ((None, None), (None, None))  ##layer 2, for the negative side pair
            aux3 = ((None, None), (None, None))  ##layer 1, the pulled back vectors
            aux = [aux1, aux2, aux3]


        T0_global = time.time()


        self.bkz1 = BKZReduction(self.M1)
        self.bkz1.lll_obj()

        self.bkz2 = BKZReduction(self.M2)
        self.bkz2.lll_obj()

        if self.layer==2:
            self.bkz3 = BKZReduction(self.M3)
            self.bkz3.lll_obj()

            self.bkz4 = BKZReduction(self.M4)
            self.bkz4.lll_obj()

        #print(self.M1.B)
        #print()
        #print(self.M2.B)
        #print()
        #print(self.M3.B)
        #print()
        #print(self.M4.B)

        if self.layer==1:
            result = self.check_for_one_layer((key1_found, key2_found))
        elif self.layer == 2: #two layers
            result = self.check_for_two_layers((key1_found, key2_found))

        if self.verbose:
            # print("group:", self.group)
            fmt = "{'initial LLL applied: layers: ':'%2d', 'total walltime': %.3f}"
            print(fmt % (self.layer, time.time() - T0_global))
            if result == "failure":
                print("failure")
            else:

                key_tuple_t = result[0]
                #print(key_tuple[0])
                #print(key_tuple[1])
                aux_t = result[1]
                if key_tuple_t[0][0] != None:
                    if self.verbose:
                        print("(Non ternary key, its norm)", key_tuple_t[0])
                    key1_found = True
                    beta[0] = 2
                    if self.layer ==1:
                        key_tuple[0] = key_tuple_t[0]
                        aux[0] = aux_t[0]

                    elif self.layer==2:
                        key_tuple = key_tuple_t
                        aux = aux_t

                if key_tuple_t[1][0] != None:
                    if self.verbose:
                        print(("(Ternary key ,its norm)", key_tuple_t[1]))
                    key2_found = True
                    beta[1] = 2
                    if self.layer==1:
                        key_tuple[1] = key_tuple_t[1]
                        aux[1] = aux_t[1]

        #print(key1_found)
        #print(key2_found)
        if not ((key1_found and key2_found) or (self.layer==2 and  key1_found)):

            for blocksize in self.blocksizes:  ##apply bkz with increasing block size
                T0_local = time.time()
                if self.verbose:
                    print("New round with block size: ", blocksize)

                for t in range(self.ntours):  # runs BKZ tours
                    par = BKZ_FPYLLL.Param(blocksize,
                                           strategies=BKZ_FPYLLL.DEFAULT_STRATEGY,
                                           max_loops=8)

                    self.bkz1(par)
                    self.bkz2(par)
                    if self.layer==2:
                        self.bkz3(par)
                        self.bkz4(par)

                    if self.layer ==1:

                        result= self.check_for_one_layer((key1_found, key2_found))


                    else:
                        result = self.check_for_two_layers((key1_found, key2_found))
                    if result == "failure":
                        print("failure")
                    else:
                        key_tuple_t = result[0]
                        aux_t = result[1]
                        if key_tuple_t[0][0] != None and not key1_found:
                            if self.verbose:
                                print("(Non ternary key, its norm)", key_tuple_t[0])

                            key1_found = True
                            beta[0] = blocksize
                            if self.layer == 1:
                                key_tuple[0] = key_tuple_t[0]
                                aux[0] = aux_t[0]

                            elif self.layer == 2:
                                key_tuple = key_tuple_t
                                aux = aux_t

                        if key_tuple_t[1][0] != None:
                            if self.verbose:
                                print(("(Ternary key ,its norm)", key_tuple_t[1]))
                            key2_found = True
                            beta[1] = blocksize
                            if self.layer == 1:
                                key_tuple[1] = key_tuple_t[1]
                                aux[1] = aux_t[1]
                    if (key1_found and key2_found) or (key1_found and self.layer==2):
                        break
                if  (key1_found and key2_found) or (key1_found and self.layer==2):
                    break
                fmt = "{'BKZ 'beta': %2d, number of layer':'%2d', 'total walltime': %.3f}"
                print(fmt % (blocksize, self.layer, time.time() - T0_local))
        print("Block size needed to find (non-ternary key, ternary key) is ({},{})".format(beta[0], beta[1]))

        t_end = time.time()-T0_global
        #print("len of key", len(key_tuple[0][0]))
        #print("f: ", key_tuple[0][0][:self.order])
        #print("len: ", len(key_tuple[0][0][:self.order]))
        #print("g: ", key_tuple[0][0][self.order:])
        #print("len: ", len(key_tuple[0][0][self.order:]))

        if key_tuple[0][1]!=None:
            #print(self.generator.ZDn_multiply(key_tuple[0][0][:self.order], self.h, self.q))
            decrypted_message1 = self.generator.decrypt(self.encrypted_message, key_tuple[0][0][:self.order])
            assert decrypted_message1 == self.m

        if key_tuple[1][1]!=None:
            decrypted_message1 = self.generator.decrypt(self.encrypted_message, key_tuple[1][0][:self.order])
            assert decrypted_message1 == self.m
        #print("original message: ", self.m)
        #print("decrypted message:", decrypted_message1)
        if self.dump:

            dump_seed(self.seed, self.layer, self.filename)
            ###[seed, f, g, norm, h, f1_prime, f1_norm, f_prime2, f2_norm, beta1, beta2, total_time]
            sample = [self.f, self.g, self.keynorm, self.h, key_tuple[0][0], key_tuple[0][1],
                      key_tuple[1][0],key_tuple[1][1], beta[0], beta[1], t_end]
            sample_cont = []
            if self.layer ==1:
                sample_cont = [aux[0][0], aux[0][1], aux[1][0], aux[1][1]]
            elif self.layer==2:
                sample_cont = [aux[0][0], aux[0][1], aux[1][0], aux[1][1], aux[2][0], aux[2][1]]
            sample = sample+sample_cont
            dump_blocksize_for_layer(self.seed, self.filename, self.layer, sample)



def main_call(params):
    attack_inst = attack(params)
    return attack_inst()


if __name__ == '__main__':
    print("main")
    all_parsed_params = parse_args()
    run_all(main_call, all_parsed_params)
