from random import Random
import math

import numpy as np
from sage.all import *
from utils import center_lift_form, get_q_no_error
"""
Implementation of GR-NTRU based on the dihedral group.

The class writes function to generate key, encryption and decryption.

"""


class DiTRU:

    """
    Create a DiTRU object with parameters
        -n: half the order  of the dihedral group (The order of the dihedral group is 2n).
        -q: the required modulo
        -seed: the seed to generate a random bit string.

    """

    def __init__(self, n, q, seed=None, h=None, layer=1):
        self.n = n
        self.q = q
        self.p = 3
        self.order = self.n*2
        if h==None:
            self.h = None
        else:
            self.h = h
        self.cache = {}
        self.seed = seed
        self.sample_fixed_type = 30 *self.order
        self.sample_key_bits = 2 * self.sample_fixed_type
        self.sample_iid_bits = 8*self.order
        self.d = math.floor((self.order/ 3))
        self.FFp = IntegerModRing(self.p)
        self.FFq = IntegerModRing(self.q)

        if seed == None:
            self.seed = randint(0, 2 ** 64)
        self.layer = layer

    """
        Input: Integer s.
        Output: Random bit array of length s.
    """

    def randomBitArray(self, s):

        random = Random()
        random.seed(self.seed)
        return [random.randrange(2) for i in range(s)]

    """
        set the seed as a random bit array of length sample_key_bits.
        To be used as seed in getKey() and getLattice().
    """

    def newSeed(self):
        self.seed = self.randomBitArray(self.sample_key_bits)
        return self.seed
    """
    Update the seed to a new seed
    """
    def update_seed(self,new_seed):
        self.seed =new_seed

    """
        Input: two elements representing two polynomials from the ring  Z_qDN =~Z_q[x,y]/(x^N-1, y^2-1, yx-x^(N-1)y)
        Output: the result of multiplication 
        """

    def ZDn_multiply(self, element1, element2, mod):
        # print("inside")
        multi_result = [0] * self.n*2
        n = self.n #half the order of the dihedral group.

        for i in range(n):
            for j in range(n):
                multi_result[(i + j) % n] = (multi_result[(i + j) % n] + element1[i] * element2[
                    j]) % mod

        for i in range(n):
            for j in range(n):
                multi_result[n + (j - i) % n] = (multi_result[n + (j - i) % n] + element1[i] * element2[
                    j + n]) % mod
        for i in range(n):
            for j in range(n):
                multi_result[n + (i + j) % n] = (multi_result[n + (i + j) % n] + element1[i + n] * element2[
                    j]) % mod
        for i in range(n):
            for j in range(n):
                multi_result[(-i + j) % n] = (multi_result[(-i + j) % n] + element1[i + n] * element2[
                    j + n]) % mod
        return multi_result

    """
               From  https://github.com/ElenaKirshanova/ntru_with_sieving
               based on https://ntru.org/f/ntru-20190330.pdf, Section 1.10.5.
               Input: A bit array b of length sample_fixed_type_bits.
               Output: A ternary polynomial with exactly d1 coefficients equal to 1 and d2 coefficients equal to −1.
    """

    def fixed_type(self, b, d1, d2):
        A = [0] * (self.order)
        v = [0] * (self.order)
        i = 0
        while i < d1:
            A[i] = 1
            for j in range(30):
                A[i] += 2 ** (2 + j) * b[30 * i + j]
            i += 1
        while i < d1 + d2:
            A[i] = 2
            for j in range(30):
                A[i] += 2 ** (2 + j) * b[30 * i + j]
            i += 1

        while i < self.order:
            for j in range(30):
                A[i] += 2 ** (2 + j) * b[30 * i + j]
            i += 1

        A.sort()

        for i in range(self.order):
            v[i] = A[i] % 4
            if v[i] == 2:
                v[i] = -1

        return v

    """
        From  https://github.com/ElenaKirshanova/ntru_with_sieving
        based on https://ntru.org/f/ntru-20190330.pdf, Section 1.10.3
        Input: A bit array b of length sample_iid_bits.
        Output: A ternary polynomial.
        
      """

    def ternary(self, b):
        v = [0] * (self.order)

        for i in range(self.order):
            coeff_i = 0
            for j in range(8):
                coeff_i += 2 ^ j * b[8 * i + j]
            v[i] = coeff_i

        for i in range(self.order):
            v[i] = v[i]%3
            if v[i]==2:
                v[i]=-1
        return v

    """
    get a random message, call the function above and samples 
    a random ternary message.
    Input: len: the length of the message(usually len is sample_iid_bits)
    Output: message: the random generated message.
    """
    def get_random_message(self, len):
        random = Random()
        random.seed(randint(0,2**64))
        random_bit_stream= [random.randrange(2) for i in range(len)]
        message = self.ternary(random_bit_stream)
        return message


    """
               Input: arr: an array , n: an integer 
               Output: shifting to left the array by n positions
    """

    def shiftLbyn(self, arr, n=0):
        return arr[n::] + arr[:n:]

    """

        Matrix representation for an element of Z_qC_n (right circulant matrix)
        It's also an auxiliary matrix for matrix representation for an element in Z_qD_n
        Input: the first row that represents an element f,g, or h.
        FF: the space over it, the matrix to be constructed either IntegerModRing(3)
        or IntgerModRing(q)


        Output: the matrix that represents the group ring element for a cyclic group
        """

    def get_A(self, first_row, FF):

        n = len(first_row)
        a = first_row
        m = []
        for i in range(n):
            m.append(a)
            a = self.shiftLbyn(a, -1)

        MS2 = MatrixSpace(FF, n, n)
        B = MS2.matrix(m)

        return B

    """
    An auxiliary matrix used to represent an element in Z_qD_n
    Input: the first row that represents an element f, g or h.
    FF: the space over it, the matrix to be constructed either IntegerModRing(3)
    or IntgerModRing(q)

    Output: the matrix that represents the group ring element for a cyclic group
            the output is a left circulant matrix.

    """

    def get_B(self, first_row, FF):
        n = len(first_row)
        a = first_row
        m = []
        for i in range(n):
            m.append(a)
            a = self.shiftLbyn(a, 1)

        MS2 = MatrixSpace(FF, n, n)
        B = MS2.matrix(m)

        return B

    """
    Input: an element of Z_qD_n.
    FF: the space over it, the matrix to be constructed either IntegerModRing(3)
    or IntgerModRing(q).
    Output: the matrix representation of the element.
    and the plus and minus upper parts of the matrices to
    be used in the reduction for one layer/ two layers of the attack.
    For one layer we return two matrices.
    For two layers of the attack we return four matrices.
    """

    def Zdn_matrices(self, element, FF, layer=1):

        n = self.n   # half the order for dihedral group.
        d = int(n/2)

        H0 = self.get_A(element[:n], FF)
        H00 = H0[0:d, 0:d]
        H01 = H0[0:d, d:n]
        H00pH01 = H00+H01
        H00mH01 = H00-H01

        H1 = self.get_B(element[n:], FF)
        H10 = H1[0:d, 0:d]
        H11 = H1[0:d, d:n]
        H10pH11 = H10+H11
        H10mH11 = H10-H11


        # print("B: ", B)
        M = block_matrix(2, 2, [H0, H1, H1, H0])
        Mh2dp_upper_right_one_layer = block_matrix(2, 2, [H00pH01, H10pH11, H10pH11, H00pH01])
        #print("h: ", element)
        #print("one layer p: \n", Mh2dp_upper_right_one_layer)
        ## For positive part, we apply ACNS
        if layer==2:
            #DCC
            #Mh2pp_upper_right_two_layers = H00pH01+ H10pH11
            #Mh2pm_upper_right_two_layers = H00pH01- H10pH11

            # ACNS
            d_p = int(d/2)
            A = H00pH01[0:d_p, 0:d_p]
            B = H00pH01[0:d_p, d_p:d]

            ApB = A+B
            AmB = A-B

            C = H10pH11[0:d_p, 0:d_p]
            D = H10pH11[0:d_p, d_p:d]

            CpD = C+D
            CmD = C-D
            Mh2pp_upper_right_two_layers = block_matrix(2, 2, [ApB, CpD, CpD, ApB])
            #print("pp two layer \n", Mh2pp_upper_right_two_layers )
            Mh2pm_upper_right_two_layers = block_matrix(2, 2, [AmB, CmD, CmD, AmB])
            #print("pm two layer \n",Mh2pm_upper_right_two_layers )


        Mh2dm_upper_right_one_layer = block_matrix(2, 2, [H00mH01, H10mH11, H10mH11, H00mH01])
        #print("one layer m: \n ",  Mh2dm_upper_right_one_layer )
        ## For minus part, we apply DCC paper
        if layer ==2:
            Mh2mp_upper_right_two_layers = H00mH01+H10mH11
            Mh2mm_upper_right_two_layers = H00mH01-H10mH11
            #print("mp two layer: \n ", Mh2mp_upper_right_two_layers )
            #print("mm two layer: \n", Mh2mm_upper_right_two_layers)
        if layer ==1:
            #print("M: ", M)
            #print("Mh2dp: ", Mh2dp_upper_right_one_layer)
            #print("Mh2dp: ", Mh2dm_upper_right_one_layer)
            toreturn = (M, Mh2dp_upper_right_one_layer, Mh2dm_upper_right_one_layer)
        elif layer==2:
            toreturn = (M,Mh2pp_upper_right_two_layers,  Mh2pm_upper_right_two_layers, Mh2mp_upper_right_two_layers, Mh2mm_upper_right_two_layers)

        return toreturn

    """
       Returns the matrix corresponding to the element,
       
    """

    def element_to_matrix(self, element, FF):
        H0 = self.get_A(element[:self.n], FF)
        H1 = self.get_B(element[self.n:], FF)
        M = block_matrix(2, 2, [H0, H1, H1, H0])
        return M

    def get_Fp(self, f):
        """
        Input: f a vector used as a decryption key.
        Output: Fp the inverse of f mod (D_n, p)
        """
        Fp_mat = self.element_to_matrix(f, self.FFp)
        return Fp_mat.inverse()[0]


    """
       Sample f, g corresponding to the initial variant of NTRU defined in  Hoffstein  book 

       Input: A bit array fg_bits of length sample_key_bits.
       Output: two vectors represent f,g.
       """



    def sample_fg(self, fg_bits):
        #print("fg_bits: ", len(fg_bits))
        f_bits = fg_bits[0:self.sample_fixed_type]
        g_bits = fg_bits[self.sample_fixed_type:]

        f = self.fixed_type(f_bits, self.d + 1, self.d)
        #print("f: ", f)
        Fp_mat = self.element_to_matrix(f, self.FFp)  ## Z_qD_n matrix
        #print("Fpmat:", Fp_mat)
        if Fp_mat.is_invertible():
            Fq_mat= self.element_to_matrix(f, self.FFq)

            if Fq_mat.is_invertible():
                Fp = Fp_mat.inverse()[0]  # inverse of f mod (p, R(p,DN))
                Fq = Fq_mat.inverse()[0]  # inverse of f mod (q, R(q,DN))
                g = self.fixed_type(g_bits, self.d, self.d)
                return (f, g, Fp, Fq)
        ## We reach here if the inverse doesn't exist for the previous seed
        ## Therefor, we generate a new seed
        # print("generating a new seed")
        # self.seed = randint(0, 2 ** 64)
        raise ValueError("Not invertible key")

    """
    Input: an element of the underlying group-ring.
    The function builds the matrix of the group ring element and 
    return True if it's invertible, otherwise, it returns False.
    """

    def is_invertible_R_p(self, element):
        Fp_mat = self.element_to_matrix(element, self.FFp)  ## Z_qD_n matrix
        if Fp_mat.is_invertible():
            return True
        return False

    """
        Input: seed
        Output: a key (f,g,h) where h = g*f^-1 mod(q, X^n-1)
    """

    def get_key(self, seed):
        if self.h != None:  ##h has been initialized in the constructor
            return (None, None, self.h)

        seedT = tuple(seed)
        if seedT in self.cache:
            return self.cache[seedT]

        else:
            f, g, Fp, Fq = self.sample_fg(seed)

            h = self.ZDn_multiply(Fq, g, self.q)
            # print("h",h)

        self.cache[tuple(self.seed)] = (f, g, h, Fp, Fq)
        return (f, g, h, Fp, Fq)

    """
        Input: seed
        Output: Coppersmith-Shamir basis for dihedral group.
        For a dihedral group  of order 2n where n is a composite number = 2d 
        for prime d, the output contains three lattices : The original lattice 
        and the other two lattices for reduction {plus, minus} lattices.
        """

    def get_lattices_one_layer(self, h):
        """
        Generate Coppersmith-Shamir lattices for one layer of the attack.
        Input: h the public key

        Output: The returned lattices are the big lattice for no reduction,
        plus basis, minus basis: corresponding to one layer of the attack.
        remember n: half the order of the lattice here is a 2d where d is a prime.

        """
        order = self.order
        n = self.n
        q = self.q
        plus_basis = None
        minus_basis = None

        B = [[0] * 2*order for i in range(order*2)]
        for i in range(order):
            B[i][i] = q
        for i in range(order):
            B[order + i][order + i] = 1

        element_mat, plus_mat, minus_mat = self.Zdn_matrices(h, self.FFq)
        for i in range(order):
            for j in range(order):
                B[order + i][j] = int(element_mat[i][j])

        plus_basis = [[0] * order for i in range(order)]
        ## remember n is half the order of the dihedral group
        for i in range(n):
            plus_basis[i][i] = q
        for i in range(n):
            plus_basis[n + i][n + i] = 1

        for i in range(n):
            for j in range(n):
                plus_basis[n + i][j] = int(plus_mat[i][j])

        # Construct minus basis
        ## The reduced lattice dimension is half of that for the original lattice.
        minus_basis = [[0] * order for i in range(order)]
        for i in range(n):
            minus_basis[i][i] = q
        for i in range(n):
            minus_basis[n + i][n + i] = 1

        for i in range(n):
            for j in range(n):
                minus_basis[n + i][j] = int(minus_mat[i][j])

        return (B, plus_basis, minus_basis)


    """
    Input: seed
    Output: Coppersmith-Shamir basis for dihedral group.
    For a dihedral group the output contains five lattices : The original lattice 
    and four lattices for reduction {plus_plus, plus_minus, minus_plus, minus_minus} lattices.
    """

    def get_lattice_two_layers(self, h):
        """
        Input: h the public key
        Generate Coppersmith-Shamir lattices, the generated lattices are four lattices for the two layers attack.
        """
        order = self.order
        # The order is 2n where n itself is a composite number n = 2d for d composite (not prime).
        n = self.n
        d = int(n/2)
        q = self.q
        plus_plus_basis  = None
        plus_minus_basis  = None
        minus_plus_basis = None
        minus_minus_basis = None


        B = [[0] * 2 * order for i in range(order * 2)]
        for i in range(order):
            B[i][i] = q
        for i in range(order):
            B[order + i][order + i] = 1

        element_mat, plus_plus, plus_minus, minus_plus, minus_minus = self.Zdn_matrices(h, self.FFq, layer=2)

        for i in range(order):
            for j in range(order):
                B[order + i][j] = int(element_mat[i][j])
        #order = 2n, therefore the original lattice dim is 4n and the reducced lattices for two layers is 2n.
        # The original lattice dimension is 4times the reduced lattices for two layers of attack.

        ################## plus plus ####################
        plus_plus_basis = [[0] * n for i in range(n)]
        ## remember n is half the order of the dihedral group
        for i in range(d):
            plus_plus_basis[i][i] = q
        for i in range(d):
            plus_plus_basis[d + i][d + i] = 1

        for i in range(d):
            for j in range(d):
                plus_plus_basis[d + i][j] = int(plus_plus[i][j])



        ################## plus minus ####################
        plus_minus_basis = [[0] * n for i in range(n)]
        ## remember n is half the order of the dihedral group
        for i in range(d):
            plus_minus_basis[i][i] = q
        for i in range(d):
            plus_minus_basis[d + i][d + i] = 1

        for i in range(d):
            for j in range(d):
                plus_minus_basis[d + i][j] = int(plus_minus[i][j])

        ################## minus plus ####################
        minus_plus_basis = [[0] * n for i in range(n)]
        ## remember n is half the order of the dihedral group
        for i in range(d):
            minus_plus_basis[i][i] = q
        for i in range(d):
            minus_plus_basis[d + i][d + i] = 1

        for i in range(d):
            for j in range(d):
                minus_plus_basis[d + i][j] = int(minus_plus[i][j])

        ################## minus minus ####################
        minus_minus_basis = [[0] * n for i in range(n)]
        ## remember n is half the order of the dihedral group
        for i in range(d):
            minus_minus_basis[i][i] = q
        for i in range(d):
            minus_minus_basis[d + i][d + i] = 1

        for i in range(d):
            for j in range(d):
                minus_minus_basis[d + i][j] = int(minus_minus[i][j])
        #print("big lattice: ", self.auxilary_func(B))
        #print("h: ",h)
        #print("plus plus ", self.auxilary_func(plus_plus_basis))
        #print("plus minus ", self.auxilary_func(plus_minus_basis))
        #print("minus plus ", self.auxilary_func(minus_plus_basis))
        #print("minus minus ", self.auxilary_func(minus_minus_basis))
        return (B, plus_plus_basis, plus_minus_basis, minus_plus_basis,  minus_minus_basis)

    def auxilary_func(self,mat):
        """
        Input mat and print the mat in friendly way for redding
        """
        length = len(mat[0])
        for i in range(length):
            print(mat[i])
            #print()

    """
    This function returns the lattices to apply the attack on.
    Input: h the public key, for which, we are generating the keys and the corresponding lattices.
    Output: Either calling the get_lattice_one_layer() and return two lattices or get_lattice_two_layers()
    and return four lattices.
    """
    def get_lattices(self,h ):
        if self.layer == 1:
            return self.get_lattices_one_layer(h)
        elif self.layer == 2:

            return self.get_lattice_two_layers(h)

    """
    Input: a message to encrypt
           h: the public key
    Output: the encrypted message
    """
    def encrypt(self, message, h):

        random = Random()
        random.seed(randint(0, 2 ** 64))
        seed_for_r = [random.randrange(2) for i in range(self.sample_fixed_type)]
        r = self.fixed_type(seed_for_r,self.d, self.d)
        e1 = self.ZDn_multiply(h, r, self.q)
        prh = list(np.multiply(self.p, e1))
        e = list(np.add(prh, message))

        return e

    """
    Input: an encrypted message 
            f: the private key 
            Fp: the inverse of the private key with respect to mod p
    Output: the decrypted message
    """
    def decrypt(self, encryptedmessage,f, Fp=None):
        if Fp==None:
            Fp= self.get_Fp(f)
        fe = self.ZDn_multiply(f, encryptedmessage, self.q)
        a =  center_lift_form(fe,self.q)
        decrypted_q = self.ZDn_multiply(Fp, a,self.p)
        message_prime = center_lift_form(decrypted_q, self.p)
        return message_prime
def main():

    n = 2*2*13
    order = 2*n
    d = int(order/3)
    q = get_q_no_error(d)
    print("q: ",q)
    #print(is_it_ternary([1, -1, 0, -1, 1, 0, 0, -1, 0, 0, 1, 1, 0, 0]))
    for i in range(100):
        try:
            seed = randint(0,2**64)
            keygen =DiTRU(n, q,seed=seed)
            keyseed = keygen.newSeed()
            f,g,h, Fp, Fq= keygen.get_key(keyseed)
            #print("f: ", f)
            #print("g: ", g)
            #print("h: ", h)
            #print("Fp: ", Fp)
            #print("Fq: ", Fq)
            m = keygen.get_random_message(keygen.sample_iid_bits)
            encrypted_message = keygen.encrypt(m,h)
            decrypted_message = keygen.decrypt(encrypted_message,f,Fp)
            assert decrypted_message == m
            #print("original message: ", m)
            #print("decrypted message: ", decrypted_message)
            L1, L2 , L3, L4, L5 = keygen.get_lattice_two_layers(keyseed)
            #print("L1:  ", L1)
            #print("L2:  ", L2)
            #print("L3: ", L3)
            #L1, L2, L3= keygen.get_lattice(keyseed)
            #print("L1: ", L1)
            #print("L2: ", L2)
            #print("L3: ", L3)
        except:
            seed+=1



if __name__ == "__main__":
    main()
