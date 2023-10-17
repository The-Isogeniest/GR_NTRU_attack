from sage.all import *
import argparse
import random

import six

import sys, os
import re
from math import sqrt, floor, ceil, log, exp, log2

from random import randint
from multiprocessing import Pool
import csv
import warnings
import  numpy as np
import copy
from decimal import Decimal
max_beta = 70





def testAndMakeDir(path):
  if not os.path.isdir(path):
      os.makedirs(path)


def get_norm(vector):
    """
    Calculates and returns the norm of the vector.
    """
    # print("inside the function: \n")
    # print("iside vector: ", vector)
    if vector== None:
        return None
    s = 0
    for v in vector:
        s += v * v
    # print("inside the norm: ", np.sqrt(s))
    return sqrt(s)

def get_key_norm(n):
    """
    Input: n the order of the group
    Output: the norm of the key of the form (f,g)
    where f in T(d+1,d) and g in T(d,d).
    """

    d = int(n/3)
    return sqrt(4*d+1)

def is_it_ternary(l):
    """
    Input: a list l.
    Check if the list is ternary and returns True otherwise returns True
    """
    for i in l:
        if i!=1 and i!=-1 and i!=0 :
            return False
    return True
def is_it_zero(l):
    """
    Input: list l
    Output: True if the list entries are all zeros, False otherwise
    """
    for i in l:
        if i!=0:
            return False
    return True

def get_q_no_error(d,p=3):
    """
    The function returns the value of q that gives no decryption failure for variant of NTRU
    that has: h = gf^-1
    Input: d = int(order of the group/3)
           p usually 3
    """
    value= p*(6*d+1)
    q= 2**(len(bin(value))-2)
    return q




def is_it_pm_2(l):
    """
    Input: a list.
    Return True if all entries are two, minus two, or zeros.
    Otherwise: False.
    """
    for i in l:
        if i!=2 and i!=-2 and i!=0 :
            return False
    return True


def divide_by_2(l):
    """
    Input: a list of {2,-2,0}
        divide the coefficients by 2 and return the resultant list.

    """
    for i in range(len(l)):
        if l[i]>0:
            l[i] =1
        elif l[i]<0:
            l[i] = -1
    return l

def center_lift_form(f,q):
    """
    Centerlifting a vector f with respect to modulo q
    Input: f is a list
           q: a modulo
    Output: the centerlifting of the vector f with respect to q
    """
    t = f[:]
    for i in range(len(f)):
        t[i] = int(t[i])
        if t[i]>int(q/2):
            t[i] = t[i]-q
    return t

def add_vectors_with_centerlifting(l1, l2, n, q):
    """
    Coefficients-wise adding of coefficients of the correspondence vectors
    with centrelifting
    Input: l1, l2 two lists representing two vectors in the lattice.
    n: the vector length.
    Output: the resultant vector after adding them.
    """
    res = [0]*n
    for i in range(n):
        res[i] = (l1[i]+l2[i])%q
        if res[i]>int(q/2):
            res[i] = res[i]-q
    return res


def substract_vectors_with_centerlifting(l1, l2, n, q):
    """
    Coefficients-wise adding of coefficients of the correspondence vectors
    with centrelifting
    Input: l1, l2 two lists representing two vectors in the lattice.
        n: the vector length
    Output: the resultant vector after adding them.
    """
    res = [0] * n
    for i in range(n):
        res[i] = (l1[i] - l2[i]) % q
        if res[i] > int(q / 2):
            res[i] = res[i] - q
    return res


def dump_basis(B, filename, seed=None):
    """
    dumps integral basis B to a file named filename_seed#seedvalue.txt in the localpath/basis_dumps/
    """

    path = "basis_dumps/"
    testAndMakeDir(path)

    if not seed == None: filename += '_seed'+str(seed)+'.txt'
    else: filename += '.txt'

    d = B.nrows
    original_stdout = sys.stdout
    with open(path+filename, 'w') as f:
        sys.stdout = f
        #print('[')
        #for i in range(d):
        print(str(B))
        #print(']')

    sys.stdout = original_stdout
    return 1


def rough_estimate_on_betas(n,q, dihedral=False):
    """
    use the output of this function in case the use did not provide us with blocksizes
    """
    if dihedral:
        n = int(n/2)
    if n<150: beta_low = 10
    else: beta_low = floor( 0.28*4.*n/( (log(q)/log(n))**2 + 1))
    return list(range(beta_low, max_beta))

def get_next(l,n,upper_bound):
    """

    Input:  a list l represents incices to increase
            upper_bound: the upper bound allowed for increasing
            n: the list dimension
    The function adds one and increases the indices and update the list l
    """
    #l = [0,0,0,0,...0] #n zero
    if(l==[]):
        return [0]*n
    l[0]+=1
    i= 0
    while(i<n):
        if(l[n-1]==upper_bound):
            return "failure"
        if(l[i]<upper_bound):
            return l
        else:
            l[i] = l[i]%upper_bound
            l[i+1]+=1
            i+=1


def process_sample(sample):
    """

    Input: a sample as  a list of info
    a sample can be a list holding infor like: [seed, f, g, norm, h, f_prime1, f1_norm, f_prime2, f2_norm,  layer, beta1, beta2, filename,  total_time]

    The function processes the sample and returns a list of strings values corresponding to the values in the list.

    """
    l = []

    for i in range(len(sample)):

        s = str(sample[i])
        string_to_write = ""
        for i in range(len(s)):
            if (s[i] != "]" and s[i] != "["):
                string_to_write += s[i]
        l.append(string_to_write)
    return l


def  estimate_blocksize(n, q,layer=0):
    """
    Input: n half the order of the dihedral group.
           q: the modulo q used to build the lattice
          layer: 0 means no dimension reduction
                `1 means 1 layer of reduction
                2 means 2 layers of reduction
    Output: the blocksize estimated according to 2016-estimation needed
    to retrieve the key.
    """
    print("Please note that 2016-estimation is accurate for beta>=50 where beta<<d: the dimension of the lattice!!")
    print("n: {}, q:{}, layer:{}".format(n,q,layer))

    s_norm = get_key_norm(n)
    if layer==0:
        d_q = 2*n
        d = 4*n

    elif layer==1:
        d_q = n
        d = 2*n
        s_norm =sqrt(2)*s_norm  ##worst case estimation for the key norm for one layer attack
    else:
        d_q = int(n/2)
        d = n
        s_norm = 2*s_norm  # worst case estimation for the key norm for two layers attack



    for beta in range(30, 1000):
        left_side = Decimal(np.sqrt(beta/d)*s_norm)
        #print("left_side: ", left_side)
        delta = ((beta / (2 * np.pi * np.e)) * (np.pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))
        right_side = delta**(2*beta-d-1)*np.sqrt(q) ##np.sqrt(q) is
        #print("right side: ", right_side)
        if left_side<right_side:
            return beta
    return None


def create_file(seed, layer,filename):
    """
    Input: the seed, the layer and the file name
    The function creates a file with the specified path and write the header to the file.
    header = [seed, f, g, norm, h, f1_prime, f1_norm, f_prime2, f2_norm, beta1, beta2, total_time]
    """
    #print("file created: ")
    org_seed = seed
    seed = seed - (seed % 10 ** 8)
    path = "keys_dumps/layer_" +str(layer) + "/records/"
    testAndMakeDir(path)
    header1 = ['f', 'g', 'key norm', 'h', 'k1 (non-ternary)', 'k1-norm', 'k2 (ternary)', 'k2-norm', 'beta1', 'beta2', 'total time (seconds)']
    if layer ==1:
        header_cont= ['aux11 (non ternary) : (key, norm)', 'aux12  (non ternary): (key, norm)', 'aux21(ternary) : (key, norm)', 'aux22(ternary) : (key, norm)']
    elif layer ==2:
        header_cont = ['aux Mpp (layer 2)', 'aux Mpm(layer 2)', 'auxMmp(layer 2)', 'auxMmm(layer 2)','auxMp(layer 1)', 'auxMm(layer 1)' ]
    header = header1+header_cont


    filename += "_" + str(seed) + ".csv"
    isExisting = os.path.exists(path+filename)
    if not isExisting:
        with open(path + filename, "w", newline='') as wfl:
           csvwriter = csv.writer(wfl, delimiter=',')
           csvwriter.writerow([val for val in header])

def dump_seed(seed, layer,filename):
    """
    Input: the seed and the file name
    Output: write the seed to the file to later add the trails and the betas.
    """
    org_seed = seed
    seed = seed - (seed % 10 ** 8)
    path = "keys_dumps/layer_" + str(layer) + "/seeds/"
    testAndMakeDir(path)

    filename += "_" + str(seed) + ".txt"
    with open(path + filename, "a+") as f:
        print("seed: ",org_seed, file=f)
def dump_blocksize_for_layer(seed, filename, layer, sample):
  """
  Input:
        seed: the seed for it the key is generated.
        filename: the file name
        layer: how many layers of attack we have.

        sample:  a list contains [seed, f,g, norm, h, f_prime1, f1_norm, f_prime2, f2_norm,  layer, beta1, beta2, filename, total_time]

  """

  seed = seed - (seed % 10**8)
  path = "keys_dumps/layer_"+str(layer)+"/records/"
  testAndMakeDir(path)

  to_write = process_sample(sample)
  filename += "_" + str(seed) + ".csv"
  with open(path+filename, "a+", newline='') as wfl:
      csvwriter = csv.writer(wfl, delimiter=',')
      csvwriter.writerow([val for val in to_write])
      #print( str(beta1) + "\t" + str(beta2) + "\t" + str(total_time) + "\t"+ datetime.now().strftime('%Y-%m-%d %H:%M:%S'), file = f )





def parse_args():
    parser = argparse.ArgumentParser(description='Parse NTRU attack params.')

    #main parameters
    parser.add_argument('n', type=int, help="half the order of the dihedral group")
    parser.add_argument('--layer', type=int, dest="layer",  help="1 or 2 layers", default=1)
    parser.add_argument('--group', type=str, help="cyclic or dihedral", default="cyclic")
    parser.add_argument('-q', type=int, dest="q", default=None, help="NTRU modulus")
    parser.add_argument('--nsamples', type=int, default=None, dest="nsamples", help="Number of samples/rows of rot(h) used")
    parser.add_argument('--seed',  type=int, dest="seed", default=None, help="randomness seed")
    parser.add_argument('--h', dest="h", default=None, help="Uses given input as h, instead of creating a random instance.")
    parser.add_argument('--dump', dest='dump', default=False, help="flag to dump intermediate bases")
    # number of runs, number of threads
    parser.add_argument('-t', '--trials', type=int, dest="trials", default=1,
                        help="number of experiments to run per dimension")
    parser.add_argument('-w', '--workers', type=int, dest="workers", default=1,
                        help="number of parallel experiments to run")
    parser.add_argument('--threads', type=int, dest="threads", default = 1, help="number of threads used by 1 worker")


    parser.add_argument('--bkz_betas', type=str, dest="blocksizes", default=None, help="bkz block sizes as string of the form: min_beta:max_beta:step")
    parser.add_argument('--bkz_tours', type=int, dest="tours", default=8, help="number of tours of bkz reduction")


    parser.add_argument('--verbose', dest="verbose", default=False, help="verbosity")
    parser.add_argument('--dry-run', dest="dry_run", default=False,
                        help="Show parameters that would be used but don't run any actual experiments.")
    parser.add_argument('--show-defaults', dest="show_defaults", action='store_true',
                        help="Show default parameters and exit.")

    parser.add_argument('--filename', dest='filename', default=None, help="prefix of the dump filenames")
    parser.add_argument('--2016_estimation', dest ="2016_estimation", default=False, help="calculating the blocksize according 2016 estimation." )


    args, unknown = parser.parse_known_args()


    fmt = "{key:%ds}: {value}"%20

    if len(unknown)>0:
        print('Parameters', unknown, 'are not recognized and will be ignored')

    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    if args.show_defaults:
        for k, v in six.iteritems(all_defaults):
            print(fmt.format(key=k, value=v))
        exit(0)

    all_params = check_parsed_params(vars(args))

    if args.dry_run:
        for k, v in six.iteritems(all_params):
            print(fmt.format(key=k, value=v))
        exit(0)


    return all_params

def check_parsed_params(params):

    if params['q'] == None:
        ## for dihedral group, calculate for cyclic the power of two that gives
        ## error less than 2**-100, then q' = sqrt(2)*q for dihedral.

        n = params['n']
        d = int((2 * n) / 3)
        q = get_q_no_error(d)
        params['q'] = q

    if params['layer']==None and params['2016_estimation']:
        ## defualt for layer when 2016 estimation is enabled is layer=0
        ##i.e., doing the estimation for the big lattice without size reduction.
        params['layer'] = 0

    if not params['layer'] in [0,1,2]:
        raise ValueError("Please enter either 1 or 2 for layer")

    if params['2016_estimation']:
        return params
    if params['nsamples']==None: params['nsamples'] = params['n']
    else: assert(params['nsamples'] > 0 and params['nsamples']<=params['n'])



    if params['blocksizes'] ==None:
        params['blocksizes'] = rough_estimate_on_betas(params['n'], params['q'], dihedral= params['group']=="dihedral")
    else: params['blocksizes'] = eval("range(%s)" % re.sub(":", ",", params['blocksizes']))


    assert(len(params['blocksizes'])>0)



    if params['seed']==None:
        params['seed'] = randint(0, 2**64)



    if params['filename']==None:
        params['filename'] = str(params['n'])+'_'+str(params['q'])
    return params

def run_all(f, params):
    jobs = []
    if params['2016_estimation']:
        beta = estimate_blocksize(params['n'],params['q'], params['layer'])
        print("beta: ", beta)
        exit()
    original_seed = params['seed']
    if params['dump']:
        create_file(original_seed, params['layer'], params['filename'])
        #dump_seed(original_seed,params['group'],params['filename'])

    for t in range(params['trials']):
        params_  = copy.deepcopy(params)
        params_['seed'] = original_seed+random.randint(1,1000)
        jobs.append(params_)
    if params['workers'] == 1:
        for job in jobs:
            res = f(copy.deepcopy(job))
    else:
        pool = Pool(params['workers'])
        pool.map(f, jobs)
