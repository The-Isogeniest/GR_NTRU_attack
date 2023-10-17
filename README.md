# Extension of Gentry's attack to non-commutative group-ring NTRU over dihedral group
This repository contains the scripts accompanying the article

## Providing one-layer attack and two-layer attack for composite parameters of group-ring NTRU over dihedral group 


# Requirements

* [fpylll](https://github.com/fplll/fpylll)

* [SageMath 9.5+](https://www.sagemath.org/) 


# Description of files
Short description of the content:
* `attack.py` main file to run lattice reduction attack on group-ring NTRU (GR-NTRU) over dihedral group.
* `DiTRU.py` contains the functions for key generation, encryption and decryption for GR-NTRU over dihedral group.
* `utils.py` contains helping functions
* folder `keys_dumps` contains two folders `layer_1` and `layer_2`: each of them contains two subfolders `records` and `seeds`
  * subfolder `records` contains records corresponding to the original keys and the retrieved keys according to the attack
    ; for one-layer-attack: the records saves `f`,`g`, `key norm`, `h` for the original key 
    `k1 (non-ternary)`: the non-ternary key retrieved by the attack, `k1-norm`: its norm,`k2 (ternary)`: the ternary key 
     retrieved by the attack,`k2-norm`: its norm,`beta1`: blocksize needed to retrieve k1, `beta2`: blocksize needed to
     retrieve the ternary key,`total time (seconds)`: total time of the attack,`aux11 (non ternary) : (key, norm)`,
    `aux12  (non ternary): (key, norm)`, `aux21(ternary) : (key, norm)`, `aux22(ternary) : (key, norm)`: the auxiliary 
     vectors in the smaller lattices that their pull-back gives the vector in the original lattice: aux11, aux12: the 
     auxiliaries for the non-ternary key and aux21, aux22: the auxiliaries for the ternary keys.
     For two-layers-attack: the records does not store info related to the related to the returned ternary key since it 
     is not returnable, and the auxiliary vectors are six, first four represent the four vectors lie in the second layer
     after applying the two-layers attack and last two vectors are the images of the previous four lie in the lattices
     of the first layer.
  * subfolder `seeds` contains the seeds for which the records have been generated; running the attack using the same 
    seed reproduces the same results.
# How to use

Run `attack.py` with the following parameters

* `n` : defines group-ring NTRU over dihedral group $Z_q D_n $. n: refers to half the order of  the dihedral group
  Integer. Obligatory parameter. The original lattice dimension before applying the attack is `4n`.
* `-q`: NTRU modulus. If not entered, it will be calculated automatically to be the first power of two to that guarantees
 no decryption failure.
* `--seed`: randomness seed to the to generate key and build the corresponding lattices.
* `--layer`: 0 for no dimension reduction, 1 to run one-layer attack, 2 to run two-layers attack.
* `--bkz_betas`: a string as 'b1:b2', where b1 indicates the first blocksize and b2 indicates the last blocksize
to try upon running progressive bkz.
* `--dump`: True to save the results into files, False otherwise.
* `--verbose`: True to print detailed output in the console while running the program, False otherwise.
* `--filename`: the name of the file where to save the results.

Parameters related to run experiments on parallel: 
* `-t` (or `--trials`): number of experiments to run per GR-NTRU dimension 
* `-w` (or `--workers`): number of parallel experiments to run
* `--threads`: number of threads used by 1 worker



# Experiments

To attack GR-NTRU over dihedral group for n=86 and one-layer of attack, you can run
```
python attack.py 86 --layer=1  --bkz_betas=3:50 --verbose=True --dump=True

```

It takes approximately less than two minutes on a laptop.
It generates a random seed and an instance corresponding to the seed and run the attack.

To generate a specific example or to reproduce one of our experiments, you can specify the seed


```
python attack.py 86 --layer=1  --bkz_betas=3:50 --seed=15721658425189707788
```


To run 100 trails and  two-layers attack against GR-NTRU with n=128, you can run


```
python attack.py 128 --layer=2  --bkz_betas=3:50 --verbose=True --dump=True --t=100
```
It takes approximately less than half a minute  on a laptop.




# Finding a short vector lies in a GR-NTRU lattice of a particular structure of dimension 1024.


Running the command 
```
python attack.py 256 --layer=2  --bkz_betas=15:50 --verbose=True --dump=True --seed=18265050100202766618 
```
finds a non-ternary key lies in a lattice of dimension 1024 with norm = 100.399 at blocksize= 23 
```commandline
0, -2, 2, 4, 2, 3, 0, 2, -1, -1, 2, -4, 1, -2, 4, -3, -1, 3, 3, 5, -3, -4, 1, -2, 1, -4, -4, 2, -4, -3, -1, 0, 0, 1, -1
, 0, -3, -3, 1, 3, 2, -4, 2, 4, -1, 1, 0, -3, -1, 4, -4, -4, -3, 1, 1, -11, 0, 2, 4, 4, 1, -1, 7, 1, 3, 2, 8, 3, -5, -4, 
6, -1, 0, -1, -4, 2, 0, 5, 6, -3, -4, 1, 4, 0, -3, 0, -6, -3, -4, -4, 2, 1, -2, -3, 4, 0, -5, -1, -3, 1, 3, 3, 2, -4, 0, 
-4, -1, 0, -2, 2, -6, -1, -1, -3, -4, -2, -4, 0, 0, -3, 0, -3, 0, 5, 3, -1, 1, -4, 0, -2, 2, -2, 4, 5, -2, -4, 1, -1, 4, 
-4, -5, 0, -8, 5, 5, 1, 3, 1, 1, -4, 1, 2, 9, 0, -4, -2, 0, -5, 3, 8, -4, 3, 1, -4, -1, -1, 3, -1, -4, 0, 2, 2, -3, 5, 2, 
-1, 1, 2, 2, 2, -3, 1, 5, 3, 6, 0, 0, 0, -1, -3, 3, 3, 5, -2, -4, 7, -1, -4, 0, -5, 4, -5, 2, 2, 0, 1, -2, -3, 0, 3, 2, 
-2, 1, -4, 0, -5, 2, -4, -2, -5, -2, -1, 2, 0, -3, -3, -1, 3, 1, 5, -6, -6, -2, 4, 1, 2, 2, 4, -4, -3, 1, 5, -2, -4, 2, -2, 
-2, 3, 2, -3, 8, -1, 1, 5, 1, 4, 4, 0, 4, -2, 1, -1, 2, 3, 1, -2, 3, 0, -1, 5, 0, 0, 2, -3, 3, -5, 1, -2, 3, 2, 0, 1, -3, 
4, -3, -1, 1, -9, 1, 5, 3, -2, 0, -3, -1, -1, -3, 2, -2, -2, -3, -1, -4, -3, 0, 1, -2, 4, 2, 1, -4, 5, 2, -5, -1, 1, -1, 
4, 3, -6, 7, -2, -8, -1, 2, -6, 4, 0, -2, -2, 7, 2, -8, -2, -4, 4, -1, -1, 0, -2, 5, -2, 0, 3, -1, 3, 1, -1, 3, 7, 4, -1,
0, 1, -1, -1, 0, 1, 0, -2, 3, 2, -3, 0, 0, 2, 8, 1, -2, 2, -2, 0, 1, -2, 5, -1, 2, 0, 1, 2, 1, -2, 1, -3, 0, -4, -4, -4, 
-1, 5, 0, -3, 3, -6, 5, 0, 1, 7, -4, 4, 0, 3, 3, -5, -3, -2, -1, 2, 0, -3, -3, 4, -7, 1, 1, -1, 1, -1, 5, -2, 2, -1, 1, 7,
-1, -2, -2, 0, -1, 3, 2, -1, -2, -1, 4, 6, 2, 1, 0, 3, 0, 1, 3, 1, 1, 2, 3, 0, 1, -2, -4, -5, 2, -6, 2, 0, 2, 2, -3, -2, 0, 
-2, 4, -4, -1, 1, -2, 4, 1, -2, 2, 1, -3, -1, -3, 5, -1, -3, 2, -1, -2, -1, 9, 1, -2, -1, 0, 0, 1, -6, -1, -6, 4, 0, -2, 7, 
-4, -2, -4, -2, 3, 0, -1, 5, 0, 0, 1, 0, -5, 0, 1, -3, -4, 2, 0, -2, -6, -1, -1, 4, 2, -5, 3, 5, 7, 3, 3, 3, -5, -2, -2, -4, 
3, -3, 0, -4, -3, -2, 1, -6, 5, 0, 2, -2, 1, 5, 7, -1, -2, 3, 0, 2, -2, 3, 4, 4, 1, 2, -1, -1, 1, 1, 2, 2, 5, -3, -2, 4, 5, 0,
0, -3, 3, -4, 5, -5, 4, -2, -1, 3, -2, -4, -3, -1, -2, 0, -1, 3, 1, 7, 1, -6, -1, -5, 0, -1, -3, -2, -1, -2, -4, 1, 0, 1, 1, 
2, 6, 5, -1, 3, -2, 0, 3, 4, -3, 0, 4, 4, -4, 4, -3, 0, -1, 1, -2, -1, 2, -3, -1, 5, 0, 3, -4, 0, -3, 0, 2, 1, -2, 1, -2, 2,
0, 2, 2, -3, -1, 0, -2, -3, 3, -1, 7, -1, -1, -7, 1, 0, -8, -2, -1, 3, 2, -4, -5, 2, 3, -4, 1, 0, 0, -2, -1, -3, 3, -3, 0, 7, 
8, -2, 0, -1, 6, 0, -1, 0, 1, 7, 3, -5, 4, 0, -5, 1, 2, 2, 3, -4, 2, -5, 3, -2, -3, -1, -6, -2, 1, 1, -2, 0, 9, -3, 2, -4, -1,
1, 1, -1, -3, -2, -3, -5, 2, -5, 1, 2, -1, -6, -4, -1, -4, -3, 1, 2, 0, -1, 1, -1, 0, -4, -5, 2, 3, 0, -2, -6, 6, 0, 3, 2, -3, 
5, 2, -3, 0, 5, 1, 1, 0, 3, -4, 4, 9, 0, 0, 1, -4, 1, 2, 0, -5, 1, -3, 2, -1, -1, 3, 1, -4, -3, 4, -2, -1, 1, -2, -2, 4, 0, -2, 
5, -2, 2, 3, -1, 2, -5, 3, 2, -1, 6, 0, -1, 2, 7, 8, -2, 3, -1, 0, 5, 0, -7, -2, 2, -7, -5, 3, -1, 0, -3, 1, -4, 0, -1, -2, 3, 
0, -2, 0, -2, 0, -1, 0, 4, 2, -4, -1, 1, 1, 0, 1, -2, -2, -1, 2, 2, 3, -2, 0, -1, -6, 1, -4, 8, -1, -4, 0, -2, 1, 2, -3, -2, 1, 
1, 2, -3, 1, 4, 1, -1, -2, -3, 6, 6, 4, 1, -2, 0, -3, -3, -4, 3, -6, 0, 1, 6, -6, 2, 3, -7, 5, -3, 1, -1, -1, -7, 4, 0, -5, 1, 
-3, 0, -1, -1, -5, 3, 0, 3, 0, 2, -5, 3, 0, 0, -2, 2, 2, 3, -4, 2, 5, -1, 4, 1, -1, -6, -7, -6, -2, -1, -2, 1, -6, 2, 3, 3, -2, 
1, 2, 1, -4, -4, -5, -1, -3, 5, 2, -1, 7, 4, 0, -3, 4, -1, -2, -4, 0, -4, 0, 3, 2, 0, 0, 2, 3, -3, 1, 2, -7, -6, 2, 1, 2, 2, -5,
2, 2, -5, 0, 5, 0, 0, -1, 0, 0, 4, 1, -2, 1, 2, 3, 3, -4, 1, 3, -4, 1, 5, 0, -3, 0, 0, 2, -3, -4, 6, 3, 5, 0, 1, 0, 0, -1, -2, 6,
-6, -1, 1, 1, 5, 3, 3, -3, -3
```
It took around 6.125 
hours on our device Intel(R) Xeon(R) CPU E3-1246 v3 @ 3.50GHz and 32 GB installed
RAM. 

## 2016 estimation

To get the blocksize estimated by [2016 estimation](https://www.usenix.org/system/files/conference/usenixsecurity16/sec16_paper_alkim.pdf)  to find the shortest vector, you can run
```commandline
python attack.py 256 --2016_estimation=True --layer=0

```
The previous command outputs the estimation of the blocksize needed to retrieve the ternary key
when no dimension reduction is applied i.e., lattice dimension is 1024.

To get the estimation for the blocksize in the case of one-layer or two-layers attack,
replace 0 by 1,2, respectively.
