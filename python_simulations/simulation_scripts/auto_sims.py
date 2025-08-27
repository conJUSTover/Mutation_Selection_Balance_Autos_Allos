#!/usr/bin/python

import numpy as np
from scipy import signal
import argparse

def random_mating(g):
    #takes in 3-element list of gamete frequencies "g" and returns a 1D convolution to produce genotype frequencies
#    print(g)
    return signal.convolve(g, g, method='direct')

def selection(G, fitness):
    #takes in 5-element list "G" and 5-element list "fitness"
    #calculates relative fitness and returns post-selection genotype frequencies 
    Gfit = [x*y for x,y in zip(G, fitness)]
    return [i/sum(Gfit) for i in Gfit]

def meiosis(G, seg):
    #takes in 5-element list of genotype frequencies "G" and 5-element list of segregation patterns "seg"
    #each segregation pattern in "seg" is the gamete frequencies produced from each genotype
    g0 = [seg[0][i] * G[0] for i in range(3)]
    g1 = [seg[1][i] * G[1] for i in range(3)] 
    g2 = [seg[2][i] * G[2] for i in range(3)]
    g3 = [seg[3][i] * G[3] for i in range(3)]
    g4 = [seg[4][i] * G[4] for i in range(3)]
    g = [sum(x) for x in zip(g0, g1, g2, g3, g4)]
    return g

def mutation(g, mu, nu):
    #takes in a 3-element list of gamete frequencies "g", forward mutation rate "mu", and backmutation rate "nu"
    #returns post-mutation gamete frequencies 
    g0 = g[0]*(1.0-mu)*(1.0-mu) + g[1]*nu*(1.0-mu) + g[2]*nu*nu
    g1 = g[0]*2.0*(1.0-mu)*mu + g[1]*(1.0-nu)*(1.0-mu) + g[1]*nu*mu + g[2]*2*(1.0-nu)*nu
    g2 = g[0]*mu*mu + g[1]*(1.0-nu)*mu + g[2]*(1.0-nu)*(1.0-nu)
    return [g0, g1, g2]

def set_seg(alpha):
    G0 = [1.0, 0.0, 0.0]
    G1 = [1/2 + alpha/4, 1/2 - alpha/2, alpha/4]
    G2 = [1/6 + alpha/3 , 2/3 - 2*alpha/3, 1/6 + alpha/3]
    G3 = [alpha/4, 1/2 - alpha/2, 1/2 + alpha/4]
    G4 = [0.0, 0.0, 1.0]
    return [G0, G1, G2, G3, G4]

def set_fitness(s, h1, h2, h3):
    return [1, 1-h1*s, 1-h2*s, 1-h3*s, 1-s]

def process_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", type=float, nargs='?', default = 0.0, const = 0.0, help="Selection Coefficient.")
    parser.add_argument("-h1", type=float, nargs='?', default = 0.0, const = 0.0, help="Dominance coefficient to modify selection acting on G1 individuals.")
    parser.add_argument("-h2", type=float, nargs='?', default = 0.0, const = 0.0, help="Dominance coefficient to modify selection acting on G2 individuals.")
    parser.add_argument("-h3", type=float, nargs='?', default = 0.0, const = 0.0, help="Dominance coefficient to modify selection acting on G3 individuals.")
   
    parser.add_argument("-mu", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Mutation Rate.")
    parser.add_argument("-nu", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Back Mutation Rate.")

    parser.add_argument("-alpha", type=restricted_float, nargs='?', default = 0.0, const=0.0, help="Rate of Double Reduction.")
 
    parser.add_argument("-N", type = int, nargs='?', default = 0, const = 0, help="Number of generations to run the simulation.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-g', action='store_true')
    group.add_argument('-G', action='store_true')

    parser.add_argument("-G0", type=restricted_float, nargs='?', default = 1.0, const = 1.0, help="Percent of Individuals with Genotype of 0 in initial generation.")
    parser.add_argument("-G1", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 1 in initial generation.")
    parser.add_argument("-G2", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 2 in initial generation.")
    parser.add_argument("-G3", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 3 in initial generation.")
    parser.add_argument("-G4", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 4 in initial generation.")

    parser.add_argument("-g0", type=restricted_float, nargs='?', default = 1.0, const = 1.0, help="Percent of Gametes with Genotype of 0 in initial generation.")
    parser.add_argument("-g1", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Gametes with Genotype of 1 in initial generation.")
    parser.add_argument("-g2", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Gametes with Genotype of 2 in initial generation.")

    parser.add_argument("-o", nargs='?', default='output_auto_sim.txt', help="Filename to print output of each generation. If not filename is provide, then output will be written to 'output_auto_sim.txt'.")

    args = parser.parse_args()
    return args

def generation(Gt0, segs, fitness, mu, nu):
    post_selection = selection(Gt0, fitness)
    post_meiosis = meiosis(post_selection, segs)
    g_t1 = mutation(post_meiosis, mu, nu)

    return g_t1

def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

def main():
    inargs = process_args()

    header = 'Generation,s,h1,h2,h3,mu,nu,alpha,g0,g1,g2,G0,G1,G2,G3,G4'

    if inargs.G:
        gams = [np.nan, np.nan, np.nan]
        genos = [inargs.G0, inargs.G1, inargs.G2, inargs.G3, inargs.G4]
        if abs(sum(genos) - 1.0) > 1e-6:
            raise Exception("Initial Genotype Frequencies must sum to 1.0, not ", sum(genos)) 

    if inargs.g:
        gams = [inargs.g0, inargs.g1, inargs.g2]
        if abs(sum(gams) - 1.0) > 1e-6:
            raise Exception("Initial Gamete Frequencies must sum to 1.0, not ", sum(gams))
        genos = random_mating(gams)

    fitness = set_fitness(inargs.s, inargs.h1, inargs.h2, inargs.h3)

    segs = set_seg(inargs.alpha)

    gen_stats = np.array([[x for y in [[0, inargs.s, inargs.h1, inargs.h2, inargs.h3, inargs.mu, inargs.nu, inargs.alpha], gams, genos] for x in y]])
    gams = generation(genos, segs, fitness, inargs.mu, inargs.nu)
    genos = random_mating(gams)

    for i in range(1, inargs.N):
        add_stack = [x for y in [[i, inargs.s, inargs.h1, inargs.h2, inargs.h3, inargs.mu, inargs.nu, inargs.alpha], gams, genos] for x in y]
#        print(add_stack)
        gen_stats = np.vstack([gen_stats, add_stack])
        gams = generation(genos, segs, fitness, inargs.mu, inargs.nu)
        genos = random_mating(gams)

    #print(generation_geno(init_genos, segs, fitness, inargs.mu, inargs.nu))
    np.savetxt(inargs.o, gen_stats, delimiter=',', header = header, fmt = ['%g', '%.3e', '%.3e', '%.3e', '%.3e', '%.3e', '%.3e', '%.3e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e'])

main()
