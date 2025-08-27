#!/usr/bin/python

import numpy as np
from scipy import signal
import argparse

def random_mating(g):
    #takes in 4-element list of gamete frequencies "g" and returns a 2D convolution to produce genotype frequencies
    g_sq = [g[0:2],g[2:]]
    G_sq = signal.convolve2d(g_sq, g_sq)
    G = [x for lst in G_sq for x in lst]
    return G

def selection(G, fitness):
    #takes in 9-element list "G" and 9-element list "fitness"
    #calculates relative fitness and returns post-selection genotype frequencies 
    Gfit = [x*y for x,y in zip(G, fitness)]
    return [i/sum(Gfit) for i in Gfit]

def meiosis(G, seg):
    #takes in 9-element list of genotype frequencies "G" and 9-element list of segregation patterns "seg"
    #each segregation pattern in "seg" is the gamete frequencies produced from each genotype
    g00 = [seg[0][i] * G[0] for i in range(4)]
    g01 = [seg[1][i] * G[1] for i in range(4)] 
    g02 = [seg[2][i] * G[2] for i in range(4)]
    g10 = [seg[3][i] * G[3] for i in range(4)]
    g11 = [seg[4][i] * G[4] for i in range(4)]
    g12 = [seg[5][i] * G[5] for i in range(4)]
    g20 = [seg[6][i] * G[6] for i in range(4)]
    g21 = [seg[7][i] * G[7] for i in range(4)]
    g22 = [seg[8][i] * G[8] for i in range(4)]
    g = [sum(x) for x in zip(g00, g01, g02, g10, g11, g12, g20, g21, g22)]
    return g

def mutation(g, mu, nu):
    #takes in a 4-element list of gamete frequencies "g" [g00, g01, g10, g11], forward mutation rate "mu", and backmutation rate "nu"
    #returns post-mutation gamete frequencies 
    g00 = g[0]*(1-mu)*(1-mu) + g[1]*(1-mu)*nu + g[2]*nu*(1-mu) + g[3]*nu*nu 
    g01 = g[0]*(1-mu)*mu + g[1]*(1-mu)*(1-nu) + g[2]*nu*mu + g[3]*nu*(1-nu)
    g10 = g[0]*mu*(1-mu) + g[1]*mu*nu + g[2]*(1-nu)*(1-mu) + g[3]*(1-nu)*nu
    g11 = g[0]*mu*mu + g[1]*mu*(1-nu) + g[2]*(1-nu)*mu + g[3]*(1-nu)*(1-nu)
    return [g00, g01, g10, g11]

def set_seg(beta, gamma):
    G00 = [1.0, 0.0, 0.0, 0.0]
    G01 = [(1/2 + beta/16), (1/2 - 5*beta/16), (3*beta/16), (beta/16)]
    G02 = [((2*beta + beta*(1-beta)*(1-gamma))/8), (1 - (6*beta + beta*(1-beta)*(1-gamma))/8), ((beta + beta*beta*(1-gamma) + beta*gamma)/8), ((2*beta + beta*(1-beta)*(1-gamma))/8)]
    G10 = [(1/2 + beta/16), (3*beta/16), (1/2 - 5*beta/16), (beta/16)]
    G11 = [(1/4 - (beta*(1-beta)*(1-gamma))/16), (1/4 + (beta*(1-beta)*(1-gamma))/16), (1/4 + (beta*(1-beta)*(1-gamma))/16), (1/4 - (beta*(1-beta)*(1-gamma))/16)]
    G12 = [(beta/16), (1/2 - 5*beta/16), (3*beta/16), (1/2 + beta/16)]
    G20 = [((2*beta + beta*(1-beta)*(1-gamma))/8), ((beta + beta*beta*(1-gamma) + beta*gamma)/8), (1 - (6*beta + beta*(1-beta)*(1-gamma))/8), ((2*beta + beta*(1-beta)*(1-gamma))/8)]
    G21 = [(beta/16), (3*beta/16), (1/2 - 5*beta/16), (1/2 + beta/16)]
    G22 = [0.0, 0.0, 0.0, 1.0]
    return [G00, G01, G02, G10, G11, G12, G20, G21, G22]

def set_fitness(s, h1, h2, h3):
    return [1, 1-h1*s, 1-h2*s, 1-h1*s, 1-h2*s, 1-h3*s, 1-h2*s, 1-h3*s, 1-s]

def process_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", type=float, nargs='?', default = 0.0, const = 0.0, help="Selection Coefficient.")
    parser.add_argument("-h1", type=float, nargs='?', default = 0.0, const = 0.0, help="Dominance coefficient to modify selection acting on G1 individuals.")
    parser.add_argument("-h2", type=float, nargs='?', default = 0.0, const = 0.0, help="Dominance coefficient to modify selection acting on G2 individuals.")
    parser.add_argument("-h3", type=float, nargs='?', default = 0.0, const = 0.0, help="Dominance coefficient to modify selection acting on G3 individuals.")
   
    parser.add_argument("-mu", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Mutation Rate.")
    parser.add_argument("-nu", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Back Mutation Rate.")

    parser.add_argument("-beta", type=restricted_float, nargs='?', default = 0.0, const=0.0, help="Average percentage of pairs of homoeologs involved in HEs per generation.")
    parser.add_argument("-gamma", type=float, nargs='?', default = 0.0, const=0.0, help="Interference/Synergy parameter of HEs.")
 
    parser.add_argument("-N", type = int, nargs='?', default = 0, const = 0, help="Number of generations to run the simulation.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-g', action='store_true')
    group.add_argument('-G', action='store_true')

    parser.add_argument("-G00", type=restricted_float, nargs='?', default = 1.0, const = 1.0, help="Percent of Individuals with Genotype of 00 in initial generation.")
    parser.add_argument("-G01", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 01 in initial generation.")
    parser.add_argument("-G02", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 02 in initial generation.")
    parser.add_argument("-G10", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 10 in initial generation.")
    parser.add_argument("-G11", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 11 in initial generation.")
    parser.add_argument("-G12", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 12 in initial generation.")
    parser.add_argument("-G20", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 20 in initial generation.")
    parser.add_argument("-G21", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 21 in initial generation.")
    parser.add_argument("-G22", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Individuals with Genotype of 22 in initial generation.")

    parser.add_argument("-g00", type=restricted_float, nargs='?', default = 1.0, const = 1.0, help="Percent of Gametes with Genotype of 00 in initial generation.")
    parser.add_argument("-g01", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Gametes with Genotype of 01 in initial generation.")
    parser.add_argument("-g10", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Gametes with Genotype of 10 in initial generation.")
    parser.add_argument("-g11", type=restricted_float, nargs='?', default = 0.0, const = 0.0, help="Percent of Gametes with Genotype of 11 in initial generation.")

    parser.add_argument("-o", nargs='?', default='output_allo_sim.txt', help="Filename to print output of each generation. If not filename is provide, then output will be written to 'output_allo_sim.txt'.")

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

def check_gamma(beta, gamma):
    gamma_min = -(0.5 - abs(0.5 - beta)) / (0.5 + abs(0.5 - beta))

    if gamma > 1.0 or gamma < gamma_min:
        raise Exception("Gamma value invalid. For your value of beta, gamma must a value between 1 and ", gamma_min)
    print(str(gamma_min))
    return


def main():
    inargs = process_args()

    header = 'Generation,s,h1,h2,h3,mu,nu,beta,gamma,g00,g01,g10,g11,G00,G01,G02,G10,G11,G12,G20,G21,G22'

    if inargs.G:
        gams = [np.nan, np.nan, np.nan, np.nan]
        genos = [inargs.G00, inargs.G01, inargs.G02, inargs.G10, inargs.G11, inargs.G12, inargs.G20, inargs.G21, inargs.G22]
        if abs(sum(genos) - 1.0) > 1e-6:
            raise Exception("Initial Genotype Frequencies must sum to 1.0, not ", sum(genos)) 

    if inargs.g:
        gams = [inargs.g00, inargs.g01, inargs.g10, inargs.g11]
        if abs(sum(gams) - 1.0) > 1e-6:
            raise Exception("Initial Gamete Frequencies must sum to 1.0, not ", sum(gams))
        genos = random_mating(gams)

    check_gamma(inargs.beta, inargs.gamma)

    fitness = set_fitness(inargs.s, inargs.h1, inargs.h2, inargs.h3)

    segs = set_seg(inargs.beta, inargs.gamma)

    gen_stats = np.array([[x for y in [[0, inargs.s, inargs.h1, inargs.h2, inargs.h3, inargs.mu, inargs.nu, inargs.beta, inargs.gamma], gams, genos] for x in y]])
    gams = generation(genos, segs, fitness, inargs.mu, inargs.nu)
    genos = random_mating(gams)

    for i in range(1, inargs.N):
        add_stack = [x for y in [[i, inargs.s, inargs.h1, inargs.h2, inargs.h3, inargs.mu, inargs.nu, inargs.beta, inargs.gamma], gams, genos] for x in y]
        gen_stats = np.vstack([gen_stats, add_stack])
        gams = generation(genos, segs, fitness, inargs.mu, inargs.nu)
        genos = random_mating(gams)

    np.savetxt(inargs.o, gen_stats, delimiter=',', header = header, fmt = ['%g', '%.3e', '%.3e', '%.3e', '%.3e', '%.3e', '%.3e', '%.3e', '%.3e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e', '%.18e'])

main()
