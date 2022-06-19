# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 12:39:48 2021

@author: ljens
"""


def freqwords(sub_genome, k):   #returns kmer:frequency dictionary
    kmers_counts = {}
    for i in range(0, len(sub_genome)-k+1):
        kmer = sub_genome[i:(i+k)]
        if kmer in kmers_counts:
            kmers_counts[kmer] += 1
        else:
            kmers_counts[kmer] = 1
    return kmers_counts




def clump(genome, k, L, t):
    clumps = []
    for i in range(0, len(genome)-L+1):
        sub_genome = genome[i:(i+L)]
        counts = freqwords(sub_genome, k)
        for key in counts:
            if counts[key] >= t:
                clumps.append(key)
    clumpsu = []
    [clumpsu.append(x) for x in clumps if x not in clumpsu]
    return clumpsu

genome = input('genome')
k = int(input('k'))
L = int(input('L'))
t = int(input('t'))

print(clump(genome, k, L, t))
    
        
        
        


#print(freqwords('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', 5))
