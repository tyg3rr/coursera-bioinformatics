# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 09:44:58 2021

@author: ljens
"""
'''
sample_input = input('input')
k = int(input('number'))
'''

def freqwords(text, k):   #returns most frequent k-mers
    kmers_counts = {}
    for i in range(0, len(text)-k+1):
        kmer = text[i:(i+k)]
        if kmer in kmers_counts:
            kmers_counts[kmer] += 1
        else:
            kmers_counts[kmer] = 1
    #if text[(len(text)-k):len(text)] in kmers_counts:
        #kmers_counts[text[(len(text)-k):len(text)]]
    max_freq = 0
    for key in kmers_counts:
        if kmers_counts[key] > max_freq:
            max_freq = kmers_counts[key]
    most_freq = []
    for key in kmers_counts:
        if kmers_counts[key] == max_freq:
            most_freq.append(key)
    return most_freq

print(freqwords('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT', 3))

'''

def number_of_kmers(window, k):
    kmers_counts = {}
    for i in range(0, len(window)-k+1):
        kmer = window[i:(i+k)]
        if kmer in kmers_counts:
            kmers_counts[kmer] += 1
        else:
            kmers_counts[kmer] = 1
    return kmers_counts
    

k refers to k-mer
L is the length of the window
t is frequency of appearances in the window
'''

def clump_finding(genome, k, L, t):
    kmers_counts = {}
    for i in range(0, int(L)-int(k)+1):
        kmer = genome[i:(i+int(k))]
        if kmer in kmers_counts:
            kmers_counts[kmer] += 1
        else:
            kmers_counts[kmer] = 1
    clumps = []
    print(kmers_counts)
    for kmer in kmers_counts:
        if kmers_counts[kmer] >= int(t):
            clumps.append(kmer)
    return clumps
    
'''
genome = input('genome ')
k = input('k ')
L = input('l ')
t = input('t ')


genome = 'ACGTACGT'
k = '1'
L = '5'
t = '2'
print(clump_finding(genome, k, L, t))
    
'''
    
    
    
    
    