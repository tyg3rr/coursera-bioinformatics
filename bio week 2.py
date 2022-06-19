# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 14:01:58 2021

@author: ljens
"""

import itertools
import os
import math

def skew(genome):    #returns list of skew values
    count = 0
    skew_list = [0]
    for n in genome:
        if n == 'G':
            count += 1
            skew_list.append(count)
        if n == 'C':
            count -= 1
            skew_list.append(count)
        if n == 'A' or n == 'T':
            count += 0
            skew_list.append(count)
    return skew_list

def min_skew(skew_list):    #returns list of min skew locations in genome
    min_value = min(skew_list)
    min_locations = []
    for i in range(0, len(skew_list)):
        if skew_list[i] == min_value:
            min_locations.append(i)

    return ' '.join(str(x) for x in min_locations)

def hamming(genome1, genome2):    #returns hamming integer for two genes
    count = 0
    for i in range(0, len(genome1)-1):
        if not genome1[i] == genome2[i]:
            count += 1
    return count

def pattern_match(genome, subgenome, d):    #returns locations of matches along genome
    matches = []
    for i in range(0, len(genome)-len(subgenome)+1):
        ham = hamming(genome[i:(i+len(subgenome))], subgenome)
        if ham <= int(d):
            matches.append(i)
    #return ' '.join(str(x) for x in matches)
    return matches

def kmer_list(genome, k):    #returns list of every possible kmer
    kmers = []
    for i in range(0, len(genome)-int(k)+1):
        kmers.append(genome[i:(i+k)])
    return kmers


def freq_mis(genome, k, d):    #returns list of most freq approx matches
    frequencies = {}
    answer = []
    invisible = []
    for thing in itertools.product('ATCG', repeat = int(k)):
        invisible.append(list(thing))
    for item in invisible:
        m = ''.join(item)
        if m in frequencies:
            pass
        else:
            item_matches = len(pattern_match(genome, m , d))
            frequencies[m] = item_matches
    maxi = 0
    for value in frequencies.values():
        if value > maxi:
            maxi = value
    for key, value in frequencies.items():
        if value == maxi:
            answer.append(key)
    return answer



def reverse_complement(sequence):    #returns list of most freq approx matches and their reverse complements
    complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    r = []
    for i in sequence:
        r.append(complements[i])
    r.reverse()
    return r

def rc_freq_mis(genome, k, d):
    frequencies = {}
    answer = []
    invisible = []
    for thing in itertools.product('ATCG', repeat = int(k)):
        invisible.append(list(thing))
    for item in invisible:
        m = ''.join(item)
        if m in frequencies:
            pass
        else:
            item_matches = len(pattern_match(genome, m , d))
            rc_item_matches = len(pattern_match(reverse_complement(genome), m, d))
            frequencies[m] = item_matches + rc_item_matches
    print(frequencies)
    '''
    for key in frequencies:
        t = reverse_complement(key)
        if t in frequencies:
            frequencies[key] = frequencies[t] + frequencies[key]
            #frequencies[t] = 0
    print(frequencies)
    '''
    maxi = 0
    for value in frequencies.values():
        if value > maxi:
            maxi = value
    for key, value in frequencies.items():
        if value == maxi:
            answer.append(key)
            #answer.append(reverse_complement(key))

    return ' '.join(x for x in answer)



def immediate_neighbors(gene):    #returns list of neighbors with hamming = 1
    gene = list(gene)
    nucleotides = ['A', 'T', 'C', 'G']
    immediate_neighbors = []
    for i in range(0, len(gene)):
        for n in nucleotides:
            immediate_neighbors.append(gene[0:i] + [n] + gene[(i+1):])
    for element in immediate_neighbors:
        if element == gene:
            immediate_neighbors.remove(element)
    immediate_neighbors.append(gene)
    return immediate_neighbors


def neighborhood(gene, d):    #returns list of neighbors with hamming = d
    immediate_n = immediate_neighbors(gene)
    neighborhood = []
    for element in immediate_n:
        neighborhood.append(element)
    if d == 1:
        return neighborhood
    else:
        i = 1
        k = len(immediate_n)
        while i <= d:
            for y in range(0, k):
                temp = immediate_neighbors(immediate_n[y])
                for thingie in temp:
                    immediate_n.append(thingie)
                i += 1
        for element in immediate_n:
            if element not in neighborhood:
                neighborhood.append(element)
        return neighborhood


def freqwords(sub_genome, k):   #returns kmer:frequency dictionary
    kmers_counts = {}
    for i in range(0, len(sub_genome)-k+1):
        kmer = sub_genome[i:(i+k)]
        if kmer in kmers_counts:
            kmers_counts[kmer] += 1
        else:
            kmers_counts[kmer] = 1
    return kmers_counts




def clump(genome, k, L, t):    #returns genes that appear t or more times
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




def brute_motif(A, k, d):
    motifs = set()
    for i in range(0, len(A)):
        for kmer in kmer_list(A[i], k):
            motifs.add(kmer)
            for u in neighborhood(kmer, d):
                motifs.add(''.join(u))
    g = []
    for motif in motifs:
        boo = True
        for i in range(0, len(A)):
            if pattern_match(A[i], motif, d):
                pass
            else:
                boo = False
        if boo == False:
            g.append(motif)
    for thingie in g:
        motifs.discard(thingie)
    return ' '.join(motifs)

def all_pos(k):
    all_pos = itertools.product('ACGT', repeat = int(k))

    return [''.join(x) for x in all_pos]

def median_string(k, a):
    motifs = dict.fromkeys(all_pos(k))
    d = 10000000000
    med = str()
    for kmer in motifs:
        total_sum = 0
        for string in a:
            baby_sum = []
            for i in range(0, len(string)-k+1):
                baby_sum.append(hamming(kmer, string[i:i+k-1]))
            total_sum += min(baby_sum)
        motifs[kmer] = total_sum
    print(motifs)
    for (key, value) in motifs.items():
        if value < d:
            d = value
            med = key
    return med

'''
file_path = input('Enter Filepath: ')
file_input1 = open(file_path, 'r')
all_lines = file_input1.readlines()
gene = all_lines[0].strip()
k = int(all_lines[1].strip())
profile = []
for j in range(2, 6):
    profile.append([float(x) for x in all_lines[j].strip().split()])


file_input1.close()

'''



def most_probable_kmer(k, profile, gene):
    dick = profile
    print(dick)
    k = int(k)
    probs = {}
    for i in range(0, len(gene)+1-k):
        kmer = gene[i:i+k]
        baby_probs = []
        for j in range(0, k-1):
            baby_probs.append(dick[(kmer[j])][j])
        probs[kmer] = math.prod(baby_probs)
    prob_flag = 0
    most_probable = None
    for key, value in probs.items():
        if value > prob_flag:
            prob_flag = value
            most_probable = key
    return most_probable

dna = ['GGCGTTCAGGCA',
       'AAGAATCAGTCA',
       'CAAGGAGTTCGC',
       'CACGTCAATCAC',
       'CAATAATATTCG']
motif = 'GGC'
profile = {'A':[], 'C':[], 'G':[], 'T':[]}
substring_kmers = []
substring = dna[0]
for i in range(0, len(dna[0])-2):
    substring_kmers.append(substring[i:i+4])
for kmer in substring_kmers:
    for j in range(0, 4):





def profile_generator(list_of_lists, length_of_element,k):
    profile = {'A':[], 'C':[], 'G':[], 'T':[]}
    for i in range(0, length_of_element):
        for (key, value) in profile.items():
            count_key = 0
            for item in list_of_lists:
                if item[i] == key:
                    count_key += 1
            profile[key].append(count_key/len(list_of_lists))
    return profile


def greedy_motif(dna, k, t):
    arbritrary_kmers = [dna[0][0:k]]
    profile = profile_generator(dna[0], 1, k)
    for i in range(1, t):
        kmer = [most_probable_kmer(k, profile, [dna[i]])]
        arbritrary_kmers.append(kmer)
        profile = profile_generator(arbritrary_kmers, len(arbritrary_kmers), k)

    return arbritrary_kmers




#print(greedy_motif(dna, 3, 5))











