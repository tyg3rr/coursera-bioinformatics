# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 11:51:23 2021

@author: ljens
"""

complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

def reverse_complement(sequence):
    r = []
    for i in sequence:
        r.append(complements[i])
    r.reverse()
    return ''.join(r)


pattern = input('Input pattern: ')
#genome = input('Input genome: ')


def positions(pattern, genome):
    pattern_complement = reverse_complement(pattern)
    positions = []
    for i in range(0, len(genome)-len(pattern)+1):
        if genome[i:(i+len(pattern))] == pattern: #or genome[i:(i+len(pattern))] == pattern_complement:
            positions.append(str(i))
    return ' '.join(positions)

#print(positions(pattern, genome))
    
print(reverse_complement(pattern))
