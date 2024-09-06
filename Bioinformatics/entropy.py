#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 20:09:30 2021

@author: Prog_ReS.S
"""
def Entropy(Motifs):
    Profile = {}
    # check that all motifs are the same length
    ListLength = len(Motifs)
    L1 = len(Motifs[1]) # length of the first motif
    print('There are {} motifs of length {}'.format(ListLength, L1))
    for i in range(len(Motifs)):
        if len(Motifs[i]) != L1:
            ShortMotif = Motifs[i]
            ShortMotifLen = len(ShortMotif)
            print('Oops, Motif {} is {} nucleotides instead of {}!'.format(ShortMotif, ShortMotifLen, L1))
            break
        
    # fill all positions with frequency of 0
    for nucleotide in 'ACGT':
        values = [0] * L1
        Profile[nucleotide] = values
        
    # iterate through each position in the motif matrix, counting nucleotide frequencies
    TotalEntropy = 0
    for key, values in Profile.items():
        for Motif in Motifs:
            for i in range(len(Motif)):
                if Motif[i] == key:
                    Profile[key][i] += 1
        
        # convert nucleotide frequencies to probabilities
        for i in range(len(values)):
            Profile[key][i] = Profile[key][i] / float(ListLength)
        
        # calculate total entropy (Sum of (Prob_value * log2 Prob_n))
        import math
        for value in values:
            if value > 0:
                TotalEntropy += abs(value * math.log(value, 2))
            else: continue
            
    return(TotalEntropy)


Motifs1 = [
"TCGGGGGTTTTT",
"CCGGTGACTTAC",
"ACGGGGATTTTC",
"TTGGGGACTTTT",
"AAGGGGACTTCC",
"TTGGGGACTTCC",
"TCGGGGATTCAT",
"TCGGGGATTCCT",
"TAGGGGAACTAC",
"TCGGGTATAACC"
]

print(Entropy(Motifs1))