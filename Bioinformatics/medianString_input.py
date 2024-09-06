#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 20:21:11 2021

@author: Prog_ReS.S
"""
k = 7
DNA = [
       "CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
       "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
       "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"
       ]


def Neighbors(Pattern,d):
    if d == 0:
        return Pattern
    if len(Pattern) == 1:
        return 'ACGT'
    Neighborhood = []
    SuffixNeighbors = Neighbors(Pattern[1:],d)
    for text in SuffixNeighbors:
        if HammingDistance(Pattern[1:],text) < d:
            for index in 'ACGT':
                Neighborhood.append(index+text)
        else:
                Neighborhood.append(Pattern[0]+text)
    return Neighborhood

def HammingDistance(Pattern,text):
    count = 0
    for i in range(len(Pattern)):
        if Pattern[i] != text[i]:
            count +=1
    return count

def MedianString(DNA,k):
    kmerlist = []
    for string in DNA:
        for i in range(0,len(string)-k+1):
            pattern = string[i:i+k]
            if pattern not in kmerlist:
                kmerlist.append(pattern)
    for kmer in kmerlist:
        score = 0
        for string in DNA:
            min_distance=1000
            for i in range(0,len(string)-k+1):
                pattern = string[i:i+k]
                distance = HammingDistance(kmer,pattern)
                if min_distance >= distance:
                    min_distance = distance
            score += min_distance
        if kmer == kmerlist[0]:
            bestscore = score
            median = kmer
        if bestscore >= score:
            bestscore = score
            median = kmer
    return (median,bestscore)



print(MedianString(DNA, k))