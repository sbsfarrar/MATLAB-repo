#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 19:37:58 2021

@author: Prog_ReS.S
"""
Dna = [
        "ATTTGGC",
        "TGCCTTA",
        "CGGTATC",
        "GAAAATT"
        ]
k = 3
d = 1
def Suffix(Pattern):
    if len(Pattern)==1:
        return ""
    suffix = Pattern[1:]
    return suffix

def HammingDistance(p, q):
    n = len (p) #since both are of equal distances
    hamming_distance = 0
    for i in range(0,n):
        if p[i] != q[i]:
            hamming_distance += 1
    return hamming_distance

def FirstSymbol(Pattern):
    x = Pattern[0]
    return x

def Neighbours(Pattern,d):
    if d==0:
        return {Pattern}
    if len(Pattern)==1:
        x = ['A','C','T','G']
        return x
    Neighbourhood = set()
    SuffixNeighbours = Neighbours(Suffix(Pattern),d)
    for Text in SuffixNeighbours:
        if HammingDistance(Suffix(Pattern),Text) < d:
            for y in "ATCG":
                string = y + Text
                Neighbourhood.add(string)
        else:
            m = FirstSymbol(Pattern) + Text
            Neighbourhood.add(m)
    Neighbourhood = list(Neighbourhood)
    return Neighbourhood

def MotifEnumeration(Dna, k, d):
    patterns = set()
    n = len(Dna)
    AllKmers = []
    for j in range(n):
        a = Dna[j]
        kmers = []
        n_1 = len(a)
        for i in range(n_1 - k + 1):
            kmers.append(a[i:i+k])
        neigh = []
        for i in range(len(kmers)):
            l = Neighbours(kmers[i],d)
            for val in l:
                neigh.append(val)
        AllKmers.append(neigh)   
    x1 = set(AllKmers[0])
    x2 = set(AllKmers[1])
    patterns = x1 & x2
    for y in range(2,len(AllKmers)):
        patterns = patterns & set(AllKmers[y])
    patterns = list(patterns)
    string = ""
    for i in patterns:
        string = string + i + " "
    return string

# if __name__ == "__main__":
#     k = int(input())
#     d = int(input())
#     Dna = []
#     for i in range(0, 6): 
#         ele = input()
#         Dna.append(ele)
#     p = MotifEnumeration(Dna, k, d)
#     print(' '.join(p))#

#print(MotifEnumeration(Dna, k, d))