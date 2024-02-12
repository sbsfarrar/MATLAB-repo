#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 19:26:53 2021

@author: Prog_ReS.S
"""
dna = [
       "TCTACATCAAAGGGAGATATGCACA",
       "TGTACACCATCCCCACCTTGTGGGA",
       "CGGCACCAGTGCCTTTAAGTCACCT",
       "CGGAAATTGGCAGAGCTAATGAGCT",
       "TAAACTGACCTGTTCAGGGAGTGTC",
       "ATCTTGGGCATGCTGCGGTGAAGAA"
       ]
k = 5
d = 2

def HammingDistance(p, q):
    count=0
    for i in range (0,len(p)):
        if p[i]!=q[i]:
            count+=1

    return count

def Neighbors(Pattern, d):
    if d==0:
        return [Pattern]
    elif len(Pattern)==1:
        return ['A','C','G','T']
    Neighborhood=[]
    Suffix_Pattern=Pattern[1:]
    FirstSymbol_Pattern=Pattern[0]

    SuffixNeighbors = Neighbors(Suffix_Pattern, d)
    for Text in SuffixNeighbors:
        if HammingDistance(Suffix_Pattern, Text) < d:
            for nucleotide in ['A','C','G','T']:
                Neighborhood.append(nucleotide+Text)
        else:
            Neighborhood.append(FirstSymbol_Pattern+Text)
    return Neighborhood

# Write your MotifEnumeration() function here along with any subroutines you need.
# This function should return a list of strings.
def MotifEnumeration(dna, k, d):
    patterns=[]
    for i in range (0,len(dna[0])-k+1):
        neighbors=Neighbors(dna[0][i:i+k],d)
        for j in neighbors:
            count=0
            for l in dna:
                for i in range(0,len(l)-k+1):
                    if HammingDistance(j, l[i:i+k])<=d:
                        count+=1
                        break
            if count==len(dna):
                patterns.append(j)
    Patterns = [] 
    [Patterns.append(x) for x in patterns if x not in Patterns] 
    Result = ""
    for item in Patterns:
        Result = Result + " " + str(item)
    return Result
print(MotifEnumeration(dna,k,d))