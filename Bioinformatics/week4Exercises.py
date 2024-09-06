#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 18:39:38 2021

@author: Prog_ReS.S
"""
##############################################################################
#Compute the probability that ten randomly selected 15-mers from the ten 
#600-nucleotide long strings in the Subtle Motif Problem capture at least 
#one implanted 15-mer. (Allowable error: 0.000001)
##############################################################################
kmer = 15
p1 = (600-kmer)/(600-kmer+1)
pC = pow(p1,10)
pAnswer = 1-pC
print(pAnswer)

#############################################################################
# Compute the probability that ten randomly selected 15-mers from ten 
# 600-nucleotide long strings (as in the Subtle Motif Problem) capture at 
# least two implanted 15-mers. (Allowable error: 0.000001)
#############################################################################
kmer = 15
p1 = (600-kmer)/(600-kmer+1)
p2 = 1-p1

from itertools import *
counter = 0
for seq in combinations(range(10),2):
    counter += 1
counter

pAns = pow(p2,2)*pow(p1,8)*counter
print(pAns)