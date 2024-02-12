#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 22:48:08 2021

@author: Prog_ReS.S
"""
#data = open('/Users/Prog_ReS.S/Downloads/Salmonella_enterica.txt').readlines()
#remove the header
#data = data[1:-1]
#print minimumSkew("".join(map(str.strip,data)))

genome = open('/Users/Prog_ReS.S/Downloads/Salmonella_enterica.txt').readlines()
#remove the header
genome = genome[1:-1]

def Skew(genome):
    Skew = {}

    Skew[0] = 0
    for i in range(1, len(genome)+1):
        if genome[i - 1] == "G":
            Skew[i] = Skew[i - 1] + 1
        elif genome[i - 1] == "C":
            Skew[i] = Skew[i - 1] - 1
        else:
            Skew[i] = Skew[i-1]
    return Skew


def MinimumSkew(genome):

    positions = [] # output variable
    s = Skew(genome)
    m = min(s.values())
    for (k,v) in s.items():
        if v == m:
            positions.append(k)
    return positions
#print(MinimumSkew("".join(map(str.strip,genome))))