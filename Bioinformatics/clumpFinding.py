#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 22:21:58 2021

@author: Prog_ReS.S
"""
import itertools
from collections import Counter


def find_frequent(string, k,t):
    words = []
    frequent = []

    for i in range(len(string)):
        word = "".join(string[i: i + k])

        if len(word) == k:
            words.append(word)

    return Counter(words).most_common()

def clump_finding_problem(string, k, L, t):
    string = string[3764856:3764856+500]
    words = []
    for i in range(len(string)):
        strings1 = string[i:i + L]
        if len(strings1) == L:
            words.append(find_frequent(strings1, k, t))

    pattern = list(itertools.chain(*words))
    print(*set([x[0] for x in pattern if x[1] >= t]))

#genome = open('/Users/Prog_ReS.S/Downloads/Salmonella_enterica.txt').readlines()
#remove the header
#genome = genome[1:-1]
#data = "".join(open('rosalind_ba1e.txt')).split()
#clump_finding_problem(data[0], int(data[1]), int(data[2]), int(data[3]))
