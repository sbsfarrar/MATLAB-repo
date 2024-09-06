#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 20:14:57 2021

@author: Prog_ReS.S
"""
def hamming(p1, p2):
    d = 0
    for i in range(len(p1)):
        if p1[i] != p2[i]:
            d += 1
    return(d)

def MedianString(infile):

    # open infile
    with open(infile, 'r') as file:
        k = int(file.readline())
        Dna = file.readlines()
        
        # iterate through each Dna string
        kmerList = []
        for line in Dna:
            String = line.strip('\n')
            for i in range(len(String) - k+1):
                pattern = String[i:i+k]
                if pattern not in kmerList:
                    kmerList.append(pattern)
            
            # pattern that minimizes hamming distance
            distance = float('inf')
            for kmer in kmerList:
                for i in range(len(String) - k+1):
                    if distance > hamming(kmer, String[i:i+k]):
                        distance = hamming(kmer, String[i:i+k])
                        Median = kmer
    return(Median)


infile = 'dataset_158_9.txt'

print(MedianString(infile))