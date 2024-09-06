#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 23:52:24 2021

@author: Prog_ReS.S
"""
"""
This function is to use recurrsion to generate a list of sequence with 'd' number of variations of the seqence. 

pattern_O is the sequence you want variations of (O standing for Original. An unaltered comparesent outside of the function is need for comparing)

d is the number of difference you accept. 
"""
import time
pattern_O = 'AGGGCTCTTAAG'
d = 2
dict_of_variation = {}
base_pair = 'A','T','C','G'
def Neighbors (pattern, d):
    start = time.process_time()
    working_pattern = pattern
    for i in range(len(pattern)):
        print ("finding prefix and suffix")
        prefix = pattern[0:i]
        suffix = pattern[i+1:len(pattern)]
        for b in base_pair:
            working_pattern = prefix+b+suffix
            print ("Going though base: "+ str(b) + " At location: " + str(i))
            if HammingDistance (pattern_O, working_pattern) <= d:
                print ("hamming distence accepted")
                if working_pattern in dict_of_variation:
                    pass
                else:
                    dict_of_variation[working_pattern]=+1
                    print ("seqence added, Recalling function")
                    Neighbors(working_pattern, d)
    end = time.process_time()
    print("Time:", end - start)
    return (dict_of_variation)

def HammingDistance(p, q):

   count = 0

   for i in range(len(p)):
       x = p[i]
       y = q[i]
       if x != y:
           count = count + 1
   return count          

print (*Neighbors (pattern_O, d))