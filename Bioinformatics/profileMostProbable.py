#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 20:49:40 2021

@author: Prog_ReS.S
"""
def compute_probability(kmer,profile_mat):
  prob=1
  for i in range(0,len(kmer)):
    if kmer[i]=='A':
      prob=prob*profile_mat[0][i]
    if kmer[i]=='C':
      prob=prob*profile_mat[1][i]
    if kmer[i]=='G':
      prob=prob*profile_mat[2][i]
    if kmer[i]=='T':
      prob=prob*profile_mat[3][i]
  return prob


def profile(dna_seq,k,profile_matrix):
  max_probab=0
  for i in range(len(dna_seq)-k+1):
    kmer=dna_seq[i:i+k]
    pro=compute_probability(kmer,profile_matrix)
    if pro>max_probab:
      max_probab=pro
      sol=kmer
  return sol

filename = 'profile.txt'
with open(filename, "r") as dataset:
    data = []
    for line in dataset:
        data.append(line.strip())
    Text = data[0]
    k = int(data[1])
    raw_profile = data[2:]
    bases = ['A', 'C', 'G', 'T']
    prof = [list(map(float, raw_profile[i].split())) for i in range(len(raw_profile))]
    prof_dict = dict(zip(bases, prof))
print(profile(Text, k, prof))