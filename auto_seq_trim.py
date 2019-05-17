#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:44:24 2018
@author: skylar
"""

import argparse
from Bio import AlignIO               # Read fasta files
from Bio import SeqIO

def get_raw(aln, ret):
    """Passes the gap percentages and Amino Acid counts back by reference."""
    rows = len(aln)
    col = ret["col"]
    gap = [0]*col
    num_gaps = 0
    fault = [False]*col
    faults = []
    for i in range(ret["start"],ret["end"]):
        for j in range(rows):
            char = (aln[j].seq)[i]
            if char in ('.', '-'):
                num_gaps += 1
        gap[i] = num_gaps/rows
        if gap[i] > ret["gap_score"]:
            fault[i] = True
        num_gaps = 0
        faults = []
        j = 0
    for i in range(col):
        if faults != [] and i < faults[len(faults)-1][1]:
            i = faults[len(faults)-1][1]+2
        fault_len = 0
        if i < len(fault) -1 and fault[i]:
            for j in range(i, col):
                if not fault[j]:
                    break
                fault_len+= 1
            if fault_len >= ret["trim_score"]:
                temp = []
                temp.append(i)
                temp.append(i+fault_len-1)
                faults.append(temp)
    for i in range(len(faults)):
        for j in aln:
            aa_ct = 0
            for k in range(faults[i][0],faults[i][1]+1):
                if j.seq[k] != '-' and j.seq[k] != '.':
                    aa_ct += 1
            if aa_ct >= ret["blame"]*(faults[i][1]-faults[i][0]):
                if not j.name in ret["out"]:
                    ret["out"][j.name] = ""
                ret["out"][j.name] += str(faults[i][0])+","+str(faults[i][1])+","
    

class Groups():
    """Handles data storage."""
    grp_dict = dict()
    dat_dict = dict()
    aln = None
    def __init__(self, alignments, gaps, trim, blame, start, end):
        for i in alignments:
            self.grp_dict[i.name] = i
        self.aln = alignments
        self.dat_dict["gap_score"] = gaps
        self.dat_dict["trim_score"] = trim
        self.dat_dict["blame"] = blame
        self.dat_dict["end"] = end
        self.dat_dict["start"] = start
        self.dat_dict["out"] = dict()
        self.dat_dict["col"] = len(alignments[0].seq)
        get_raw(self.aln, self.dat_dict)
        for i in self.dat_dict["out"]:
            self.grp_dict.pop(i)
    def get_seq(self):
        l = list(self.grp_dict.values())
        SeqIO.write(l, "out.fasta", "fasta")
        

def main():
    """handles IO and calls necessary functions from the group class"""
    parser = argparse.ArgumentParser(
        description='Calculate group entropy')
    parser.add_argument('Alignment', type=str,
                        help='Alignment file in fasta format')
    parser.add_argument('gapScore', type=float,
                        help='Percent of gaps at a position to consider it a gap position')
    parser.add_argument('numConGap', type=float,
                        help='Number of consecutive gaps to run trimming')
    parser.add_argument('blame', type=float,
                        help='Proportion of gaps created for blame to fall')
    parser.add_argument('start', type=int,
                        help='What length from the ends to consider part of the front end terminal')
    parser.add_argument('end', type=int,
                        help='What length from the ends to consider part of the back end terminal')
    
    gaps = parser.parse_args().gapScore
    trim = parser.parse_args().numConGap
    blame = parser.parse_args().blame
    end = parser.parse_args().end
    start = parser.parse_args().start
    aln = AlignIO.read(parser.parse_args().Alignment, "fasta")
    grp = Groups(aln, gaps, trim, blame, start, end)
    grp.get_seq()
    out = open("trim.csv", 'w')
    
    for i in grp.dat_dict["out"]:
        out.write(i+","+grp.dat_dict["out"][i]+"\n")
if __name__ == '__main__':
    main()
