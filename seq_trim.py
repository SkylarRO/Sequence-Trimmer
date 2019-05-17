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
    for j in aln:
        aa_ct = 0
        for k in range(ret["start"],ret["end"]+1):
            if j.seq[k] != '-' and j.seq[k] != '.':
                aa_ct += 1
        if aa_ct >= ret["blame"]*(ret["end"]-ret["start"]):
            ret["out"][j.name] = ""
    

class Groups():
    """Handles data storage."""
    grp_dict = dict()
    dat_dict = dict()
    aln = None
    def __init__(self, alignments, blame, start, end):
        for i in alignments:
            self.grp_dict[i.name] = i
        self.aln = alignments
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
    parser.add_argument('alignment', type=str,
                        help='Alignment file in fasta format')
    parser.add_argument('blame', type=float,
                        help='Proportion of gaps created for blame to fall')
    parser.add_argument('s_pos', type=int,
                        help='Start position')
    parser.add_argument('e_pos', type=int,
                        help='End position')
    blame = parser.parse_args().blame
    end = parser.parse_args().e_pos
    start = parser.parse_args().s_pos
    aln = AlignIO.read(parser.parse_args().alignment, "fasta")
    grp = Groups(aln, blame, start, end)
    grp.get_seq()
    out = open("trim.csv", 'w')
    
    for i in grp.dat_dict["out"]:
        out.write(i+"\n")
if __name__ == '__main__':
    main()
