#! /usr/bin/env python3

import sys
import json
import argparse
'''
This program determines how the position of a nucleotide within a context window changes its prediction.
'''

def parseArgs():
    parser = argparse.ArgumentParser(description='Get the prediction values at a single position within the sequence. Output values represent the prediction within the context window (i.e., is the nucleotide at the beginning of the prediction sequence or the end?). Requires output from parseWindow.py')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input file from parseWindow.py output.')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output text file')
    parser.add_argument('-n', '--number_of_tokens_per_seq', type=int, default=4096,required=False,help='Number of tokens per sequence')
    parser.add_argument('-t', '--token_size', type=int, default=6,required=False,help='Token size. Most models currently use a token size of 6. Do not change unless a model uses a different token length')
    parser.add_argument('-p', '--position', type=int, default=850,required=False,help='Position within the sequence to query')
    args = parser.parse_args()
    return args

def getOnePosition(args):
    window= (args.number_of_tokens_per_seq *args.token_size) -1 #4096*6 minus 1 since one nucleotide overlaps the window from the gene.
    num=0
    inArrays=False
    intron=False
    exon=False
    exon_list=[]
    needed_pos=args.position #The position being queried
    intron_list=[]
    with open(args.input) as inputF:
        for line in inputF:
            line=line.strip()
            if line =="":
                inArrays=False
                continue
            elif line == "Intron":
                num=0
                exon=False
                intron=True
                inArrays=True
                continue
            elif line == "Exon":
                num=0
                exon=True
                intron=False
                inArrays=True
                continue
            if inArrays:
                num +=1
                yes=False
                if exon:
                    exon_list.append(json.loads(line.strip())[window-num+needed_pos])
                elif intron:
                    intron_list.append(json.loads(line.strip())[window-num+needed_pos])
    
    with open(args.output,'w') as output:
        output.write("exon=" + str(exon_list)+"\n")
        output.write("intron=" + str(intron_list)+"\n")

def main():
    args=parseArgs()
    getOnePosition(args)

if __name__ == '__main__':
    main()

