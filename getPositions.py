#! /usr/bin/env python3

import sys
import json
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(description='Get the prediction values across the entire gene sequence for the first, middle, and last nucleotide. Requires output from parseWindow.py')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input file from parseWindow.py output.')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output text file')
    parser.add_argument('-n', '--number_of_tokens_per_seq', type=int, default=4096,required=False,help='Number of tokens per sequence')
    parser.add_argument('-t', '--token_size', type=int, default=6,required=False,help='Token size. Most models currently use a token size of 6. Do not change unless a model uses a different token length')
    #parser.add_argument('-g', '--gene_length', type=int, default=3612,required=False,help='The number of nucleotides in the gene sequence. Default is size of APOE') #This argument was once used as a quality control. Commented lines below pertain to this argument
    args = parser.parse_args()
    return args

def writeToFile(outputFile,exon_lists,intron_lists): 
    with open(outputFile,'w') as output:
        output.write("{\n\"Exon\":\n")
        output.write("{\"first\": " +str(exon_lists[0]) +",\n")
        output.write("\"middle\": " +str(exon_lists[1]) +",\n")
        output.write("\"last\": "   +str(exon_lists[2]) +"\n},\n")
        output.write("\"Intron\":\n")
        output.write("{\"first\": " +str(intron_lists[0]) +",\n")
        output.write("\"middle\": " +str(intron_lists[1]) +",\n")
        output.write("\"last\": "   +str(intron_lists[2]) +"\n}\n}\n")

def getPositions(args):
    inputF=open(args.input)
    window=int(inputF.readline().split(" ")[0])
    #gene_length=args.gene_length
    num=0
    inArrays=False
    intron=False
    exon=False
    exon_lists=[]
    intron_lists=[]
    #expected_length=args.gene_length + (args.number_of_tokens_per_seq *args.token_size)-1 #Minus 1 since there will always be one overlapping nucleotide
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
        #Get first, middle, last
        if inArrays:
            num +=1
            valuesNeeded=False
            if num ==1:
                line=json.loads(line.strip())
                '''
                if len(line) != expected_length:
                    line_length=len(line)
                    raise ValueError(f"Line has unexpected length of {line_length}. Expected {expected_length}.") 
                '''
                line=line[window-1:] #-1 since 1 nt overlaps
                valuesNeeded=True
            elif num ==(window/2):
                line=json.loads(line.strip())
                '''
                if len(line) != expected_length:
                    line_length=len(line)
                    raise ValueError(f"Line has unexpected length of {line_length}. Expected {expected_length}.") 
                '''
                line=line[int((window/2)-1):(-1*int((window/2)))]
                valuesNeeded=True
            elif num ==window:
                line=json.loads(line.strip())
                '''
                if len(line) != expected_length:
                    line_length=len(line)
                    raise ValueError(f"Line has unexpected length of {line_length}. Expected {expected_length}.") 
                '''
                line=line[0:((-1*window)+1)] 
                valuesNeeded=True
            if valuesNeeded:
                if exon:
                    exon_lists.append(line)
                elif intron:
                    intron_lists.append(line)
    writeToFile(args.output,exon_lists,intron_lists)

def main():
    args=parseArgs()
    getPositions(args)

if __name__ == '__main__':
    main()
