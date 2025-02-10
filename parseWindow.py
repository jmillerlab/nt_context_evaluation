#! /usr/bin/env python3
'''
The purpose of this program is to read through the output from run_nt.py. run_nt.py takes a snapshot of a single window. Here, we combine all of those predictions to create a list
of all predictions at the nth position of an input sequence. It makes it easier to query predictions at the beginning or end of the sequence given to SegmentNT. 
'''

import sys
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description='Parse output from run_nt.py to make it easier to find predictions at any position in the input sequence for exons and introns')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input file from run_nt.py output.')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output text file')
    args = parser.parse_args()
    return args

def reformatRunNTResults(args):
    intron_pred=[]
    exon_pred=[] #list of lists with each position in the list being an array of all predictions at that context level
    inArray=False
    intron=False
    exon=False
    pos=0
    lastRow=False
    lastlastRow=False
    
    #Read each line of the run_nt.py output
    with open(args.input) as inputF: ###open(output_from run_nt.py)
        for line in inputF:
            line=line.strip()
            if not inArray:
                if not inArray and not line.startswith("Probabilities,"):
                    continue
            if line.startswith("Probabilities, "):
                line=line[15:] #remove "Probabilities, " from the line
                pos=0
                if line.startswith('intron: [['):
                    inArray=True
                    intron=True
                    line=line[10:] #remove 'intron: [[' from the line
                elif line.startswith('exon: [['):
                    inArray=True
                    line=line[8:] #remove 'exon: [[' from the line
                    exon=True
            elif line.startswith('['):
                line=line[1:]
                pos=0
            if line.endswith("]]"):
                lastlastRow=True
                line=line[0:-2] #remove the bracket from the end of the line
            elif line.endswith("]"):
                line=line[0:-1] #remove the bracket from the end of the line
                lastRow=True
        
            predictions=map(float,line.split()) #get predictions and store as floats. Ensures that all values are floats.
        
            for num in predictions:
                if intron:
                    if len(intron_pred)==pos:
                        intron_pred.append([])
                    intron_pred[pos].append(num)
                elif exon:
                    if len(exon_pred)==pos:
                        exon_pred.append([])
                    exon_pred[pos].append(num)
                pos+=1
            if lastRow:
                pos=0
                lastRow=False
            if lastlastRow:
                inArray=False
                intron=False
                exon=False
                lastlastRow=False
    
    with open(args.output,'w') as output:
        output.write(str(len(intron_pred)) +" " +str(len(exon_pred)) +"\n")
        
        output.write("Intron\n")
        for array in intron_pred:
            output.write(str(array[1:-1]) +"\n") #original program includes one extra nt before and after sequence. Get rid of this
            #output.write(str(array) +"\n") 
        
        output.write("\nExon\n")
        for array in exon_pred:
            output.write(str(array[1:-1])+"\n")#original program includes one extra nt before and after sequence. Get rid of this
            #output.write(str(array)+"\n")
            
def main():
    args=parseArgs()
    reformatRunNTResults(args)

if __name__ == '__main__':
    main()

