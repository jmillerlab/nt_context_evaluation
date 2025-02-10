#! /usr/bin/env python

'''
The purpose of this program is to run SegmentNT on a gene sequence.
A sliding window of 1 nt is used to determine how hyperparameters impact exon and intron predictions
'''

import os
import re
import math
import sys
import numpy
import haiku as hk
import jax
import jax.numpy as jnp
import argparse
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"]="false"
from nucleotide_transformer.pretrained import get_pretrained_segment_nt_model


def parseArgs():
    parser = argparse.ArgumentParser(description='Run SegmentNT multiple times using a sliding window of 1 nt. Specify a gene sequence and context size to iterate over all possible contexts.')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input fasta file with the gene sequence. Only 1 header and sequence are allowed, although the sequence can be split on multiple lines.')
    parser.add_argument('-o', '--output', type=str, required=True,help='(output) Required output text file to write the predictions.')
    parser.add_argument('-n', '--number_of_tokens_per_seq', type=int, default=4096,required=False,help='Number of tokens per sequence')
    parser.add_argument('-t', '--token_size', type=int, default=6,required=False,help='Token size. Most models currently use a token size of 6. Do not change unless a model uses a different token length')
    parser.add_argument('-m', '--max_seqs_per_run', type=int, default=2,required=False,help='The number of sequences to supply to SegmentNT for each run. Dependant on GPU capabilities')
    parser.add_argument('-p', '--padding', type=int, default=50000,required=False,help='The number of nucleotides before and after the gene sequence that are included in the fasta file')

    args = parser.parse_args()
    return args


def runSegmentNT(args):
    '''
    This function runs SegmentNT multiple times so that each position in a gene sequence will have the change to be the beginning,
    middle, end, and anywhere in between of a sequence passed to SegmentNT. Predictions are made for all possibilities.
    '''
    numpy.set_printoptions(threshold=sys.maxsize)
    max_num_nucleotides=args.number_of_tokens_per_seq
    token_length=args.token_size
    max_sequences_per_run=args.max_seqs_per_run 
    #Use the pre-trained segment_nt model
    parameters, forward_fn, tokenizer, config = get_pretrained_segment_nt_model(
        model_name="segment_nt",
        max_positions=max_num_nucleotides+1,
        verbose=True
    )
    
    forward_fn = hk.transform(forward_fn)
    random_key = jax.random.PRNGKey(0) # Initialize random key
    padding=args.padding
    gene_seq=""
    with open(args.input) as inputF:
        header=inputF.readline()
        seq = inputF.readline().strip()
        while seq !="" and seq[0]!=">":
            gene_seq +=seq
            seq=inputF.readline().strip()
    #Take out the padding so that only the area of interest is analyzed. Default padding is 50000 nt.
    gene_seq=gene_seq[padding-(max_num_nucleotides*token_length)+1:(-1*(padding-(max_num_nucleotides*token_length+1)))] #the +1 is needed so the window can overlap with one nt in the gene seq
    sequences = []
    for x in range(len(gene_seq)-(max_num_nucleotides *token_length)): 
        seq = gene_seq[x:x+(max_num_nucleotides*token_length)]
        seq=seq.upper()
        seq=seq.replace("\n","")
        seq = re.sub(r'[^ATCG]','N',seq)
        sequences.append(seq)
    
    with open(args.output,'w') as output:
        output.write("Number of windows: " + str(len(sequences)) +"\n")
        for x in range(0,len(sequences),max_sequences_per_run):
            output.write(f"Start: {x}\n".format(x))
            #Run SegmentNT in batches. Default 2, specified by user
            tokens_ids = [b[1] for b in tokenizer.batch_tokenize(sequences[x:x+max_sequences_per_run])]
            tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)
            outs = forward_fn.apply(parameters, random_key, tokens)
            logits = outs["logits"]
            probabilities = jnp.asarray(jax.nn.softmax(logits, axis=-1))[...,-1]
            output.write(f"Probabilities shape: {probabilities.shape}\n")
            output.write(f"Features inferred: {config.features}\n")
            # Get probabilities associated with intron
            idx_intron = config.features.index("intron")
            probabilities_intron = probabilities[..., idx_intron]
            output.write(f"Intron probabilities shape: {probabilities_intron.shape}\n\n")
            output.write("Probabilities, intron: " +str(probabilities_intron) + "\n")
            # Get probabilities associated with exon
            idx_exon = config.features.index("exon")
            probabilities_exon = probabilities[..., idx_exon]
            output.write(f"Exon probabilities shape: {probabilities_exon.shape}\n\n")
            output.write("Probabilities, exon: " + str(probabilities_exon) +"\n")
    
def main():
    args=parseArgs()
    runSegmentNT(args)

if __name__ == '__main__':
    main()
