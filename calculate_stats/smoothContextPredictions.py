import json
import argparse

'''
The purpose of this program is to show that if we average predictions across 24 contexts then the results no longer oscillate.
'''

def parseArgs():
    parser = argparse.ArgumentParser(description='Smooth predictions across sliding window')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required input json file containing a dictionary with key:value pairs of probabilities.')
    parser.add_argument('-o', '--output', type=str, required=True, help='(output) Required output json file containing a dictionary with key:value pairs of smoothed probabilities.')
    parser.add_argument('-w', '--window', type=int, default=24,required=False,help='Window size to smooth predictions. Default: 24 nt')
    args = parser.parse_args()
    return args



def smoothPreds(preds,window):
    smoothedPreds=[]

    for x in range(0,len(preds)-window+1):
        smoothedPreds.append(sum(preds[x:x+window])/window)
    return smoothedPreds


def main():
    args=parseArgs()
    with open(args.input) as inputF, open(args.output,'w') as outputF:
        line1=inputF.readline().strip()
        assert line1.startswith("exon=["),'First line must start with exon=['
        preds= json.loads(line1[5:])
        exons=smoothPreds(preds,args.window)
        outputF.write("exon=" +json.dumps(exons) +"\n")
        line2=inputF.readline().strip()
        assert line2.startswith("intron=["),'First line must start with exon=['
        preds= json.loads(line2[7:])
        introns=smoothPreds(preds,args.window)
        outputF.write("intron=" +json.dumps(introns) +"\n")

        

if __name__ == '__main__':
    main()

