import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from argparse import ArgumentParser
import numpy as np

parser=ArgumentParser()
parser.add_argument("-p","--pfam", help="pfam families for histogram")
parser.add_argument("-b","--bins", help="pfams bp covered in ccr %ile bins")
parser.add_argument("-o","--output", help="histogram pic")
args=parser.parse_args()
pfams=open(args.pfam,"r")
bins=open(args.bins,"r")

for pfam in pfams:
    for bin in bins:
        if pfam == bin:
            plt.savefig(args.output)
