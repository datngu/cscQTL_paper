#!/usr/bin/env python

import sys
import glob
import os
import argparse
import gzip
import pandas as pd
import numpy as np



def get_attr(info, key = 'gene_id', pick = 0):
    x = info.split(';')
    idx = np.argmax([i.find(key) for i in x])
    x = x[idx].replace('"', "")
    res = x.split(' ')[-1].split(',')[pick]
    return res


def ciri_quant_convert_bed(path):
    data = pd.read_csv(path, comment ='#', sep = '\t', header = None)
    data[9] = [get_attr(i) for i in data[8]]
    data[10] = [get_attr(i, 'bsj') for i in data[8]]
    order = [0,3,4,9,5,10]
    data = data.iloc[:, order]
    data.columns = ['chr', 's', 'e', 'id', 'cpm', 'strand']
    out = path.replace('.gtf', '.bed')
    data.to_csv(out, sep = '\t', header = False, index = False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Convert ciriquant gtf result to bed format (1-base coordinate format)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--out", help = "output file surfix, replace .gtf by the provided surfix", default = "_ciriquant.bed")
    parser.add_argument("--ciriquant", nargs='+', required = True, help = "list of Circall results")
    args = parser.parse_args()

    # args.ciriquant
    [ciri_quant_convert_bed(fi) for fi in args.ciriquant]
