import numpy as np 
import math as m
import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = 'Iris_Ortega-Ferran_Rodriguez.py', description = 'Computes')
    parser.add_argument('inpath', help = 'path')
    arguments = parser.parse_args()

# Directoris
    inpath = arguments.inpath
    with open(inpath) as f:
        for line in f:
            print(line.strip('\n').split(' '))

    data = pd.read_csv(inpath, sep = ' ', names=['AA','N1','step','a0','ALPHA','time','iter'], header=None)
    print(data)