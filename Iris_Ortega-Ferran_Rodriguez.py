import numpy as np 
import math as m
import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = 'Iris_Ortega-Ferran_Rodriguez.py', description = 'Computes')
    parser.add_argument('inpath', help = 'Path of the input file. Example: /home/iris/Desktop/UNI/8è/MQ_N-cossos/Entrega/MQNCSU_Gross-Pitaevskii/')
    arguments = parser.parse_args()

# Directoris
    inpath = arguments.inpath
    data = pd.read_csv(inpath, sep=" ", header=None)
    data.columns = ["a0", "n1", "step", "aa", "time", "alpha", "iter"]
    print(data)