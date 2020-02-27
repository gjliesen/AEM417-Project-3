# Gregory Liesen
# AEM 417 Project 3

# Imports
import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import signal
from scipy import linalg as la
import georinex as gr

def read_rinex():
    # This function will read in the satellite data from the text files and return them
    print('Read rinex data')
    rinex_dat = gr.load('brdc2930.11n')
    rinex_df = rinex_dat.to_dataframe()
    print(rinex_df)



def read__novatel():
    print('Read novetel data')


def gen_least_sqr():
    print('Generalized Least Squares Position')


def vdop():
    print('VDOP')


def hdop():
    print('HDOP')

def pdop():
    print('PDOP')

def main():
    read_rinex()
    print('Done')
main()
